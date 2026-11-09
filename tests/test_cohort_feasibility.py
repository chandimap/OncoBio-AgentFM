"""Challenge unsupported endpoint inference with clinically consequential source states."""

import csv
import hashlib
import io
import json

import pytest

from oncobio_agentfm.cohort_feasibility import endpoint_feasibility, private_endpoint_queries
from oncobio_agentfm.cohort_sources import load_inventory
from oncobio_agentfm.cohort_validation import synthetic_bundle


@pytest.fixture
def report_for(tmp_path):
    """Invent a single-patient release with completely consistent candidate intervals."""

    def report(changes):
        paths, lock_path = synthetic_bundle(tmp_path, count=1)
        reader = csv.DictReader(io.StringIO(paths["clinical"].read_text()))
        rows = list(reader)
        rows[0].update(
            {
                "Survival Status": "Dead",
                "CT Date": "1/1/2020",
                "PET Date": "1/1/2020",
                "Days between CT and surgery": "9",
                "Date of Death": "1/20/2020",
                "Date of Last Known Alive": "1/19/2020",
                "Time to Death (days)": "10",
            }
        )
        rows[0].update(changes)
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        raw = stream.getvalue().encode()
        paths["clinical"].write_bytes(raw)
        lock = json.loads(lock_path.read_bytes())
        lock["sources"]["clinical"].update(
            sha256=hashlib.sha256(raw).hexdigest(), size_bytes=len(raw)
        )
        lock_path.write_text(json.dumps(lock), encoding="utf-8")
        return endpoint_feasibility(load_inventory(paths, lock_path))

    return report


def diagnostics(report):
    return report["source_populations"]["clinical_release"]["arithmetic_diagnostics"]


def test_even_perfect_source_agreement_cannot_create_clinically_eligible_patients(report_for):
    report = report_for({})
    counts = diagnostics(report)
    assert counts["time_to_death_matches_death_minus_candidate_surgery"] == 1
    assert counts["time_to_death_disagrees_with_death_minus_candidate_surgery"] == 0
    assert report["source_populations"]["all_four_source_types"]["clinical_rows"] == 1
    for setting in ("pretreatment", "postoperative_day_30"):
        assert report["analysis_populations"][setting] == {
            "status": "NOT_IDENTIFIABLE_FROM_LOCKED_SOURCES",
            "eligible_patients": None,
            "survival_events": None,
        }
    assert "first_treatment_resection" in report["unresolved_required_facts"]
    assert "OS_cutoff" in report["unresolved_required_facts"]


def test_an_apparently_matching_CT_origin_does_not_trigger_adaptive_origin_selection(report_for):
    report = report_for({"Time to Death (days)": "19"})
    assert diagnostics(report)["time_to_death_matches_death_minus_CT"] == 1
    assert diagnostics(report)["time_to_death_disagrees_with_death_minus_candidate_surgery"] == 1
    assert report["analysis_populations"]["pretreatment"]["survival_events"] is None


def test_missing_death_date_in_living_patient_is_not_an_invalid_censoring_record(report_for):
    report = report_for(
        {"Survival Status": "Alive", "Date of Death": "N/A", "Time to Death (days)": "N/A"}
    )
    counts = diagnostics(report)
    assert counts["dead_without_parseable_death_date"] == 0
    assert counts["candidate_surgery_recorded_alive"] == 1
    assert counts["time_to_death_surgery_comparison_pairs"] == 0
    assert report["analysis_populations"]["postoperative_day_30"]["survival_events"] is None


@pytest.mark.parametrize("token", ["", "Unknown", "Not Collected", "N/A"])
def test_missing_durations_never_become_zero_or_a_reconstructed_surgery(report_for, token):
    report = report_for({"Days between CT and surgery": token})
    states = report["source_populations"]["clinical_release"]["field_states"]
    assert states["Days between CT and surgery"]["source_missing_token"] == 1
    assert states["Days between CT and surgery"]["zero"] == 0
    assert diagnostics(report)["candidate_surgery_records"] == 0


@pytest.mark.parametrize("token", ["0.5", "NaN", "inf", " 9", "9 ", "99999999999"])
def test_invalid_durations_are_reported_without_rounding_or_overflow(report_for, token):
    report = report_for({"Days between CT and surgery": token})
    states = report["source_populations"]["clinical_release"]["field_states"]
    assert states["Days between CT and surgery"]["invalid"] == 1
    assert diagnostics(report)["candidate_surgery_records"] == 0


def test_negative_gap_is_not_repaired_by_absolute_value_or_sign_reversal(report_for):
    report = report_for({"Days between CT and surgery": "-9"})
    states = report["source_populations"]["clinical_release"]["field_states"]
    assert states["Days between CT and surgery"]["negative"] == 1
    assert diagnostics(report)["candidate_surgery_records"] == 0


def test_zero_day_durations_remain_visible_without_epsilon_or_validated_events(report_for):
    report = report_for({"Days between CT and surgery": "0", "Time to Death (days)": "0"})
    states = report["source_populations"]["clinical_release"]["field_states"]
    assert states["Days between CT and surgery"]["zero"] == 1
    assert states["Time to Death (days)"]["zero"] == 1
    assert diagnostics(report)["time_to_death_disagrees_with_death_minus_candidate_surgery"] == 1
    assert report["analysis_populations"]["pretreatment"]["survival_events"] is None


@pytest.mark.parametrize("token", ["2/30/2020", "2020-01-01", "1/1/20", "private-date-token"])
def test_bad_or_unreviewed_date_formats_cannot_supply_a_time_origin(report_for, token):
    report = report_for({"CT Date": token})
    states = report["source_populations"]["clinical_release"]["field_states"]
    assert states["CT Date"]["invalid"] == 1
    assert diagnostics(report)["death_CT_pairs"] == 0
    assert diagnostics(report)["candidate_surgery_records"] == 0
    assert token not in json.dumps(report)


@pytest.mark.parametrize("death", ["1/9/2020", "1/10/2020"])
def test_death_before_or_on_candidate_surgery_is_a_visible_query_not_an_event(report_for, death):
    report = report_for({"Date of Death": death})
    assert diagnostics(report)["candidate_death_on_or_before_surgery"] == 1
    assert report["analysis_populations"]["postoperative_day_30"]["survival_events"] is None


def test_living_contact_before_candidate_surgery_cannot_imply_landmark_survival(report_for):
    report = report_for(
        {
            "Survival Status": "Alive",
            "Date of Death": "N/A",
            "Date of Last Known Alive": "1/5/2020",
            "Time to Death (days)": "N/A",
        }
    )
    assert diagnostics(report)["candidate_alive_contact_on_or_before_surgery"] == 1
    assert report["analysis_populations"]["postoperative_day_30"]["eligible_patients"] is None


@pytest.mark.parametrize(
    ("changes", "field"),
    [
        ({"Survival Status": "Alive"}, "alive_with_death_date"),
        ({"Survival Status": "Alive"}, "alive_with_numeric_time_to_death"),
        ({"Date of Death": "N/A"}, "dead_without_parseable_death_date"),
        ({"Date of Last Known Alive": "1/21/2020"}, "last_alive_after_death"),
        ({"Date of Death": "1/1/2020"}, "death_on_or_before_CT"),
    ],
)
def test_status_and_chronology_contradictions_are_preserved_for_review(report_for, changes, field):
    assert diagnostics(report_for(changes))[field] == 1


def test_unrecognized_status_is_not_inferred_from_a_death_date_or_duration(report_for):
    report = report_for({"Survival Status": "private-status-token"})
    population = report["source_populations"]["clinical_release"]
    assert population["unrecognized_status"] == 1
    assert population["recorded_dead"] == population["recorded_alive"] == 0
    assert diagnostics(report)["candidate_surgery_recorded_deaths"] == 0
    assert "private-status-token" not in json.dumps(report)


def test_large_integer_gap_does_not_overflow_a_calendar_or_silently_clip_dates(report_for):
    report = report_for({"CT Date": "12/31/9999", "Days between CT and surgery": "9999999999"})
    assert diagnostics(report)["candidate_death_on_or_before_surgery"] == 1
    assert report["analysis_populations"]["pretreatment"]["eligible_patients"] is None


def test_report_contains_no_patient_rows_dates_or_derived_survival_targets(report_for):
    report = report_for({})
    encoded = json.dumps(report)
    for private_value in ("SYN-cohort-0", "1/1/2020", "1/20/2020", "GSM900000"):
        assert private_value not in encoded
    for target in ('"event":', '"duration_days":', '"censored_on":', '"index_on":'):
        assert target not in encoded


def test_private_endpoint_queries_keep_source_row_lineage_and_never_become_exclusions(tmp_path):
    inputs = synthetic_bundle(tmp_path, count=2)
    inventory = load_inventory(*inputs)
    report = private_endpoint_queries(inventory)
    assert report["privacy"] == "PRIVATE_DO_NOT_COMMIT"
    assert (
        report["clinical_source_sha256"]
        == hashlib.sha256(inputs[0]["clinical"].read_bytes()).hexdigest()
    )
    assert report["records"] == [
        {
            "patient_id": "SYN-cohort-1",
            "locator": "csv/data-row/1",
            "query_codes": ["dead_without_parseable_death_date"],
        }
    ]
    public = endpoint_feasibility(inventory)
    assert public["source_populations"]["clinical_release"]["clinical_rows"] == 2
    assert public["analysis_populations"]["pretreatment"]["eligible_patients"] is None
