"""Locked source inventory; repository presence never establishes clinical availability."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import date
from math import isfinite
from pathlib import Path
from types import MappingProxyType

from oncobio_agentfm.prognosis import ProtocolViolation
from oncobio_agentfm.source_audit import REVIEWED_SOURCE_COLUMNS

SOURCE_URLS = {
    "clinical": "https://www.cancerimagingarchive.net/wp-content/uploads/"
    "NSCLCR01Radiogenomic_DATA_LABELS_2018-05-22_1500-shifted.csv",
    "imaging": "https://services.cancerimagingarchive.net/nbia-api/services/v1/"
    "getSeries?Collection=NSCLC%20Radiogenomics&format=json",
    "rna_metadata": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/"
    "GSE103584/soft/GSE103584_family.soft.gz",
    "rna_expression": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/"
    "GSE103584/suppl/GSE103584_R01_NSCLC_RNAseq.txt.gz",
}
IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}\Z")
UID = re.compile(r"[0-9]+(?:\.[0-9]+)+\Z")


def implementation_fingerprint() -> str:
    """Fingerprint every package Python source, including inherited clinical policy."""

    digest = hashlib.sha256()
    for path in sorted(Path(__file__).parent.glob("*.py")):
        raw = path.read_bytes().replace(b"\r\n", b"\n")
        digest.update(path.name.encode("ascii") + b"\x00" + hashlib.sha256(raw).digest())
    return digest.hexdigest()


def require(condition: bool, code: str) -> None:
    """Use fixed diagnostics rather than interpolate local patient/source values."""

    if not condition:
        raise ProtocolViolation(code)


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject ambiguous JSON, including duplicate keys in nested source objects."""

    result: dict[str, object] = {}
    for key, value in pairs:
        require(key not in result, "source_duplicate_json_key")
        result[key] = value
    return result


def identifier(value: object) -> str:
    """Keep exact source identity; do not trim, case-fold, pad or fuzzy-match it."""

    require(
        isinstance(value, str) and IDENTIFIER.fullmatch(value) is not None,
        "source_patient_identifier_invalid",
    )
    return value


@dataclass(frozen=True, slots=True)
class SourceReceipt:
    """Public source identity and exact bytes; a checksum is not authentication."""

    kind: str
    version: str
    retrieved_on: str
    sha256: str
    size_bytes: int


def locked_bytes(path: Path, kind: str, spec: dict[str, object]) -> tuple[bytes, SourceReceipt]:
    """Require reviewed byte identity before interpreting any source record."""

    require(
        type(spec) is dict and set(spec) == {"version", "retrieved_on", "sha256", "size_bytes"},
        "source_lock_shape_invalid",
    )
    require(kind in SOURCE_URLS, "source_kind_invalid")
    require(
        type(spec["version"]) is str
        and re.fullmatch(r"[A-Za-z0-9._-]{1,80}", spec["version"]) is not None,
        "source_version_invalid",
    )
    require(type(spec["retrieved_on"]) is str, "source_retrieval_date_invalid")
    try:
        parsed = date.fromisoformat(spec["retrieved_on"])
    except ValueError as exc:
        raise ProtocolViolation("source_retrieval_date_invalid") from exc
    require(parsed.isoformat() == spec["retrieved_on"], "source_retrieval_date_invalid")
    require(
        type(spec["sha256"]) is str and re.fullmatch(r"[0-9a-f]{64}", spec["sha256"]) is not None,
        "source_hash_invalid",
    )
    require(type(spec["size_bytes"]) is int and spec["size_bytes"] > 0, "source_size_invalid")
    raw = path.read_bytes()
    require(len(raw) == spec["size_bytes"], "source_size_mismatch")
    require(hashlib.sha256(raw).hexdigest() == spec["sha256"], "source_hash_mismatch")
    return raw, SourceReceipt(kind, **spec)


@dataclass(frozen=True, slots=True)
class Series:
    """One canonical NBIA series; dates remain on the unverified source calendar."""

    patient_id: str
    study_uid: str
    series_uid: str
    modality: str
    image_count: int
    source_row: int
    series_date: str | None
    study_date: str | None

    def __post_init__(self) -> None:
        identifier(self.patient_id)
        require(
            all(
                type(value) is str and len(value) <= 64 and UID.fullmatch(value)
                for value in (self.study_uid, self.series_uid)
            ),
            "imaging_uid_invalid",
        )
        require(self.modality in {"CT", "PT", "SEG"}, "imaging_modality_unreviewed")
        require(type(self.image_count) is int and self.image_count > 0, "imaging_count_invalid")
        require(type(self.source_row) is int and self.source_row >= 0, "imaging_source_row_invalid")
        for value in (self.series_date, self.study_date):
            if value is not None:
                require(type(value) is str, "source_imaging_date_invalid")
                try:
                    valid = date.fromisoformat(value).isoformat() == value
                except ValueError:
                    valid = False
                require(valid, "source_imaging_date_invalid")


def _source_date(value: object) -> str | None:
    """Preserve unknown dates; never infer a clinical result-release date."""

    if value in (None, ""):
        return None
    require(
        type(value) is str and re.fullmatch(r"\d{4}-\d{2}-\d{2} 00:00:00\.0", value) is not None,
        "source_imaging_date_invalid",
    )
    try:
        return date.fromisoformat(value[:10]).isoformat()
    except ValueError as exc:
        raise ProtocolViolation("source_imaging_date_invalid") from exc


def read_clinical(raw: bytes) -> dict[str, dict[str, str]]:
    """Retain all patients and every raw column without redefining baseline variables."""

    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), strict=True)
    require(tuple(reader.fieldnames or ()) == REVIEWED_SOURCE_COLUMNS, "clinical_schema_changed")
    records = {}
    for row in reader:
        require(
            None not in row and all(value is not None for value in row.values()),
            "clinical_row_width_invalid",
        )
        patient = identifier(row["Case ID"])
        require(patient not in records, "clinical_duplicate_patient")
        records[patient] = row
    require(bool(records), "clinical_source_empty")
    return records


def read_series(raw: bytes) -> tuple[tuple[Series, ...], int]:
    """Collapse byte-equivalent metadata rows; quarantine conflicting UID ownership."""

    rows = json.loads(raw, object_pairs_hook=unique_object)
    require(type(rows) is list and bool(rows), "imaging_source_invalid")
    canonical: dict[str, tuple[dict[str, object], Series]] = {}
    study_owners: dict[str, str] = {}
    duplicates = 0
    required = {
        "PatientID",
        "StudyInstanceUID",
        "SeriesInstanceUID",
        "Modality",
        "ImageCount",
        "Collection",
    }
    for offset, row in enumerate(rows):
        require(type(row) is dict and required <= row.keys(), "imaging_schema_changed")
        require(row["Collection"] == "NSCLC Radiogenomics", "imaging_collection_mismatch")
        patient = identifier(row["PatientID"])
        for key in ("StudyInstanceUID", "SeriesInstanceUID"):
            require(
                type(row[key]) is str
                and len(row[key]) <= 64
                and UID.fullmatch(row[key]) is not None,
                "imaging_uid_invalid",
            )
        study, uid = row["StudyInstanceUID"], row["SeriesInstanceUID"]
        require(
            type(row["Modality"]) is str and row["Modality"] in {"CT", "PT", "SEG"},
            "imaging_modality_unreviewed",
        )
        require(
            type(row["ImageCount"]) is int
            or (
                type(row["ImageCount"]) is str
                and re.fullmatch(r"[0-9]+", row["ImageCount"]) is not None
            ),
            "imaging_count_invalid",
        )
        count = int(row["ImageCount"])
        require(count > 0, "imaging_empty_series")
        require(
            study not in study_owners or study_owners[study] == patient,
            "imaging_study_owner_conflict",
        )
        study_owners[study] = patient
        item = Series(
            patient,
            study,
            uid,
            row["Modality"],
            count,
            offset,
            _source_date(row.get("SeriesDate")),
            _source_date(row.get("StudyDate")),
        )
        if uid in canonical:
            require(canonical[uid][0] == row, "imaging_series_metadata_conflict")
            duplicates += 1
        else:
            canonical[uid] = (row, item)
    return tuple(
        sorted((pair[1] for pair in canonical.values()), key=lambda item: item.series_uid)
    ), duplicates


def read_rna_metadata(raw: bytes) -> dict[str, str]:
    """Link GEO sample titles exactly to TCIA identities; keep aliquots patient-grouped."""

    lines = gzip.decompress(raw).decode("utf-8").splitlines()
    require(lines.count("!Series_geo_accession = GSE103584") == 1, "rna_series_mismatch")
    samples: dict[str, dict[str, str]] = {}
    current = None
    fields = {"!Sample_title", "!Sample_library_strategy", "!Sample_type"}
    for line in lines:
        if line.startswith("^"):
            current = None
        if line.startswith("^SAMPLE = "):
            current = line.removeprefix("^SAMPLE = ")
            require(re.fullmatch(r"GSM[0-9]+", current) is not None, "rna_sample_invalid")
            require(current not in samples, "rna_duplicate_sample")
            samples[current] = {}
        key, separator, value = line.partition(" = ")
        if current is not None and separator and key in fields:
            require(key not in samples[current], "rna_duplicate_sample_field")
            samples[current][key] = value
    result = {}
    for accession, sample in samples.items():
        require(set(sample) == fields, "rna_sample_fields_missing")
        require(
            sample["!Sample_library_strategy"] == "RNA-Seq" and sample["!Sample_type"] == "SRA",
            "rna_assay_type_mismatch",
        )
        result[accession] = identifier(sample["!Sample_title"])
    require(bool(result), "rna_metadata_empty")
    return result


@dataclass(frozen=True, slots=True)
class SourceInventory:
    """Private patient-level inventory; only aggregate projections may be published."""

    receipts: tuple[SourceReceipt, ...]
    clinical: Mapping[str, Mapping[str, str]]
    series: tuple[Series, ...]
    rna_samples: Mapping[str, str]
    rna_expression_patients: frozenset[str]
    expression_rows: int
    expression_missing_cells: int
    expression_observed_per_patient: Mapping[str, int]
    duplicate_series_rows: int
    clinical_source_rows: Mapping[str, int]
    expression_source_columns: Mapping[str, int]

    def __post_init__(self) -> None:
        require(
            type(self.receipts) is tuple
            and all(isinstance(receipt, SourceReceipt) for receipt in self.receipts)
            and len(self.receipts) == len(SOURCE_URLS)
            and {receipt.kind for receipt in self.receipts} == set(SOURCE_URLS),
            "inventory_receipts_incomplete",
        )
        require(
            type(self.series) is tuple and all(isinstance(item, Series) for item in self.series),
            "inventory_series_must_be_immutable",
        )
        require(type(self.rna_expression_patients) is frozenset, "inventory_rna_patients_invalid")
        require(
            type(self.duplicate_series_rows) is int and self.duplicate_series_rows >= 0,
            "inventory_duplicate_count_invalid",
        )
        require(bool(self.clinical), "clinical_source_empty")
        for patient, row in self.clinical.items():
            identifier(patient)
            require(
                set(row) == set(REVIEWED_SOURCE_COLUMNS)
                and row["Case ID"] == patient
                and all(type(value) is str for value in row.values()),
                "inventory_clinical_row_invalid",
            )
        require(
            len({item.series_uid for item in self.series}) == len(self.series),
            "inventory_duplicate_series",
        )
        owners = {}
        for item in self.series:
            require(
                item.study_uid not in owners or owners[item.study_uid] == item.patient_id,
                "imaging_study_owner_conflict",
            )
            owners[item.study_uid] = item.patient_id
        for sample, patient in self.rna_samples.items():
            require(
                type(sample) is str and re.fullmatch(r"GSM[0-9]+", sample) is not None,
                "rna_sample_invalid",
            )
            identifier(patient)
        for patient in self.rna_expression_patients:
            identifier(patient)
        require(
            set(self.rna_samples.values()) == self.rna_expression_patients,
            "rna_matrix_metadata_patient_set_mismatch",
        )
        for coordinates, expected, minimum in (
            (self.clinical_source_rows, set(self.clinical), 0),
            (self.expression_source_columns, self.rna_expression_patients, 1),
        ):
            require(
                set(coordinates) == expected
                and all(type(value) is int and value >= minimum for value in coordinates.values())
                and len(set(coordinates.values())) == len(coordinates),
                "inventory_source_coordinates_invalid",
            )
            require(
                set(coordinates.values()) == set(range(minimum, minimum + len(coordinates))),
                "inventory_source_coordinates_invalid",
            )
        require(
            len({item.source_row for item in self.series}) == len(self.series)
            and all(
                item.source_row < len(self.series) + self.duplicate_series_rows
                for item in self.series
            ),
            "inventory_series_coordinates_invalid",
        )
        require(
            type(self.expression_rows) is int and self.expression_rows > 0,
            "inventory_expression_rows_invalid",
        )
        require(
            all(
                type(count) is int and 0 <= count <= self.expression_rows
                for count in self.expression_observed_per_patient.values()
            ),
            "inventory_expression_counts_invalid",
        )
        require(
            set(self.expression_observed_per_patient) == self.rna_expression_patients,
            "inventory_expression_patient_set_mismatch",
        )
        require(
            type(self.expression_missing_cells) is int
            and self.expression_missing_cells
            == self.expression_rows * len(self.rna_expression_patients)
            - sum(self.expression_observed_per_patient.values()),
            "inventory_expression_missingness_mismatch",
        )
        object.__setattr__(
            self,
            "clinical",
            MappingProxyType(
                {patient: MappingProxyType(dict(row)) for patient, row in self.clinical.items()}
            ),
        )
        object.__setattr__(self, "rna_samples", MappingProxyType(dict(self.rna_samples)))
        object.__setattr__(
            self, "clinical_source_rows", MappingProxyType(dict(self.clinical_source_rows))
        )
        object.__setattr__(
            self,
            "expression_source_columns",
            MappingProxyType(dict(self.expression_source_columns)),
        )
        object.__setattr__(
            self,
            "expression_observed_per_patient",
            MappingProxyType(dict(self.expression_observed_per_patient)),
        )

    @property
    def patient_ids(self) -> frozenset[str]:
        """Use the union so failed joins cannot silently disappear from the denominator."""

        return (
            frozenset(self.clinical)
            | {item.patient_id for item in self.series}
            | frozenset(self.rna_samples.values())
            | self.rna_expression_patients
        )


def read_expression(raw: bytes) -> tuple[dict[str, int], int, int]:
    """Audit every matrix cell; NA is missing expression, never numeric zero."""

    reader = csv.reader(
        io.StringIO(gzip.decompress(raw).decode("utf-8")), delimiter="\t", strict=True
    )
    header = next(reader, [])
    require(len(header) > 1 and header[0] == "", "rna_matrix_header_invalid")
    patients = [identifier(value) for value in header[1:]]
    require(len(set(patients)) == len(patients), "rna_duplicate_patient_column")
    observed = dict.fromkeys(patients, 0)
    genes: set[str] = set()
    missing = 0
    for row in reader:
        require(len(row) == len(header), "rna_matrix_row_width_invalid")
        require(bool(row[0]) and row[0] not in genes, "rna_gene_identifier_invalid_or_duplicate")
        genes.add(row[0])
        for patient, value in zip(patients, row[1:], strict=True):
            if value == "NA":
                missing += 1
                continue
            try:
                numeric = float(value)
            except ValueError as exc:
                raise ProtocolViolation("rna_expression_nonnumeric") from exc
            require(isfinite(numeric), "rna_expression_nonfinite")
            observed[patient] += 1
    require(bool(genes), "rna_matrix_empty")
    return observed, len(genes), missing


def load_inventory(paths: dict[str, Path], lock_path: Path) -> SourceInventory:
    """Read only local, version-locked inputs; no network or identity repair occurs."""

    lock = json.loads(lock_path.read_bytes(), object_pairs_hook=unique_object)
    require(
        type(lock) is dict and set(lock) == {"schema_version", "sources"}, "source_lock_invalid"
    )
    require(
        lock["schema_version"] == "1.0" and type(lock["sources"]) is dict,
        "source_lock_version_invalid",
    )
    require(set(paths) == set(lock["sources"]) == set(SOURCE_URLS), "source_set_incomplete")
    source = {
        kind: locked_bytes(paths[kind], kind, lock["sources"][kind]) for kind in sorted(SOURCE_URLS)
    }
    clinical = read_clinical(source["clinical"][0])
    series, duplicates = read_series(source["imaging"][0])
    samples = read_rna_metadata(source["rna_metadata"][0])
    observed, rows, missing = read_expression(source["rna_expression"][0])
    require(set(observed) == set(samples.values()), "rna_matrix_metadata_patient_set_mismatch")
    return SourceInventory(
        tuple(item[1] for item in source.values()),
        clinical,
        series,
        samples,
        frozenset(observed),
        rows,
        missing,
        observed,
        duplicates,
        {patient: index for index, patient in enumerate(clinical)},
        {patient: index + 1 for index, patient in enumerate(observed)},
    )
