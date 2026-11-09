# NSCLC source-cohort and endpoint feasibility

The source-cohort and endpoint-feasibility analysis constructs a reproducible
multimodal **source cohort** and evaluates whether the locked source evidence
supports the existing prognostic protocol. Its scope encompasses source
availability, observed selection differences, and limitations in endpoint
ascertainment. Clinically eligible pretreatment and postoperative prognostic
populations remain unestablished.

## Evidence and observed endpoint queries

The four [locked sources](../configs/cohort/nsclc_source_lock_v1.json) link 211 clinical
patients, CT metadata for 211, PT metadata for 201 and linked RNA-seq for 130. All four source
types overlap in 130 patients. These are repository-evidence counts, not clinical eligibility.

The data descriptor [Bakr et al., 2018](https://doi.org/10.1038/sdata.2018.202), Methods / Subject
Demographics and Clinical Data, documents shifted clinical dates with preserved relative
chronological ordering. Its RNA Sequencing section describes tissue-availability/quality
selection. Neither statement supplies individual pretreatment assessments or clinical RNA
result-release dates. Collection-level preoperative imaging and resection descriptions cannot
establish each patient's diagnosis confirmation, untreated status or first-treatment chronology.

The [generated aggregate report](evidence/source_availability.json), `endpoint_feasibility`,
checks the actual source values. `CT Date`, last-known-alive date and CT-to-surgery interval
parse for all 211 patients. The 63 dead-labelled rows have parseable death dates and integer
`Time to Death (days)`; the 148 alive-labelled rows have no recorded death date or numeric
death duration. No recorded last-known-alive date is after recorded death. These useful
consistency checks do not define an eligible analysis population or a censoring policy.

For a deliberately labelled **diagnostic hypothesis**, let S* = recorded CT date + recorded
CT-to-surgery interval. S* is not an adjudicated first-treatment surgery date. We compare
recorded `Time to Death` with recorded death date minus S*. This does not select a time origin.

| Arithmetic check | All clinical source rows | Four-source overlap |
|---|---:|---:|
| Source patients | 211 | 130 |
| Recorded dead labels | 63 | 45 |
| Death-duration comparison pairs | 63 | 45 |
| Agreement with death minus S* | 47 | 30 |
| Disagreement with death minus S* | 16 | 15 |
| S* after recorded death | 1 | 1 |

None of the 63 recorded death durations equals recorded death date minus CT date. That
rejects exact CT-origin equivalence in these fields; it does not identify the true origin.
The inconsistencies may involve the interval, a date, duration or their source definitions;
the public fields do not identify which is wrong. There is no sign reversal, majority-origin
assignment, date repair, duration replacement or removal of inconvenient records. The
matching 47 are **not** a validated event cohort. Raw event labels are not interchangeable
with adjudicated outcomes. These findings are reproducible arithmetic queries, not an
allegation that all TCIA follow-up is invalid.

`tools/reproduce_source_audit.py` independently recalculates these counts from the four
locked files without importing project modules. With `--local-manifest`, the cohort CLI
also writes local `PRIVATE_endpoint_queries.json`: checksum-bound CSV row references for
the 16 queried patients, including the chronology query. Keep that artifact private.

## Original objective assessment

| Objective | Completed evidence | Unverified requirement for original completion |
|---|---|---|
| Real reproducible multimodal cohort | Locked source union, exhaustive observed availability matrix, private manifest and reproducible aggregates | A real **eligible prognostic** cohort is not constructed. Source presence must satisfy patient-specific clinical criteria. |
| Patient and modality linkage | Exact CSV/NBIA/GEO/matrix keys; study/series ownership; immutable source coordinates; duplicate/conflict gates | Patient identity beyond deposited keys, diagnostic image selection, actual acquisition and explicit usable PET/CT components require source review. |
| Survival and censoring integrity | Date/status syntax and chronology diagnostics; inherited endpoint/cutoff gates; visible duration inconsistencies | A justified prediction origin, resolved endpoint queries, clinical follow-up adjudication and the protocol's reviewed cutoff policy are absent. |
| Prediction-time availability | Existing timing contracts and reviewed-source linkage reject unverified or late evidence | Actual earliest assessment, diagnosis/baseline access, untreated status, image result access and patient-level molecular release are not in the locked bundle. |
| Inclusion/exclusion accounting | Complete source denominator; 130 complete versus 81 incomplete; no missing-modality inner-join deletion; reviewed-exclusion interface | The source flow is not a clinical eligibility flow. Real clinical reasons and resulting clinical attrition remain unassessed. |
| Eligible patients and outcome events | Raw label/partition budgets are explicitly separate from validated outcomes | Both setting-specific eligible and event counts are `null`; neither zero nor 211/63 or 130/45 is justified. Positive landmark follow-up and events after the correct origin/cutoff must be established. |
| Cohort selection differences | Pooled and R01/AMC-stratified source comparisons, effect sizes and explicit missing denominators | These describe source selection, not selection into a clinically eligible population, a causal mechanism, or correction of bias. |
| Provenance, privacy and leakage | Byte locks, private lineage, aggregate-only outputs, fixed patient splits and adversarial source/timing tests | Clinical truth, undisclosed aliases, pixel duplicates, quantitative PET/image quality and cross-source time-axis alignment remain unverified. Native Windows execution is not established by local Linux checks. |

## Why the clinical counts are not identifiable

Identical public fields can be consistent with different qualifying histories. A patient
may have confirmed diagnosis and accessible CT at an untreated MDT assessment, or may have
received induction treatment before such an assessment. The CSV does not distinguish these
histories. A result in GEO can likewise arise from a clinical assay accessible on day 20 or
research processing after day 30. Both histories leave the deposited source presence unchanged
while changing eligibility under the existing protocol. No deterministic loader can recover
the missing facts from those fields.

The required evidence is patient-specific documentation of the earliest qualifying t0,
diagnosis and baseline release, affirmative untreated status, selected image timing/access,
first-treatment resection and specimen linkage, molecular clinical release, source-axis
alignment, endpoint-query resolution and an administrative-cutoff policy. A normalized
assertion ledger makes supplied reviews inspectable; it does not create that evidence.

The cutoff is a requirement of the current protocol, **not a universal requirement that all
retrospective survival studies use one common calendar date**. A separately justified design
could censor at recorded last contact. That would change the endpoint-construction policy;
it would not recover the missing pretreatment t0 or timed clinical molecular availability.
Similarly, untimed research RNA can support a research association question, but cannot be
presented as a clinically available day-30 predictor.

## Supported research claim

The source-complete population differs on measured characteristics, partly reflecting R01/AMC
recruitment composition. For example, pooled recorded Female proportions are 34/130 versus
42/81; within R01 they are 34/130 versus 9/32. All AMC patients lack the linked RNA source.
That structural composition prevents a within-AMC complete/incomplete comparison and does
not establish a gender-specific missingness mechanism. Descriptive SMDs support concern about
selection; they do not quantify the resulting model-performance bias or establish MAR/MNAR.

The supported scientific contribution is **auditable NSCLC source-cohort
construction, observed selection differences, and a negative feasibility finding
for the locked prognostic settings**. Clinically eligible populations and
validated outcome-event counts remain unestablished. The existing prognostic
protocol and patient identity, survival, and evidence contracts remain unchanged
and applicable. No survival model, predictive benefit, clinical validation,
general-NSCLC representativeness, or clinical-grade certification is claimed.