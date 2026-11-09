# Auditable NSCLC multimodal cohort construction

This audit asks whether observed multimodal source completeness selects a different
population from patients with incomplete source evidence. It links patient identities,
exposes unresolved imaging selection and separates repository presence from clinical
prognostic eligibility. It reports no model performance or clinical validation.

The [scope and endpoint-feasibility decision](cohort_feasibility.md) specifies which original
cohort objectives the locked sources can establish. This is a constructed **source cohort**;
construction of an eligible prognostic analysis cohort remains incomplete.

## Observed sources and denominator

The audit inspected the public TCIA NSCLC-Radiogenomics clinical release, an NBIA
collection-wide series response, GEO GSE103584 sample metadata and its processed RNA
expression matrix on **2026-11-07**. The [source lock](../configs/cohort/nsclc_source_lock_v1.json)
records exact compressed/raw byte fingerprints, sizes, retrieval dates and snapshot versions.
The reports also fingerprint all package Python sources, including inherited policy, using filename-framed hashes
with CRLF normalized to LF, so the implementation identity is portable across operating systems.
The clinical CSV checksum is identical to the source inspected in the existing protocol.
The public [aggregate report](evidence/source_availability.json) is generated from these
files; the [flow](evidence/source_flow.md) describes its source denominator.

Exact identifiers are joined without trimming, case folding, zero padding or fuzzy matching.
The source union is the denominator, so unmatched imaging or molecular patients cannot
disappear through an inner join. In this snapshot, all imaging and RNA patients match a
clinical row, yielding 211 unique patients. This verifies exact source-key agreement;
it does not independently authenticate patient identity or discover undocumented aliases.

| Observed source population | Patients | Recorded dead labels |
|---|---:|---:|
| CT metadata + clinical row | 211 | 63 |
| CT + PT metadata + clinical row | 201 | 63 |
| Imaging + linked RNA metadata and expression | 130 | 45 |
| CT + PT + clinical + RNA source evidence | 130 | 45 |
| Same-study PT and CT + clinical row | 199 | 63 |
| Same-study PT and CT + clinical + RNA | 129 | 45 |

These are **source-population counts and raw vital-status labels**, not adjudicated OS
events, day-730 events, clinically eligible participants or selected usable image volumes.
`clinical` means a linked row; it does not mean a complete or temporally valid clinical
baseline. `rna` specifically means the GSE103584 RNA-seq sample/matrix linkage, not any
mutation entry in the clinical CSV and not a clinically accessible molecular result.

The exhaustive CT/PT/clinical/RNA availability matrix has three nonempty patterns:

| CT | PT | Clinical | RNA | Patients |
|---|---|---|---|---:|
| Present | Present | Present | Present | 130 |
| Present | Present | Present | Absent from inventory | 71 |
| Present | Absent from inventory | Present | Absent from inventory | 10 |

All 16 patterns, including zero counts, are retained in the JSON evidence. No patient is
excluded from the inventory because an optional modality is absent.

## Imaging identity and duplicate reconciliation

The NBIA response contains 1,351 canonical series in 395 studies: 727 CT, 480 PT and
144 SEG series. Exact repeated metadata rows with the same SeriesInstanceUID collapse
once. Conflicting records for a series UID, including different owners, modalities or
metadata, stop construction. A StudyInstanceUID cannot belong to two patient identifiers.
There are no exact duplicate rows in the inspected response.

Distinct series UIDs remain distinct candidates. There are 201 patients with multiple CT
candidates and 148 with multiple PT candidates. A first row, largest volume or favorable
outcome does not select an image. Diagnostic chest CT, attenuation-correction CT, localizers,
reconstructions and repeated acquisitions are not distinguished by modality code alone.
Pixel-level duplication, image geometry and volumetric quality are unassessed by this
metadata audit. DICOM SEG objects do not create CT/PT availability.

PT and CT exist in the same StudyInstanceUID for 199 patients, including 129 linked RNA
patients. This is a metadata co-occurrence check; it does not establish attenuation
correction, registration, quantitative PET validity or an appropriate diagnostic CT.
The broader 201-patient PT count must not be labelled a verified 201-patient paired PET/CT
prediction cohort. Conversely, shared study identity is not a clinical eligibility rule:
clinical PET/CT admission requires an explicitly reviewed PT-to-CT component relationship,
not an automatic same-StudyInstanceUID filter. Both components require correct patient/source
identity and actual acquisition/availability by the original pretreatment index. Cross-study
source IDs can be reconciled only through explicit reviewed pairing evidence. Co-occurrence
alone never proves pairing, spatial registration, CT quality or clinical usability.

Twenty-five CT series have no SeriesDate. SeriesDate `(0008,0021)` describes a series,
not verified acquisition timing. For PET it may be an implementation-specific temporal
reference unrelated to acquisition start. Actual acquisition uses separate DICOM acquisition
attributes or independently reviewed clinical source evidence. NBIA getSeries does not
establish those actual acquisition dates, even when SeriesDate is populated.

The inventory preserves SeriesDate as `series_date_on_source_axis` and leaves
`acquired_on_source_axis` null. Image timing provenance labels its raw metadata tag as
`SeriesDate_0008_0021`, retains the original value (including null), and separately cites the
reviewed actual `acquired_on`, acquisition/timing reference and common-calendar alignment
reference. It never equates acquisition with SeriesDate, StudyDate or repository release.
A missing SeriesDate does not invalidate an independently sourced acquisition date. The
reviewed acquisition and clinical availability still must pass the original index gate.
The software verifies assertions and references; it cannot establish their clinical truth.

## Molecular linkage and missingness

The 130 GEO sample titles match 130 unique expression-matrix patient columns. All columns
contain numeric observations. The matrix contains 22,126 gene rows and 1,192,713 `NA`
cells among 2,876,380 cells. `NA` is retained as missing expression; numeric zero remains
an observed value. Duplicate gene identifiers, repeated patient columns, malformed rows,
nonnumeric values and infinities are rejected. Gene-level missingness does not mean that
the patient's entire assay was unacquired.

GEO aliquots sharing an exact patient title remain grouped under that patient. The source
inventory does not fit an imputer, normalize expression, select genes or establish molecular
biological quality. The source publication's tissue-selection and retrospective research
assay context remains material to interpretation.

Absence from these snapshots has clinical reason `unknown`. It is not automatically
`not_acquired`, a negative test, a normal result or lack of clinical indication. Existing
`MissingnessReason` and `EvidenceRecord` contracts govern adjudicated evidence; the
inventory does not replace them. Known `unusable`, `unavailable`, `not_recorded` or
`not_acquired` reasons require supporting local source review. Clinical access to a result
at the prediction time is separate from GEO deposition or specimen collection.

The report also retains recorded field states rather than collapsing all source tokens.
For example, the CSV has 49 PET dates labelled `Not Collected`; 39 of those patients
nevertheless have PT series in NBIA. That field state therefore cannot establish that PET
was never acquired. EGFR has five `Not collected` and 34 `Unknown` labels, with 172 other
entries; these states remain distinct. The `N/A` death dates for 148 living-labelled patients
are not themselves missing OS endpoints. Actual endpoint integrity still requires the existing
follow-up adjudication. Only date-field presence is compared; source calendars are not
assumed numerically comparable.

## Selection hypothesis and observed evidence

The descriptive hypothesis is that complete-source patients differ from clinical-source
patients with at least one source type absent. This is evaluated inside the 211 clinical
records: 130 complete-source patients versus 81 incomplete-source patients. The comparison
uses age at histological diagnosis, source-recorded gender, institution and smoking status.
These variables are descriptive source fields; they are not automatically eligible predictors.
Pathological TNM, recurrence, treatment and outcomes are not used to choose the cohort.

| Descriptive variable | Complete sources, n=130 | Incomplete sources, n=81 | SMD, complete minus incomplete |
|---|---:|---:|---:|
| Age at histological diagnosis, mean (SD), years | 68.51 (9.08) | 67.09 (11.33) | 0.138 |
| Recorded Female, n (%) | 34 (26.15%) | 42 (51.85%) | -0.546 |
| VA source, n (%) | 73 (56.15%) | 20 (24.69%) | 0.677 |
| Recorded nonsmoker, n (%) | 20 (15.38%) | 28 (34.57%) | -0.454 |

The pooled population mixes two source recruitment cohorts. All 130 complete-source patients
are R01 patients; the incomplete group contains 32 R01 and all 49 AMC patients. The report
therefore also gives R01/AMC-stratified descriptions without changing the analysis denominator.
Within R01, recorded Female is 34/130 (26.15%) versus 9/32 (28.13%); VA source is 73/130
(56.15%) versus 20/32 (62.50%). AMC has no complete-source comparison group, so its within-
cohort contrasts are unassessed. Pooled differences must not be interpreted as a gender-
specific acquisition mechanism or as an effect independent of recruitment cohort.

These observed differences support concern that source-complete analysis represents a
different measured population. They do not establish a causal selection mechanism, quantify
unmeasured selection, identify MAR/MNAR, or remove bias. No p-values or universal SMD
threshold are used. Empty groups and zero-variance denominators produce null contrasts,
not infinity or a fabricated zero effect.

Continuous SMD uses the difference in observed group means divided by the square root
of the mean of the two sample variances. A categorical level uses the difference in observed
proportions divided by `sqrt((p1*(1-p1) + p0*(1-p0))/2)`. Category denominators exclude
explicitly missing/unrecognized source values and report those counts separately. Age
denominators and missing counts are also explicit. The public report uses fixed category
vocabularies; unrecognized source text is not exported.

## Clinical eligibility, survival and event budget

The inventory is not a prognostic cohort. The existing pretreatment and day-30 postoperative
gates retain their original definitions. No assessment date, diagnosis availability, clinical
stage edition, ECOG value, endpoint origin, last-known-alive date or clinical molecular
availability is invented to make a source patient pass. The source schema's existing
chronology limitations remain unresolved. Eligible real participant and event counts are
therefore **unassessed (`null`)**, not zero.

The source publication documents chronological ordering among shifted clinical dates.
That is positive evidence about clinical chronology; it does not verify alignment with NBIA
dates, clinical access to predictors, an earliest qualifying assessment or the protocol's
administrative cutoff. The new `endpoint_feasibility` section independently checks recorded
date/status consistency and the diagnostic calculation `CT Date + Days between CT and surgery`.
Of 63 recorded deaths, 47 `Time to Death` values match that candidate surgery-to-death interval
and 16 disagree; one candidate surgery falls after recorded death. Among the 130 four-source
patients, the corresponding counts are 30 agreements and 15 disagreements among 45 recorded
deaths, including the same chronology query. These are source queries, **not** automatic
exclusions or 47/30 validated events. Clinical death/last-contact dates themselves are present
where expected; lack of all outcome data is not the reason the existing estimands are blocked.

The fixed original hash split applied administratively to the 130 complete-source patients
gives 87/21/22 train/validation/test records and 32/7/6 recorded dead labels. These are raw
source feasibility counts, not valid landmark events or evidence of adequate precision.
The split seed, fractions, endpoints, landmark, horizons and metric policies are unchanged.

`audit_linked_cases` accepts existing `AdjudicatedCase` objects only after verifying their
reviewed values, source lineage and complete inventory accounting. It then calls the existing
`audit_eligibility`, `validate_predictors` and follow-up `anchor` functions. Predictor gates
run before survival joins. The returned reports retain separate primary, PET-available,
molecular-available and complete-modality landmark event/partition counts. Complete-modality
status does not move the prediction time, and landmark survivors retain the original
day-30 origin. No fitting, seed search or performance evaluation occurs.

Malformed sources fail as a source-level error; they are not silently dropped as rows.
Preconstruction exclusions require a verified, patient-specific source-review record and
one first-failure reason. Cases plus documented exclusions must cover the entire source
union exactly once. Optional-modality absence alone is not an exclusion. Clinical gate
failures and follow-up failures remain separate and preserve the existing stable reason codes.

## Source lineage and privacy boundary

The [private manifest schema](../schemas/cohort_manifest.schema.json) defines one record
per patient, the stable patient partition, modality presence, unknown source missingness and
source links. Source links retain snapshot hash and original coordinates: zero-based clinical
CSV data-row number, zero-based NBIA JSON-array index, GEO sample accession, or one-based
expression patient-column number after the gene-label column. Image links retain original
source-axis series dates, StudyInstanceUID and SeriesInstanceUID locally. Available-on fields stay
null until clinically adjudicated; this source manifest never authorizes model input.

The clinical review API uses schema **1.1**. Its local JSON fields are exactly:

- `schema_version: "1.1"`;
- `case`: the complete `case_payload(case)` representation;
- `references`: every evidence/timing/adjudication/endpoint reference, plus explicit paired-CT
  timing/pairing references, each with a source SHA-256 and locator;
- `image_date_mappings`: each selected CT/PT and paired-CT source, with `source_axis_date`,
  `source_date_kind: "SeriesDate_0008_0021"`, independently reviewed `acquired_on`,
  `acquisition_reference` and `alignment_reference`;
- `pet_ct_pairings`: one entry per selected PT source, containing `ct_evidence` under the
  existing EvidenceRecord shape, `timing_source`, `pairing_source` and
  `availability_basis: "clinical_record"`;
- `molecular_provenance`: one entry per selected RNA sample, recording the expression-file
  hash and original patient-column locator, resection source ID, specimen-collection date,
  clinical result-availability date/basis and timing reference.

The complete review file is hash-locked. `verify_documents` retains read-only snapshots of
verified supporting bytes rather than return a set of unsubstantiated hash strings.
Clinical references must resolve through object-only RFC 6901 pointers to the exact patient
assertion bundle in an independently hashed **normalized source-adjudication document**:

```json
{
  "schema_version": "1.0",
  "evidence_type": "source_adjudication_assertions",
  "cases": {
    "LOCAL_PATIENT_KEY": {
      "case": {},
      "pet_ct_pairings": {},
      "molecular_provenance": {}
    }
  }
}
```

The empty objects above show structure only; real entries must equal the complete reviewed
case and its explicit pairing/molecular assertions. A clinical reference such as
`/cases/LOCAL_PATIENT_KEY` cannot resolve to another patient, a nonexistent entry or a
scope-only file. An image row or GEO record cannot establish clinical diagnosis/access.
The normalized document is an inspectable source-review assertion ledger, not proof that
an underlying clinical record was authentic or correctly interpreted.

Image references still resolve to the exact NBIA patient/modality/row. Pair components have
separate timing gates. RNA references resolve to the exact sample and processed column;
multiple GEO aliquots for one patient are not automatically assigned to one patient-named
column. Clinical access to a retrospective research assay is never inferred from GEO
existence, deposition, processing date or collection of the specimen. The existing
AvailabilityBasis and specimen/resection/landmark rules remain mandatory.

Sources plus documented exclusions must account for all inventory patients. Source coordinates
are captured when parsing the original bytes and stored independently of mapping iteration
order. Collections and nested mappings are immutable, with unique series/study ownership and
coordinate domains checked. Supporting exclusion records use the previously verified bytes.

Human review must establish source meaning, diagnosis, acquisition, result release, earliest
assessment, source selection, date alignment and cutoff policy. Structured assertion matching
makes these claims inspectable and blocks accidental cross-patient/nonexistent citations;
it cannot authenticate them, resolve undocumented aliases, certify image quality or establish
clinical truth. Low-level preexisting clinical constructors retain their documented trusted-
adjudication assumptions. The new source-validation entry point is `audit_linked_cases`.
Patient-level records, review payloads and normalized assertion ledgers remain private.

Public outputs are deliberately constructed aggregate projections. They contain no patient
IDs, image UIDs, sample accessions, patient dates, source rows or arbitrary source free text.
Aggregate disclosure is not a formal anonymity guarantee. Public-source attribution,
source URLs, dataset-level dates, byte hashes and coarse aggregate counts remain inspectable.
Local source/audit data belong under ignored `data/` and `outputs/`; private review filename
patterns are also ignored. Git ignore rules are assistance, not a substitute for staged-file review.

## Verification and reproducibility

The cohort modules use the Python standard library and the existing research contracts.
They do not change dependency pins, modeling code, clinical schemas or the protocol file.
`cohort_validation` creates invented source files and reviewed cases, exercises linkage,
missingness, duplicate and chronology guards, and produces deterministic expected/observed
rejection evidence. The tests additionally use a separately authored three-patient fixture to challenge source tampering, orphan identities,
RNA NA/zero semantics, private/public projections, selection denominators and attrition.
CI checks the synthetic snapshot on the existing Linux/Windows matrix and uploads synthetic
evidence. It does not fetch real patient data or independently recertify the source snapshot.

The real audit reads four local files against the reviewed lock. A byte or size mismatch
stops it; even a live API ordering change requires source-snapshot review rather than silently
changing the expected hash. Outputs are exclusively created beneath local `outputs/`.
Successful source inventory exits with code 2 because clinical readiness remains blocked.
No public patient-level manifest is generated or automatically promoted.
With `--local-manifest`, a separate `PRIVATE_endpoint_queries.json` traces diagnostic queries
to exact CSV row coordinates and the locked clinical checksum. It contains no inferred clinical
eligibility or survival targets and must remain under ignored local `outputs/`.

## Source discrepancy and interpretation limit

The inspected clinical CSV records 43 Female labels among 162 R01 rows, whereas the original
publication describes 38 females in that cohort. These are source-recorded fields, not an
independently resolved demographic truth. The discrepancy is retained as an unresolved
source-description issue; no labels are recoded or clinical fairness claim inferred from it.

## Source attribution

- Bakr et al. Data for NSCLC Radiogenomics, Version 4. TCIA.
  [Dataset DOI](https://doi.org/10.7937/K9/TCIA.2017.7hs46erv),
  [collection](https://www.cancerimagingarchive.net/collection/nsclc-radiogenomics/),
  [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/).
- Bakr et al. A radiogenomic dataset of non-small cell lung cancer. *Scientific Data* (2018).
  [DOI](https://doi.org/10.1038/sdata.2018.202).
- [NBIA Search REST API guide](https://wiki.cancerimagingarchive.net/display/Public/NBIA+Search+REST+API+Guide).
- [NCBI GEO GSE103584](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103584).

These are research-source and software-validity audits. No medical-device, clinical-grade
certification, prospective validation, general-NSCLC representativeness or clinical utility
claim is made.

- [DICOM PET series timing, PS3.3 C.8.9.1.1.2](https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.8.9.html).
- [DICOM acquisition attributes, PS3.3 C.7.10.1](https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.7.10.html).

The independent standard-library calculator `tools/reproduce_source_audit.py` imports no
project modules. Given the four locked files, it checks raw overlap, matrix cells, recorded
status/partition budgets, endpoint arithmetic and pooled/stratified descriptive calculations against the aggregate
report. `source_reproduction.json` records that comparison. It cannot reproduce clinical
eligibility, valid OS events or image quality from sources that do not establish those facts.
