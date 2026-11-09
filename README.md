# OncoBio-AgentFM

Research Software for Oncology AI.

The [multimodal cohort audit](docs/cohort_audit.md) links locked TCIA/NBIA and GEO sources
at patient level, reports observed modality overlap and complete-source selection differences,
and distinguishes pooled selection from R01/AMC recruitment composition. It keeps source
inventory separate from clinical prognostic eligibility. SeriesDate remains series metadata;
actual acquisition/availability, PET/CT component linkage and assay provenance require review. The published
[aggregate evidence](docs/evidence/source_availability.json) contains no patient rows or image
identifiers. Reviewed clinical cases require source lineage before the existing temporal and
overall-survival gates can be applied; source presence alone does not establish validity.

The [endpoint-feasibility assessment](docs/cohort_feasibility.md) defines the
supported research scope: reproducible source-cohort construction, observed
selection differences, and evaluation of prognostic endpoint feasibility.
The locked sources do not establish eligible patient or validated event counts
for either prediction setting. Comparing recorded death durations with death
dates minus a candidate surgery date—calculated from CT date plus the recorded
CT-to-surgery interval—yields 16 disagreements among 63 recorded deaths.
Agreement in the remaining 47 cases does not validate a surgery-based time
origin. These findings do not establish a clinically eligible prognostic cohort.
Unknown eligible-patient and event counts remain `null`.

The [NSCLC research protocol](docs/research_protocol.md) defines all-cause survival
prognosis at a documented pretreatment assessment and a separate exploratory day-30
postoperative landmark. Its [versioned specification](configs/research/nsclc_os_v1.json)
fixes the estimands, hypotheses, permitted evidence and evaluation rules. Patient-level
temporal gates reject post-index evidence, unverified clinical availability and unsupported
historical staging. Overall-survival follow-up is joined after predictor eligibility.
The source bundle requires a reviewed administrative-cutoff policy shared by the cohort
and rejects documented cancer treatment or surgery after death. Original CT/PET roles
remain uniquely declared when frozen evidence is reused at the postoperative landmark.
Metric-specific outcome support, fixed-grid/group prerequisites and paired modality-cohort
feasibility are audited with statistical precision explicitly unestablished, separately
from predictive performance. The [evidence directory](docs/evidence/README.md) distinguishes
synthetic software validation from inspection of the TCIA source schema. Real-data temporal
eligibility remains subject to source adjudication; no model performance is reported.

The implemented NSCLC clinical baseline checks patient identity, evidence timing, event and
censor dates, AJCC stage version, ECOG status, and explicit missingness. It fits a regularized
Cox model on complete training risk sets after patient-level splitting. Clinical source
adjudication is required before constructing these records. The repository includes synthetic
software tests; it reports no cohort performance or clinical validation. Baseline evidence can
be scored without a known outcome; survival follow-up is joined only for model fitting.

## Runtime

- Python 3.13.15
- PyTorch 2.14.0
- MONAI 1.6.0
- NumPy 2.5.3

## Safety Boundary

This repository is research-only software, not a medical device and not for clinical decision-making. Do not commit patient data, credentials, or identifiable clinical artifacts. See `SECURITY.md` before introducing external serialized artifacts.

## License

Copyright (c) 2026 Chandima Liyana. Use is governed by the OncoBio-AgentFM Research Evaluation License 1.0 in `LICENSE`.
