# NSCLC prognostic research protocol

Protocol `nsclc-os-v1`, schema 1.2, revision `2026-10-15`. The [machine-readable specification](../configs/research/nsclc_os_v1.json)
contains the estimands, hypotheses and evaluation choices. This is repository-versioned
prespecification, not external preregistration. No predictive comparison is reported here.

## Clinical question and population

Among adults with source-confirmed NSCLC in the surgery-referred source cohort,
alive and untreated at a documented pretreatment assessment, does CT evidence
observable at that assessment add prognostic information for all-cause overall
survival beyond the same eligible clinical evidence, under observed care?

The statistical unit is the patient, not a scan, lesion, slice, specimen or assay.
The clinical purpose being studied is prognostic risk characterization, not treatment selection.

The likely source is [TCIA NSCLC-Radiogenomics][1]. The 211-subject collection
contains pretreatment images and tissue-derived molecular data. Its source publication
describes surgery referral and selection, including retrospectively selected molecular-tested
cases [2]. It is not a representative sample of all untreated NSCLC or all surgical candidates.
Selection associated with later resection is already present in the source; an index-date
gate cannot undo it. No primary eligibility rule additionally requires subsequent completion
of surgery. The source's existing surgical selection must accompany any interpretation.

## Prediction time and eligibility

**t0 is the end of the calendar day of the earliest qualifying documented
pretreatment surgical/MDT assessment**, after NSCLC diagnosis confirmation is
available and CT is accessible, before any cancer-directed treatment. A qualifying
assessment must satisfy the sourced baseline criteria below; later optional PET or
molecular availability never moves t0. `index_selection_source` identifies the review
of candidate assessments and the earliest-event decision. This must be an actual sourced
assessment, not the maximum feature date, a convenient CT date, registry export,
postoperative diagnosis date or model-selected date. Fix t0 without consulting outcome values.

Primary eligibility requires:

1. A unique patient identity across evidence and follow-up.
2. Age at t0 of 18–120 years under the existing clinical record contract. Age recorded
   at another clinical event requires explicit adjudication.
3. NSCLC diagnosis demonstrably **available by t0**. Surgical histology cannot retroactively
   establish a known pretreatment diagnosis.
4. Surgical referral and affirmative untreated-at-assessment evidence available by t0.
   Unknown treatment history is not proof of being untreated. The status must be
   verified at t0; an earlier untreated observation cannot establish the intervening history.
   A declared first-treatment date later than documented cancer surgery is contradictory
   and is rejected in both prediction settings.
5. CT acquired and accessible by t0. PET/CT is optional with explicit missingness.
6. A verified common patient time axis for imaging, assessment, treatment and follow-up.
   Independently shifted dates are not directly comparable.
7. For labelled analysis, a source-adjudicated all-cause endpoint with positive follow-up
   from t0. Predictor scoring does not require a known outcome.

Treatment or death on the index calendar day is excluded because within-day ordering
cannot be resolved with this day-level contract. Sources available on t0 are allowed at the
end-of-day decision point. These conventions do not equate acquisition with availability.

Malformed raw records, contradictory vital status and missing mandatory dates require
resolution or exclusion before constructing validated objects. Record those aggregate
source-adjudication exclusions separately. The executable audit reports the first failure
at each gate, not every overlapping reason.

## Permitted evidence and timing

| Evidence | Primary t0 setting | Separate day-30 molecular setting |
|---|---|---|
| Confirmed diagnosis, age at t0 | Sourced and known by t0 | Same frozen baseline information |
| Clinical stage, ECOG | Known by t0, otherwise explicit missingness | Same frozen baseline information |
| CT, optional PET/CT | Acquired and accessible by t0 | Same baseline imaging; later scans cannot substitute |
| Resection-derived molecular result | Prohibited | Specimen linked to first-treatment resection; result accessible by landmark |
| Pathological stage, recurrence, adjuvant treatment, survival fields | Prohibited as predictors | Prohibited in this specified comparison |

`EvidenceRecord` supplies acquisition and availability clocks. For molecular evidence,
acquisition means **specimen collection**, and availability means the usable result was
actually accessible at the clinical landmark. Excision, assay processing, research sequencing,
repository deposition and clinical result release are distinct events. Never assign availability
equal to collection because the release date is absent. Retrospective research assays without
credible clinical availability cannot establish clinical prediction validity.

Retrospective image-feature computation is acceptable only from eligible source images using
transformations fitted without held-out outcomes. Outcome-informed annotations or retrieval
do not become eligible merely because the scan is old. Research segmentations and semantic
annotations require their own source review; this protocol does not certify them as routine
clinical measurements. Image acquisition-to-assessment gaps should be described, not assumed
clinically interchangeable.

Present modalities must carry `AvailabilityBasis.CLINICAL_RECORD`. A result that exists only
in a retrospective research workflow, or whose availability basis is unknown, fails the
clinical prediction gate even if a plausible early date is supplied. `clinical_timing_source`
references review of baseline diagnosis, age and clinical evidence, and `age_reference_date`
must equal t0. These declarations require inspectable source adjudication.

An agent obeys the same evidence whitelist and chronology rules as a statistical model.
It may organize eligible sources; it may not retrieve outcome-bearing records, introduce
postoperative facts at t0 or fabricate missing facts. Adjudication-reference strings identify
the source review required; nonempty strings do not prove clinical truth by themselves.

## Historical staging

The inspected CSV supplies **pathological T/N/M**, without an AJCC edition or a pretreatment
clinical-stage field; see the [source schema report](evidence/source_schema_audit.json).
Raw labels and calendar years do not establish the edition of an individual observation.
Do not relabel historical values as AJCC 8/9, substitute pathological for clinical staging,
or assign a version to satisfy the encoder.

`StageProvenance` retains the raw label, clinical/pathological/unknown basis, original
system and edition, optional adjudicated group and source reference. Only source-supported
**clinical AJCC 8/9** enters `ClinicalStage`, with exact evidence/group/version agreement.
Historical/unknown values stay outside the stage encoder. A retained unusable raw stage is
represented by baseline `MissingnessReason.UNUSABLE`; the original observation stays in
local context, outside eligible model predictors. Quarantine does not claim that it was never
observed. No cross-edition reclassification algorithm is implemented.

The original clinical-stage, cohort and encoder contracts remain strict and unchanged.
A declared cohort encoder edition is not evidence of the source patients' historical edition.
If no eligible clinical stage or ECOG is observed, their missingness indicators are constant;
the effective baseline may be age-only. Describe the actual comparator and its missingness,
and do not call this adjustment for observed stage or performance status. A usable known
clinical stage cannot be deliberately quarantined to alter the comparator.

## Overall survival, time origin and censoring

The event is **death from any cause**, not recurrence, progression or cancer-specific death.
Let D be the death date, A the last date confirmed alive, and C the declared administrative
cutoff on the same verified patient time axis.

- Verified death on or before C is an event at D.
- Verified death after C is censored at C: the later death establishes survival through
  that earlier cutoff.
- For last-confirmed-alive follow-up without an observed death, censor at `min(A, C)`.
  A download or export date does not establish survival through that date.
- Unknown vital status, missing required endpoint dates, last-alive dates after death,
  different time axes and non-OS endpoint definitions are rejected.

Duration is `(endpoint - prediction_date).days`; the event indicator is 1 for death and
0 for right censoring. This protocol requires positive follow-up. The original general-purpose
`SurvivalRecord` continues to permit nonnegative durations. No epsilon or invented day repairs
an invalid record. C must be source-adjudicated, not manufactured from today's date.

`administrative_cutoff_source` must identify a separately reviewed, prespecified cutoff
policy, including any centre-specific cutoffs and their conversion onto verified patient
time axes. A cohort must use one shared policy reference. Equal cutoff dates are not required:
patient date shifts or an explicitly documented centre-specific policy can yield different
calendar values. The shared reference does not authenticate arbitrary dates; that remains
source adjudication. Missing or mixed policy references block analysis. This policy must be
fixed without inspecting model predictions or selecting censoring to improve results.

Before analysis, source-bundle construction rejects a documented cancer-treatment date or
cancer surgery later than verified death. These are contradictory source records to resolve
or count separately, not evidence for a new predictive feature. Predictor validation remains
an outcome-independent function of the original clinical context.

The CSV `Time to Death (days)` field is not reused as a model duration: its original origin
and consistency with the chosen assessment require adjudication. Recurrence or subsequent
treatment does not censor OS. The target is survival under observed subsequent care, not a
counterfactual treatment regime.

Right-censoring methods require adequate independent-censoring assumptions. A marginal
training Kaplan–Meier censoring estimator also assumes censoring is sufficiently independent
of prognosis for that weighting method. Loss to follow-up may violate this assumption.
Timing checks do not establish it; report follow-up patterns and limit interpretation.

## Separate postoperative estimand

**L = documented first-treatment resection date + 30 calendar days.** This is a fixed
design choice for an exploratory conditional question, not an established optimal clinical
landmark. Do not move L to accommodate late assays or improve performance.

The context must satisfy the baseline requirements, with resection verified as the first
cancer treatment and follow-up beyond L. Death or censoring on/before L excludes only this
setting, not an otherwise valid primary record. Clinical and imaging evidence stay frozen at t0.
Present molecular evidence needs a specimen-collection date matching resection and a result
accessible by L. A later scan cannot replace baseline CT at the molecular landmark. The adjudicated case
bundle also requires the same chosen CT/PET source evidence, dates, presence and clinical
access basis in both feature lists. A different pre-index scan or newly introduced PET
record cannot silently replace that frozen baseline. This consistency rule does not
redefine the timing gate or allow a later clinical update.
The original CT and PET roles must each be uniquely declared, including explicit missingness;
an absent or duplicated original role cannot bypass the frozen-selection check. Callers using
the low-level predictor function alone must supply the frozen source selection. The
`AdjudicatedCase.features_for()` path enforces cross-setting source consistency; the low-level
function can enforce dates and allowed roles but cannot compare against a source list it was
not given.

Molecular absence does not exclude a survivor from the landmark cohort flow. H4's added-modality
contrast uses the **same molecular-available survivors in both arms**, while reporting missing
and excluded survivors separately. The survival clock restarts at L; t0-to-L survival is not
counted as predictive event-free person-time. The conditional population differs from the
primary population. Higher landmark performance cannot establish a pretreatment improvement.
Landmarking does not identify a causal effect or remove all selection bias [4].

If molecular availability cannot be verified, only explicitly labelled molecular associations
with independently valid outcomes/time origins are defensible; these do not validate clinical
prediction at either decision point.

## Estimands and hypotheses

The primary estimand is
`P(T - t0 > u | T > t0, primary eligibility, X known at t0)` under observed care.
The secondary estimand replaces t0 with L and conditions on eligible landmark survivors and
evidence known by L. Eligibility is fixed before model fitting or evaluation.

| ID | Setting and paired comparison | Interpretation |
|---|---|---|
| H1 | Primary: clinical + CT versus the same eligible clinical baseline | Δ Uno concordance at 730 days; null Δ ≤ 0, directional alternative Δ > 0 |
| H2 | Exploratory: clinical + CT + PET/CT versus clinical + CT on the same PET-available patients | Increment conditional on PET availability; report selection and coverage |
| H3 | Exploratory: complete available imaging versus masks hiding PET/CT, or both CT and PET/CT, retaining clinical evidence | Paired concordance/Brier changes; imposed absence differs from natural missingness |
| H4 | Exploratory: clinical + CT + molecular versus clinical + CT on the same molecular-available day-30 survivors | Conditional postoperative increment; no pretreatment validity claim |

H3 applies deterministic input ablations **after freezing source eligibility**. It does not
relabel an acquired CT as clinically unacquired or change the cohort. The PET-only ablation
uses the paired CT+PET-available subset. No mask rate is searched. Unsupported masked predictions
count as withheld predictions in coverage, rather than disappearing from the analysis.

H1 is the single primary hypothesis; H2–H4 are exploratory. Report the paired effect and interval,
including null or negative results. A positive point estimate alone is not evidence of improvement;
even a positive interval does not establish clinical benefit. No clinically meaningful threshold is
asserted. Neither Cox coefficients, molecular associations nor agent explanations identify a
causal treatment effect or biological mechanism.

## Evaluation and uncertainty

Reuse `split_patient_ids`: seed `oncobio-agentfm-v1`, nominal 70/15/15 train/validation/test.
All scans, lesions, assays and landmarks of one patient retain the same partition. Hash assignment
is stable when an eligible subset changes; fractions do not guarantee exact counts. Do not search
split seeds for favorable event counts or performance.

Fitting, imputation, scaling, feature selection, harmonization, baseline-hazard estimation and
censoring estimation use training patients only. Tuning and checkpoint selection use validation
patients only. Freeze choices before the one held-out test analysis. Comparison arms share patients,
outcomes, partitions and times. Reuse the existing full-risk-set regularized Efron Cox baseline.
Category/age support and convergence failures must be reported rather than silently scored.

| Measure | Role and required input |
|---|---|
| Uno IPCW concordance truncated at 730 days | Primary discrimination; larger is better; larger risk predicts earlier death [5] |
| Harrell concordance | Descriptive discrimination with censoring-dependent interpretation |
| Cumulative/dynamic AUC at 365, 730, 1095 days | Discrimination at supported horizons |
| IPCW Brier score at those horizons | Probability error; smaller is better; requires survival probabilities |
| Integrated Brier score on a fixed daily grid, days 1–730 | Probability error over the supported interval |
| Censoring-aware calibration at supported horizons | Observed versus predicted survival; report sparse support |
| Coverage, exclusions, deaths/censoring, at-risk counts | Applicability and feasibility, not substitutes for performance |

Raw Cox log-risk is not a survival probability. Probability metrics require a training-estimated
baseline hazard and valid survival function. Concordance truncated at 730 days evaluates ordering of eligible event times and is not
the same estimand as a binary two-year mortality classifier. Event contributions require observed death time strictly less than
tau; deaths exactly on tau do not contribute as index events. Equal-duration outcome
pairs are not ordered; later-observed patients remain eligible comparison partners. Risk ties
receive half credit under the declared metric definition. Calibration uses mean predicted
survival and Kaplan–Meier observed survival in fixed risk thirds with boundaries fitted on
training predictions only; empty/tied groups and sparse confidence intervals are reported.
No pooled calibration statistic is substituted for those curves.

The present baseline exposes log-risk; this protocol
does not assert implementation or evaluation of probability metrics. No test-outcome-selected
high/low-risk cutoff or naive binary label for early-censored patients is permitted.

`assess_horizon_support` is retained as the earlier conservative common screen. It is
not authorization for every metric. `assess_metric_support` governs metric-specific
outcome prerequisites; the two reports state their different scopes explicitly.
Training-only reverse Kaplan–Meier estimates G. At tied training times, deaths leave the
censoring risk set before its censoring jump. The prespecified G >= 0.10 and at-risk floor
10 remain numerical stability policies, not statistical sample-size calculations.

| Metric | Distinct outcome support and fixed design |
|---|---|
| Uno at 730 days | Observed event time < 730 and at least one strictly later observed comparison time; report distinct event patients and comparable pairs separately. More pairs do not create more independent patients. |
| Harrell, untruncated | At least one ordered event/later-time pair; no IPCW G requirement is imposed on this descriptive metric. Equal observed times are excluded under the declared day-level tie rule. |
| Cumulative/dynamic AUC | Cases are observed deaths <= u, including deaths on u. Controls require observed follow-up > u. Both must be nonempty. Early-censored patients are not controls. |
| Brier at a fixed horizon | At least one known event/survivor contribution and supported training censoring weights. It does not require a death on every early day. Survival probabilities remain required and are not validated by this outcome audit. |
| Integrated Brier | Keep every day 1–730, step 1. Check known outcome contributions throughout; monotonic G and risk-set tail bounds cover their full grid. An unsupported tail withholds the fixed integral; it does not justify shortening it. A complete aligned probability grid is independently required. |
| Calibration thirds | Require explicit training/evaluation risk inputs for grouping. Use training-only inclusive linear tertiles; equality at a cutpoint enters the upper group. Tied cutpoints or empty evaluation groups fail the grouping design. Apply the unchanged at-risk floor separately in each group; pooled support is insufficient. |

Calibration observed survival uses within-group Kaplan–Meier; independent censoring within
those risk groups remains an assumption, not something established by a pooled G check.
The per-group floor is a conservative reporting restriction in this amendment, not a
universal recommendation that ten patients ensure accurate calibration. No group merging,
test-derived boundaries or alternative calibration estimator is selected to rescue sparse results.
When grouping risks are absent, calibration group feasibility is `NOT_ASSESSED`.

Metric reports distinguish `NOT_ESTIMABLE` (undefined required denominators/support),
`WITHHELD_BY_PROTOCOL` (a stability policy fails), `NOT_ASSESSED` (required grouping
inputs absent), and `OUTCOME_PREREQUISITES_MET`. The last status does not certify an
actual score: fitted-model support/convergence, finite aligned predictions, probability
validity and complete prediction coverage remain independently necessary. Every report
retains `precision_status: NOT_ESTABLISHED`, `prediction_outputs_validated: false` and
`model_evaluation: NOT_RUN`. No real predictive metric is calculated by this milestone.

H1 uses all source-eligible primary patients. H2 and the PET-only H3 mask use the
prespecified PET-available cohort. H4 uses molecular-available landmark survivors.
The audit applies the original hash split within each declared comparison cohort;
both model arms use the same training, validation and test patients, including the same
training cohort for censoring estimation. Stable hash assignment preserves each patient's
partition across subsets and clinical settings. Separate cohort flow retains missing-modality
patients for coverage. Do not silently intersect successfully scored test patients or replace
withheld predictions by a more favorable evaluation population.

Use 2,000 paired patient bootstrap resamples, seed 20261014, and percentile 95% intervals with linear quantile interpolation for
held-out differences, holding fitted models and the training censoring estimator fixed. Sample the
same patients in both arms. Report undefined replicates and withhold an interval if more than 5%
are undefined. This measures conditional evaluation-sample uncertainty, not training, full population
or external-generalization uncertainty. No family of unadjusted confirmatory claims across H2–H4
is permitted. Report effects, uncertainty and feasibility rather than ranking by p-value.

## Feasibility budget and documented amendment

The inspected full-release CSV has 211 records, 63 `Dead` labels and 148 `Alive` labels.
Those labels are not independently adjudicated OS endpoints. The nominal 15% allocation
is 31.65 people and 9.45 dead labels across all available follow-up, not at 730 days.
The original, unchanged hash split applied administratively to that release gives:

| Raw release partition | Records | Recorded dead | Recorded alive |
|---|---:|---:|---:|
| Train | 145 | 45 | 100 |
| Validation | 32 | 10 | 22 |
| Test | 34 | 8 | 26 |

These are **full-release administrative counts**, not the clinically eligible cohort,
landmark survivors, 730-day event counts or observed performance. Clinical and modality
restrictions can make the informative samples smaller. No calculation establishes adequate
precision for H1, AUC, integrated Brier, calibration or molecular increments. Reliable
precision assessment would need the actual eligible follow-up and prediction distribution;
those are not established. A fixed dataset and numerical gate are not a sample-size justification.
Precision-oriented evaluation methodology cautions against treating a rule of thumb as
sample-size adequacy [7].

Revision 2026-10-15 is a documented methodological amendment to the 19-file implementation
identified by archive SHA-256 `b3dc28a437f64d541e6663d251d42c297bfa21c545820e6ffc4c08c5dc4aa61a`.
The repository base is `a559a99b6a68b7a767d93d7923f4e7dc033e8cc9`.

| Finding | Amendment and justification |
|---|---|
| Frozen CT/PET date gate already correct | Preserve its implementation; add explicit post-assessment/pre-surgery timing probes and clinical-contract regression tests. The case bundle also freezes the chosen CT/PET evidence and presence across settings. |
| Shared horizon screen insufficient for all metrics | Add distinct pair/case/control/grid/group checks and separate undefined, withheld and unassessed states. Retain split, seed, horizons and G/at-risk constants. |
| Pooled counts cannot support paired modality claims | Audit each prespecified modality cohort with identical patients in both arms and every partition. |
| Small raw event budget | Preserve the primary estimation target; make lack of precision/power evidence explicit. Primary means the selected target, not an adequately powered confirmatory claim. |
| Clinical timing contradictions | Require untreated verification at t0 and reject first-treatment dates after documented cancer surgery. |
| Specification integrity | Fingerprint schema 1.1 and recursively reject wrong nested policy types, including bool substituted for integer grid steps. |

The final schema-1.2 refinement starts from the recovered 19-file archive SHA-256
`868891f0108f4a25303a5ad65d9783de0d71c3a91e19889eb1812e8e8f41f4cf`.
It preserves that implementation and adds unique original imaging-role checks, a shared
source-adjudicated administrative-cutoff policy, and rejection of treatment/surgery after
death at source-bundle construction. The chronology policy is explicit and strictly checked
in the JSON specification. No split, horizon, metric definition, threshold, model, dependency
pin or earlier source module is changed. Synthetic evidence now includes 58 temporal checks,
10 metric-feasibility probes and five source-consistency probes. Historical schema-1.1
amendments above describe the preceding revision rather than new work in this refinement.

Administrative source counts, including fixed-partition label totals, were inspected for
feasibility before any predictive model comparison. No model scores, tuning results or
performance-based choices informed this amendment. No test-population, seed, endpoint,
horizon or model-selection substitution was performed. This is versioned amendment disclosure,
not a claim of external preregistration or of never having inspected aggregate outcome counts.
Unsupported comparisons are withheld transparently; H1 remains the prespecified primary
estimation target and all numerical results would remain exploratory without precision evidence.

## Reporting alignment and evidence limits

This protocol supports TRIPOD+AI reporting principles [3, 8]; it is not a completed
prediction-model evaluation or a certification of full checklist compliance.

| Reporting topic | Present evidence and limit |
|---|---|
| Intended purpose, population, care-pathway point | Explicit prognostic question, surgery-selected population, t0/L and research-only purpose. |
| Data and participant eligibility | Named source, original selection, source schema/budget and explicit missing chronology. |
| Outcome and predictors | All-cause death, cutoff/censoring, separate time origins and acquisition/access/source-basis gates. |
| Sample size and partitioning | Fixed available-data budget and patient hash split; adequate precision is explicitly unestablished. |
| Missing data and analysis | Explicit source missingness, paired modality subsets, train-only fitting and unchanged baseline safeguards. |
| Evaluation and uncertainty | Fixed metrics, denominator/group/grid prerequisites, conditional bootstrap rules and unsupported-result disclosure. |
| Reproducibility and results | Versioned specification, tests and synthetic report; no real model result, external validation, fairness result or deployment claim. |

Source-reference strings attest to a review; the software cannot establish their authenticity,
that a complete assessment history was supplied, or real clinical availability. Date consistency
is conditional on a correctly adjudicated common patient calendar. Patient/source IDs,
review references and timing metadata are for audit/linkage; they are not permitted predictor
values. Agentic retrieval follows the same frozen clinical information whitelist.

## Inspectable evidence and acceptance boundary

- [Source schema audit](evidence/source_schema_audit.json): inspected CSV checksum, headers and stage
  vocabulary and aggregate raw release counts; no patient rows, identifiers or dates.
- [Temporal validation](evidence/temporal_validation.json): deterministic synthetic probes,
  expected/observed rejection codes and separate setting audits.
- Executable horizon-support reports distinguish feasible metric inputs from unsupported
  horizons; passing support gates does not establish performance, precision, power or censoring validity.
- Tests challenge forbidden/post-index evidence, historical staging, identity, censoring,
  mismatched time axes and landmark time origins.
- GitHub Actions checks the existing canonical runtime, lint/format, complete tests and the
  committed synthetic snapshot, and uploads a synthetic report on Linux and Windows.

A synthetic pass demonstrates the specified software guard behavior. The actual CSV inspection
is `BLOCKED_PENDING_SOURCE_ADJUDICATION`: mandatory assessment, diagnosis/feature availability,
clinical staging provenance and verified date alignment have not been established. A header audit
is not clinical temporal validation. Real patient adjudication/audits remain local and require
reviewed sources; neither predictive performance nor successful molecular prediction is claimed.
The software is research-only, with no screening, diagnosis, treatment-benefit, deployment,
general-NSCLC, prospective or external-validation claim.

## References

[1]: https://www.cancerimagingarchive.net/collection/nsclc-radiogenomics/
[2]: https://doi.org/10.1038/sdata.2018.202
[4]: https://doi.org/10.1111/tri.13081
[5]: https://doi.org/10.1002/sim.4154

1. Bakr et al. Data for NSCLC Radiogenomics, TCIA Version 4. [Dataset DOI](https://doi.org/10.7937/K9/TCIA.2017.7hs46erv).
2. Bakr et al. A radiogenomic dataset of non-small cell lung cancer. *Scientific Data*, 2018. [DOI][2]; [author-hosted paper](https://web.stanford.edu/group/rubinlab/pubs/Bakr-2018-RadiogenomicDataset.pdf).
3. Collins et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models. *BMJ*, 2024. [DOI](https://doi.org/10.1136/bmj-2023-078378).
4. Gleiss et al. An unjustified benefit: immortal time bias in the analysis of time-dependent events. *Transplant International*, 2018. [DOI][4].
5. Uno et al. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Statistics in Medicine*, 2011. [DOI][5].
6. [scikit-survival: Evaluating survival models](https://scikit-survival.readthedocs.io/en/stable/user_guide/evaluating-survival-models.html). Metric/censoring definitions; no new dependency is introduced here.

[7]: https://doi.org/10.1002/sim.9275
[8]: https://www.tripod-statement.org/wp-content/uploads/2024/04/TRIPODAI-Supplement.pdf

7. Riley et al. Minimum sample size calculations for external validation of a clinical prediction model with a time-to-event outcome. *Statistics in Medicine*, 2022. [Author manuscript](https://pure-oai.bham.ac.uk/ws/portalfiles/portal/187642222/Statistics_in_Medicine_2021_Riley_Minimum_sample_size_calculations_for_external_validation_of_a_clinical_prediction.pdf).
8. TRIPOD+AI expanded checklist, 2024. [Official supplement][8].
