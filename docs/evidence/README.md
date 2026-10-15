# Evidence scope

`temporal_validation.json` is generated from invented cases by
`python -m oncobio_agentfm.protocol_validation`. It includes the protocol's exact
SHA-256 fingerprint, expected/observed rejection codes, aggregate synthetic
counts, training-only horizon screens and metric-specific paired-cohort feasibility checks. It reports no patient performance. Recompute and compare it with:

```bash
python -m oncobio_agentfm.protocol_validation --check-report docs/evidence/temporal_validation.json
```

`source_schema_audit.json` records an inspection on 2026-10-15 of the public
[TCIA clinical CSV](https://www.cancerimagingarchive.net/wp-content/uploads/NSCLCR01Radiogenomic_DATA_LABELS_2018-05-22_1500-shifted.csv).
It contains its checksum, headers, stage vocabulary and aggregate raw-release budget, not the CSV or patient rows.
The raw hash partitions have 145/32/34 records and 45/10/8 dead labels. Those are
administrative release counts, not eligible OS events, landmark counts or 730-day outcomes.
Clinical data attribution: Bakr et al.,
[TCIA dataset DOI](https://doi.org/10.7937/K9/TCIA.2017.7hs46erv),
[CC BY 3.0](https://creativecommons.org/licenses/by/3.0/). Schema metadata were
summarized into this report; row-level dates, IDs and outcomes are excluded.

The CSV inspection is **not** a temporal eligibility validation. Its status remains
`NOT_ESTABLISHED` / `BLOCKED_PENDING_SOURCE_ADJUDICATION`. The ordinary test suite
does not redownload or authenticate TCIA data. Unknown local provenance or changed
headers require renewed source review. Acquisition dates do not prove result availability.

New reports use exclusive creation and cannot overwrite an existing file. Local
source/adjudication reports belong under ignored `outputs/`; aggregate reports are
not automatically safe for public release. Any intentionally revised public snapshot
must be separately reviewed with its protocol version and evidence scope.
CI artifacts contain synthetic evidence only.

The revised synthetic report records 58 temporal checks, 10 metric-feasibility probes
and five source-consistency probes. The latter check cutoff-policy evidence/cohort
consistency and contradictory treatment/surgery after death. Missing or duplicated original
imaging roles cannot bypass the frozen-baseline selection check.
A PASS means the observed guard/status behavior matches the expected synthetic behavior.
It does not mean every metric is estimable. `NOT_ESTIMABLE`, `WITHHELD_BY_PROTOCOL`,
`NOT_ASSESSED` and `OUTCOME_PREREQUISITES_MET` are distinct outcomes. Precision remains
`NOT_ESTABLISHED`; prediction outputs have not been validated or scored. The protocol
documents the revision and its justification. Calibration group support is unassessed
without frozen training/evaluation risk inputs.
