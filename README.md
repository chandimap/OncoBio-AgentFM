# OncoBio-AgentFM

Research Software for Oncology AI.

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
