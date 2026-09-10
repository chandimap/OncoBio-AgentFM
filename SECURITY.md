# Security

OncoBio-AgentFM is research-only software.

No patient data, credentials, private keys, model checkpoints, or identifiable medical images belong in this repository.

Treat externally supplied serialized artifacts and configuration files as untrusted. The code uses MONAI only for in-memory array transforms. It does not load MONAI bundles, executable configuration expressions, pickle files, model checkpoints, or file-backed NumPy objects. 

Report suspected security issues privately to the repository owner rather than publishing exploit details in an issue.
