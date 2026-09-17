# Documentation

Start with the root [README](../README.md) for setup and execution.

## Use and develop

- [Local development](LOCAL_DEVELOPMENT.md): build and verification commands.
- [Workflow configuration](P5_WORKFLOW.md): local/Slurm execution, budgets and resume.
- [Data contracts](DATA_CONTRACT_V1.md): canonical manifests and output schemas.
- [Alignment backends](P4_BACKENDS.md): tool boundaries and counting policy.
- [Contributing](../CONTRIBUTING.md): code layout, conventions and change checks.

## Validation evidence

- [Implementation status](IMPLEMENTATION_PROGRESS.md): completed gates and limits.
- [Report-data validation](P6_REPORT_DATA.md): public study provenance and reduced comparisons.
- [Python-free runtime](P7_RUNTIME.md): runtime audits, performance and clean staging evidence.
- [Legacy recovery](LEGACY_RECOVERY.md): recover retired code in a separate directory.

## Design and history

- [Maintainability design](CODEBASE_REFACTOR_DESIGN.md) and
  [execution plan](superpowers/plans/2026-09-17-maintainability.md).
- [Historical migration plans](history/README.md).

Phase identifiers in existing guide filenames are retained to keep established
links usable. Historical plans describe decisions at the time; current usage and
validation limits are in the guides above.
