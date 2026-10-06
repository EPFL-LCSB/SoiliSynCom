# SoilSynCom research code

We organize our research code by scientific task:

- `reconstruction/` contains the RAVEN draft-model reconstruction and
  annotation workflows.
- `model_curation/` contains strain-specific, BioLog-guided, and shared model
  curation scripts and notebooks.
- `thermodynamic_conversion/` contains the FBA-to-TFA conversion workflows.
- `analyses/` contains the iMM, DiMEs/NOI, NIS, and essentiality materials.

`soilSynComPaths.m` resolves all repository data and model locations relative
to this checkout. External software locations are read from the environment
variables `MATTFA_ROOT`, `REMIND_ROOT`, `REDGEM_ROOT`, `NIS_ROOT`,
`COBRA_ROOT`, `RAVEN_ROOT`, `CPLEX_ROOT`, and `GUROBI_ROOT`. Generated files
are written under `outputs/` by default, or under `SOILSYNCOM_OUTPUT_ROOT` when
that variable is set.

Set `SOILSYNCOM_ARCHIVE_ROOT` to the extracted companion Zenodo package for
genome and protein sequences, raw experimental data, the full thermodynamic
database, and other large workflow assets.

The compact scientific inputs are organized in
`../data/workflow_inputs/` and indexed by a SHA-256 manifest. The final model
families are available under `../models/final/`.

The root [`NOTICE.md`](../NOTICE.md) lists the scientific software, solvers,
and biological resources used throughout our study.
