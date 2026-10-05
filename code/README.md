# SoilSynCom research code

We organize our research code by scientific task:

- `reconstruction/` contains the RAVEN draft-model reconstruction and
  annotation workflows.
- `model_curation/` contains strain-specific, BioLog-guided, and shared model
  curation scripts and notebooks.
- `thermodynamic_conversion/` contains the FBA-to-TFA conversion workflows.
- `analyses/` contains the iMM, DiMEs/NOI, NIS, and essentiality materials.

We use project-root variables such as `${SOILSYNCOM_DATA_ROOT}`,
`${MATTFA_ROOT}`, `${REMIND_ROOT}`, `${REDGEM_ROOT}`, `${NIS_ROOT}`,
`${COBRA_ROOT}`, `${RAVEN_ROOT}`, `${CPLEX_ROOT}`, and `${GUROBI_ROOT}` to
identify the software and database locations used by each workflow.

The root [`NOTICE.md`](../NOTICE.md) lists the scientific software, solvers,
and biological resources used throughout our study.
