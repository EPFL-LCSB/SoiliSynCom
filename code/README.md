# SoilSynCom research code

This directory contains the selected research code supporting the SoilSynCom
manuscript. It is organized by scientific role rather than as an automated
end-to-end reproduction system:

- `reconstruction/` — RAVEN draft-model reconstruction commands and annotation helper.
- `model_curation/` — strain-specific, BioLog-guided, and shared curation notebooks.
- `thermodynamic_conversion/` — conversion of FBA models to matTFA structures.
- `analyses/` — iMM, DiMEs/niche overlap, NIS, and essentiality documentation.

Scientific calculations in retained source files are unchanged. Only source
selection, duplicate removal, and path portability were addressed. Personal
absolute paths in publication copies are replaced by descriptive tokens such
as `${SOILSYNCOM_DATA_ROOT}`, `${MATTFA_ROOT}`, `${REMIND_ROOT}`,
`${REDGEM_ROOT}`, `${NIS_ROOT}`, `${COBRA_ROOT}`, `${RAVEN_ROOT}`,
`${CPLEX_ROOT}`, `${GUROBI_ROOT}`, and `${USER_HOME}`. Set the corresponding
local locations or adapt the setup cells before using an archival notebook.

Third-party toolboxes and their repositories are not vendored. See the root
`NOTICE.md` for dependency and attribution information. Large data, models,
workspaces, and generated outputs are kept outside this code tree.
