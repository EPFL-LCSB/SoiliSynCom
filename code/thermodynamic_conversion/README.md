# Thermodynamic conversion

We provide the Live Scripts used to convert the curated rich PTYG and M9 FBA
models into matTFA structures, apply the recorded model-specific corrections,
and update gene identifiers. Companion MATLAB functions and the R summary
notebook support the conversion and quality summaries.

Compact Alberty data, compound-level thermodynamic annotations,
relaxed-reaction definitions, and strain-specific thermodynamic records are
provided under
[`data/workflow_inputs/thermodynamics/`](../../data/workflow_inputs/thermodynamics/).
The full `thermo_data.mat` database is stored under
`workflow_inputs/thermodynamics/` in the companion Zenodo package and is
resolved through `SOILSYNCOM_ARCHIVE_ROOT`.
The resulting TFA models are organized by medium and exchange repertoire under
`models/final/curated_exchanges/` and `models/final/full_exchanges/`.

The workflows use matTFA, COBRA Toolbox, and a CPLEX-compatible solver.
