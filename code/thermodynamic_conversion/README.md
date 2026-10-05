# Thermodynamic conversion

We provide the Live Scripts used to convert the curated rich PTYG and M9 FBA
models into matTFA structures, apply the recorded model-specific corrections,
and update gene identifiers. Companion MATLAB functions and the R summary
notebook support the conversion and quality summaries.

The workflows use matTFA, COBRA Toolbox, thermodynamic parameter data, and a
CPLEX-compatible solver. The resulting TFA models are provided under
`models/final/rich/TFA_long/` and `models/final/M9/TFA_long/`.
