# Thermodynamic conversion

The retained Live Scripts convert the curated rich/PTGY and M9 FBA models to
matTFA structures, apply the recorded model-specific corrections, and update
gene identifiers. Companion MATLAB functions and the R summary notebook are
included where they were part of the current conversion directory.

Large workspaces, converted model binaries, thermodynamic databases, solver
reports, plots, temporary files, and the scratch `test.mlx` notebook are not
duplicated here. The scripts require matTFA, COBRA Toolbox, CPLEX-compatible
solver support, and the thermodynamic input data referenced by the notebooks.
