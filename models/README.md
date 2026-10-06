# Genome-scale metabolic models

The `final/` directory contains the original final MATLAB models used in the
SoilSynCom study.

| Directory | Models | Meaning |
|---|---:|---|
| `final/rich/FBA_short/` | 21 | Flux balance analysis models constrained for the rich PTYG medium with vitamins |
| `final/rich/TFA_long/` | 21 | Thermodynamics-based flux analysis versions of the rich PTYG models |
| `final/M9/FBA_short/` | 19 | Flux balance analysis models constrained for M9 minimal medium with vitamins |
| `final/M9/TFA_long/` | 19 | Thermodynamics-based flux analysis versions of the M9 models |

`FBA_short` denotes the conventional stoichiometric model structure used for
flux balance analysis. `TFA_long` denotes the corresponding model extended
with the variables and constraints used for thermodynamics-based flux
analysis.

The M9 collection contains 19 strains because `cohm` and `mip` are not part of
that model set. In the M9 collection, `cxy` and `lur` use xylose as the carbon
source despite the historical `Succinate` wording in their TFA filenames.
