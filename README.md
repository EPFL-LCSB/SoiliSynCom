# SoilSynCom

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Research code, curated genome-scale metabolic models, figures, and source data
supporting the manuscript **“Genome-scale metabolic models for a synthetic
soil microbial community as a path for understanding community functioning.”**

The study uses a species-resolved metabolic framework for a defined 21-member
soil synthetic community. Manually curated genome-scale metabolic models are
used to study growth phenotypes, minimal nutrient requirements, feasible
uptake and secretion profiles, niche overlap, biosynthetic similarity, gene
essentiality, and candidate metabolic exchanges.

![SoilSynCom study workflow](figures/main/Figure_2_Workflow.png)

## Repository contents

| Path | Contents |
|---|---|
| [`models/`](models/) | 80 corrected final models and a checksum manifest |
| [`code/reconstruction/`](code/reconstruction/) | Model reconstruction entry points and project helpers |
| [`code/model_curation/`](code/model_curation/) | Manual curation scripts and the model-repair utility |
| [`code/thermodynamic_conversion/`](code/thermodynamic_conversion/) | FBA-to-TFA conversion scripts |
| [`code/analyses/`](code/analyses/) | iMM, DiMEs/NOI, NIS, and essentiality materials |
| [`plotting/`](plotting/) | Cleaned manuscript notebooks and supplementary plotters |
| [`figures/`](figures/) | Canonical main figures, Supplementary V5, and additional artwork |
| [`tables/`](tables/) | Graphing data and Supplementary V5 source tables |

## Final model collection

The MATLAB `.mat` collection contains four model families:

| Medium | Formulation | Folder | Models |
|---|---|---|---:|
| Rich/PTGY | FBA (short) | `models/final/rich/FBA_short` | 21 |
| Rich/PTGY | TFA (long) | `models/final/rich/TFA_long` | 21 |
| M9 | FBA (short) | `models/final/M9/FBA_short` | 19 |
| M9 | TFA (long) | `models/final/M9/TFA_long` | 19 |

The M9 set excludes `cohm` and `mip`. The `cxy` and `lur` M9 models use xylose
as their carbon source even where retained historical TFA filenames mention
succinate. File-level model type, medium, carbon source, dimensions, and
SHA-256 checksums are recorded in
[`models/model_manifest.tsv`](models/model_manifest.tsv).

The public copies repair annotation-vector alignment and rebuild gene-reaction
associations from the models' `grRules`. Scientific network content outside
those targeted fields is unchanged. Earlier model generations are preserved
separately by the authors and are not part of this repository.

## Code and software

The repository contains selected first-party research scripts used during
model construction, curation, analysis, and figure preparation. It is a clean
record of the code and publication assets supporting the study; it is **not a
turnkey end-to-end reconstruction pipeline**. Several workflows depend on
licensed solvers, external databases, and software versions configured on the
original workstation.

Principal external tools include MATLAB, RAVEN, the COBRA Toolbox, matTFA,
ReMIND, redGEM/redGEMX, NIS, NICEgame, IBM ILOG CPLEX, and Gurobi. Python
plotting notebooks use common scientific packages including NumPy, pandas,
SciPy, Matplotlib, Seaborn, Biopython, ETE Toolkit, and openpyxl. External
packages are not vendored here; see [`NOTICE.md`](NOTICE.md) and the folder
READMEs for context.

## Figures and tables

- `figures/main` contains the eight canonical manuscript PDFs plus lightweight
  PNG previews.
- `figures/supplementary` contains all 31 current Supplementary V5 sheets in
  PDF, editable SVG, and PNG.
- `figures/additional` and `tables/graphing` retain the complete scientific
  graph/table asset collection from the project graphing folder.
- `tables/supplementary` contains the source tables used by Supplementary V5,
  including the Supplementary Table S1 inputs.

The canonical PDFs include final manual artwork assembly. Plotting code and
editable assets document the underlying scientific panels but are not claimed
to reproduce every Illustrator-level layout edit exactly.

## Citation

Please cite the accompanying manuscript. Machine-readable software citation
metadata and the complete author order are provided in
[`CITATION.cff`](CITATION.cff). The article DOI will be added after assignment.

## License and contact

First-party software in this repository is distributed under the
[GNU General Public License v3.0 only](LICENSE). Third-party software and data
remain subject to their respective terms; see [`NOTICE.md`](NOTICE.md).

Requests for resources should be directed to Vassily Hatzimanikatis at
<vassily.hatzimanikatis@epfl.ch>.
