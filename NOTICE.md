# Third-party software and data notice

This repository contains selected first-party SoilSynCom research code. It
does not redistribute the external toolboxes, solvers, or bulk reference
databases used by the original workstation workflows.

## External scientific software

Depending on the workflow, the code refers to:

- [MATLAB](https://www.mathworks.com/products/matlab.html)
- [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)
- [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)
- [matTFA](https://github.com/EPFL-LCSB/matTFA)
- [ReMIND](https://github.com/EPFL-LCSB/remind)
- [redGEM, redGEMX, and lumpGEM](https://github.com/EPFL-LCSB/redgem)
- [NIS](https://github.com/EPFL-LCSB/nis)
- [NICEgame](https://github.com/EPFL-LCSB/NICEgame)
- [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio)
- [Gurobi Optimizer](https://www.gurobi.com/)

Install these projects from their official sources and follow their own
license and citation requirements. Their source trees are not included here.

Python notebooks additionally use packages from the NumPy, SciPy, pandas,
Matplotlib, Seaborn, Biopython, ETE Toolkit, Jupyter, openpyxl, and XlsxWriter
ecosystems. Consult each package for its current license and installation
instructions.

## Databases and identifiers

Models and scripts may contain identifiers derived from KEGG, ModelSEED, NCBI,
and other biological resources. No bulk KEGG or HMM database is distributed in
this repository. Access to and reuse of external database content are governed
by the respective provider's terms.

## Retained historical filenames

Some model filenames are retained so that they remain traceable to the
manuscript analysis. In particular, `cxy` and `lur` M9 TFA filenames may refer
to succinate even though those models use xylose. This is documented metadata,
not a recommendation to infer experimental conditions from filenames alone.

Any third-party copyright or license header already present in an adapted
source file remains controlling for that file. Unless otherwise stated,
first-party SoilSynCom software is licensed under GPL-3.0-only.
