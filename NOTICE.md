# Software and data sources

We used the following scientific software and optimization tools in the
SoilSynCom workflows:

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

Our Python notebooks use packages from the NumPy, SciPy, pandas, Matplotlib,
Seaborn, Biopython, ETE Toolkit, Jupyter, openpyxl, and XlsxWriter ecosystems.
Each software package retains its respective copyright, license, and citation
terms.

The models and data use biological identifiers from resources including KEGG,
ModelSEED, NCBI, and RAVEN. Access and reuse of these resources follow the
terms established by their providers.

## Model filename interpretation

The `cxy` and `lur` M9 models use xylose as their carbon source. Their exact
medium and carbon-source assignments are recorded in
`models/model_manifest.tsv`.

Any third-party copyright or license header present in an adapted source file
remains controlling for that file. We license first-party SoilSynCom software
under GPL-3.0-only.
