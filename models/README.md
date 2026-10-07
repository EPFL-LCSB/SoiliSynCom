# SoilSynCom TFA models

The `final/` directory contains the 80 thermodynamic flux analysis (TFA)
models used in the SoilSynCom study. Models are supplied as MATLAB MAT files,
with one `ttmodel_<strain>` structure in each file.

## Model families

| Directory | Models | Exchange reactions | Meaning |
|---|---:|---:|---|
| `final/curated_exchanges/PTYG/` | 21 | 392 | Models constrained for PTYG medium with the curated exchange repertoire |
| `final/curated_exchanges/M9/` | 19 | 392 | Models constrained for M9 medium with the curated exchange repertoire |
| `final/full_exchanges/PTYG/` | 21 | 1,184 | Models constrained for PTYG medium with the broader exchange repertoire |
| `final/full_exchanges/M9/` | 19 | 1,184 | Models constrained for M9 medium with the broader exchange repertoire |

The curated repertoire is based on *Escherichia coli* exchanges together with
additional study-relevant and metabolomics compounds. The full repertoire
represents the alternative model collection with expanded extracellular
compound coverage. All four collections are TFA models.

PTYG uses glucose as its carbon source and contains all 21 strains. The M9
collections contain 19 strains because `cohm` and `mip` are not represented in
that condition. Seventeen M9 models use succinate; `cxy` and `lur` use xylose.
The curated `cxy` and `lur` filenames retain the historical word `Succinate`,
but their model bounds correctly specify xylose uptake.
