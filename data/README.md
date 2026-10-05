# SoilSynCom data

We provide the experimental measurements, processed analysis data, graphing
workbooks, and supplementary source tables used in the SoilSynCom study.

## Directory map

- `graphing/` contains the study workbooks, delimited tables, and tree input
  used for the main analyses and manuscript figures.
- `supplementary/` contains the source data for Supplementary V5 and
  Supplementary Table S1.

## Dataset groups

| Path | Data |
|---|---|
| `graphing/ExperimentalGrowthData/` | PTYG and succinate OD600 measurements, BioLog data, growth rates, substrate utilization, and confusion metrics |
| `graphing/ExperimentalValidation/` | Growth-validation measurements and community experiments |
| `graphing/GEMStats/` | Model dimensions and genome-scale model statistics |
| `graphing/iMM/` | In-silico minimal-media solutions, identifiers, classifications, and summary statistics |
| `graphing/DiMEs/` | Compound annotations, DiMEs summaries, uptake/secretion matrices, conditioned-media profiles, and niche-overlap matrices |
| `graphing/NIS/` | Biosynthetic-similarity matrices for amino acids, cofactors, lipids, and nucleotides |
| `graphing/GeneEssentiality/` | Categorized essential-gene results for rich PTYG and M9 conditions |
| `graphing/Curation/` | Model-curation and gap-filling data |
| `supplementary/source_data/` | Figure-level TSV, CSV, JSON, and Newick source data for Supplementary V5 |
| `supplementary/GFUnified.xlsx` | Unified gap-filling source workbook for Supplementary Table S1 |

The plotting scripts reference these paths relative to the repository root.
