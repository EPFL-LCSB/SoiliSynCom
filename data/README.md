# SoilSynCom data

We provide the scientific workflow inputs, experimental measurements,
processed analysis data, graphing workbooks, and supplementary source tables
used in the SoilSynCom study.

## Directory map

- `workflow_inputs/` contains compact thermodynamic, reconstruction, curation,
  and analysis inputs used directly by the research code, together with the
  GitHub/Zenodo input manifest.
- `graphing/` contains the study workbooks, delimited tables, and tree input
  used for the main analyses and manuscript figures.
- `supplementary/` contains the source data for Supplementary V5 and
  Supplementary Table S1.
- `companion_archive_manifest.tsv` and `companion_archive_collections.tsv`
  catalogue the genome and protein sequences, raw experimental workbooks,
  thermodynamic database, reconstruction databases, essentiality inputs, and
  generated analysis collections in the companion Zenodo package.

## Dataset groups

| Path | Data |
|---|---|
| `workflow_inputs/thermodynamics/` | Compact Alberty and compound-level annotations, relaxed-reaction definitions, and strain-specific thermodynamic records |
| `workflow_inputs/curation/` | PTYG and M9 media, shared curation data, metabolic tasks, BioLog inputs, and strain-specific false-task tables |
| `workflow_inputs/analyses/` | Compact input sets for DiMEs, iMM, and NIS |
| `workflow_inputs/reconstruction/` | Reconstruction and genome metadata |
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

The plotting and analysis scripts reference these paths through the repository
root and `SOILSYNCOM_ARCHIVE_ROOT`. `workflow_inputs/input_manifest.tsv`
provides the GitHub or Zenodo package, package path, size, packaged and source
SHA-256 checksums, scientific role, packaging status, and source-archive path
for each input.
