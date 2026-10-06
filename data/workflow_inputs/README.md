# Workflow inputs

We organize the compact scientific inputs used by the SoilSynCom research code
in this directory. `input_manifest.tsv` indexes the complete input collection
across GitHub and the companion Zenodo package, recording each package path,
scientific role, workflows, byte size, packaged and source SHA-256 checksums,
packaging status, and source-archive location. For workbooks,
`metadata_sanitized` means that personal absolute paths were removed from
Excel metadata while the scientific workbook content was preserved.

## Contents

| Folder | Contents |
|---|---|
| `thermodynamics/` | Compact Alberty data, metabolite thermodynamic annotations, reactions treated without thermodynamic constraints, and strain-specific thermodynamic records |
| `curation/` | PTYG and M9 media, shared curation structures, metabolic tasks, BioLog data, and the 21 strain-specific false-task tables |
| `analyses/dimes/` | Conditioned-environment, metabolomics, exchange, and growth-rate inputs for DiMEs and niche-overlap analyses |
| `analyses/imm/` | Candidate nutrient pools and experimental growth-rate inputs for in-silico minimal media |
| `analyses/nis/` | Subsystem mappings, precursor definitions, trade-off settings, and NIS mapping data |
| `reconstruction/` | Reconstruction and genome metadata |

The companion Zenodo package contains `workflow_inputs/reconstruction/FAAs/`,
`workflow_inputs/reconstruction/genomes/`, `workflow_inputs/experimental/`,
the full `thermo_data.mat` database, and the curation reference model. Set
`SOILSYNCOM_ARCHIVE_ROOT` to the extracted package root so MATLAB can resolve
these paths automatically.

The analysis code reuses the final models under `models/final/`, and the
essentiality summaries and figure-ready analysis tables under `data/graphing/`.
MATLAB code can resolve these locations with `code/soilSynComPaths.m`.
