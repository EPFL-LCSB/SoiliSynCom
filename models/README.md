# Genome-scale metabolic models

We provide the final genome-scale metabolic model collection for the
21-member SoilSynCom as MATLAB MAT files.

## Model families

| Directory | Models | Formulation | Medium |
|---|---:|---|---|
| `final/rich/FBA_short/` | 21 | FBA | Rich PTYG |
| `final/rich/TFA_long/` | 21 | TFA | Rich PTYG |
| `final/M9/FBA_short/` | 19 | FBA | M9 minimal medium |
| `final/M9/TFA_long/` | 19 | TFA | M9 minimal medium |

“Short” denotes the flux-balance-analysis model and “long” denotes the
thermodynamics-based flux-analysis formulation. In the M9 collection, `cxy`
and `lur` use xylose as their carbon source. The carbon-source field in the
[`model manifest`](model_manifest.tsv) provides the model-level interpretation.

## Model curation and validation

Each model contains metabolite and reaction annotations aligned with its model
entities. The `rules`, `grRules`, and `rxnGeneMat` fields provide consistent
gene-protein-reaction mappings. All 80 MAT files pass annotation-dimension,
GPR-consistency, matrix-dimension, and serialization checks.

We provide the collection metadata and validation results in:

- [`model_manifest.tsv`](model_manifest.tsv), with the strain, medium, carbon
  source, formulation, model dimensions, and SHA-256 checksum;
- [`VALIDATION.md`](VALIDATION.md), with a concise validation summary; and
- [`validation_summary.json`](validation_summary.json), with machine-readable
  validation results.

We provide the associated model-quality and thermodynamic-coverage
measurements in the [`supplementary source data`](../data/supplementary/source_data/).
