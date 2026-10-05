# Genome-scale metabolic models

This directory contains the publication model set for the 21-member SoilSynCom.
The models are provided as MATLAB MAT files and retain the variable names used by
the original MATLAB workflows.

## Final model families

| Directory | Models | Formulation | Medium |
|---|---:|---|---|
| `final/rich/FBA_short/` | 21 | FBA | Rich PTYG |
| `final/rich/TFA_long/` | 21 | TFA | Rich PTYG |
| `final/M9/FBA_short/` | 19 | FBA | M9 minimal medium |
| `final/M9/TFA_long/` | 19 | TFA | M9 minimal medium |

`cohm` and `mip` are not included in the M9 families because those strains were
excluded from the final M9 model set. In M9, `cxy` and `lur` use xylose rather
than succinate. Some historical TFA filenames still contain `Succinate`; the
carbon-source column in `model_manifest.tsv` records the correct interpretation.

“Short” denotes the flux-balance-analysis (FBA) model and “long” denotes the
thermodynamics-based flux-analysis (TFA) formulation.

## Publication curation

Corrected copies were generated from the untouched model archive. The curation:

- realigned ten metabolite-annotation vectors to `mets`;
- realigned KEGG, Rhea, MIRIAM, and SBO reaction annotations when needed;
- rebuilt `rules` and `rxnGeneMat` from `grRules`;
- retained every gene referenced by a GPR and removed only genes unused by all
  `grRules` expressions; and
- preserved all fields outside the explicitly repaired annotation and GPR
  fields.

Five distinct rich-medium and four distinct M9 GPR genes that were present in
`grRules` but disconnected in the old `rxnGeneMat` were restored. The same genes
occur in the corresponding FBA and TFA formulations.

`model_manifest.tsv` records the strain identity, formulation, source and donor
lineage, dimensions, SHA-256 checksums, and per-model repair statistics.
`removed_unused_genes.tsv` records every removed unused gene. Validation results
are available in both `VALIDATION.md` and `validation_summary.json`.

The 21 rich/TFA models were checked for semantic identity against the previously
validated Python-fixed collection. All 80 MAT files were reloaded and checked
for annotation lengths, GPR/matrix agreement, absence of orphan gene columns,
and preservation of non-target fields.

## Historical models

Earlier and donor generations are preserved outside this Git repository in the
companion `SoilSynCom_Legacy_Models` archive. They are retained for historical
reference and should not be substituted for the final models here. The legacy
archive also identifies an incomplete historical `pzuV8.mat` file explicitly.

Historical MEMOTE reports remain the manuscript's quality-control record. No
new MEMOTE score is claimed for these annotation-only corrected copies.
