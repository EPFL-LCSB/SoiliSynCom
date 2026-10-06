# Model curation

We provide the MATLAB and R materials used for strain-specific,
phenotype-guided, and shared model curation. MATLAB Live Scripts preserve the
ordered curation steps and their scientific annotations.

- `strain_specific/` contains the staged curation notebooks for the final
  strain models.
- `biolog/` contains the BioLog confusion-matrix calculations and follow-up
  curation stages.
- `shared/` contains annotation, exchange, reaction-balance, and metabolic-task
  utilities.
- `shared_templates/` contains common script templates used across strain
  workflows.
- `repair_model_families.py` aligns model annotations and gene-reaction
  matrices across the four publication model families.
- `model_families.tsv` defines the rich PTYG and M9 FBA/TFA model collections.

The medium definitions, shared curation structures, metabolic tasks, BioLog
inputs, and strain-specific false-task tables are provided under
[`data/workflow_inputs/curation/`](../../data/workflow_inputs/curation/). The
associated experimental measurements are under
`workflow_inputs/experimental/` in the companion Zenodo package. The curation
reference model is stored under `workflow_inputs/curation/shared/` in the same
package.

The curation workflows use MATLAB, RAVEN, COBRA Toolbox, matTFA, and compatible
optimization solvers. The publication models produced from these workflows are
provided under `models/final/`.
