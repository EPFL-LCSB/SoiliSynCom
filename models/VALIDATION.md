# Model package validation

We validated all 80 publication models under `models/final/` for structure,
annotation alignment, gene-reaction consistency, and MAT-file serialization.

| Family | Models | Status |
|---|---:|---|
| Rich PTYG FBA | 21 | Passed |
| Rich PTYG TFA | 21 | Passed |
| M9 FBA | 19 | Passed |
| M9 TFA | 19 | Passed |

Our validation checks confirm that:

- All ten metabolite annotation arrays match `mets`.
- KEGG, Rhea, MIRIAM, and SBO reaction annotation arrays match `rxns`.
- `rxnGeneMat` has shape reactions-by-genes and has no orphan columns.
- Every numeric `rules` index exactly matches its `rxnGeneMat` row.
- Every gene referenced by `grRules` is represented in the GPR matrices.
- The biochemical-network fields retain their values throughout model
  packaging.
- Every model reloads successfully from its MAT file.

KEGG metabolite snapshot SHA-256: `636aea94e2475e5819edaad781f244db3cc265c7472c12fd66cb46a124529809`

KEGG reaction snapshot SHA-256: `932f24189bf777e8b3a5b483c055811013c1b5be3ce3bc963f47771b3200469b`
