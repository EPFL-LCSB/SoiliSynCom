# Model package validation

All 80 publication models passed structural and serialization validation.
The untouched source archive was read only; corrected copies are stored under
`models/final/`.

| Family | Models | Unused genes removed | GPR genes reconnected | Reaction annotation rebuilds |
|---|---:|---:|---:|---:|
| rich_fba | 21 | 977 | 5 | 6 |
| rich_tfa | 21 | 977 | 5 | 6 |
| m9_fba | 19 | 900 | 4 | 5 |
| m9_tfa | 19 | 900 | 4 | 5 |

Validation checks:

- All ten repaired metabolite annotation arrays match `mets`.
- KEGG, Rhea, MIRIAM, and SBO reaction annotation arrays match `rxns`.
- `rxnGeneMat` has shape reactions-by-genes and has no orphan columns.
- Every numeric `rules` index exactly matches its `rxnGeneMat` row.
- Every field outside the explicitly repaired metadata/GPR fields is unchanged.
- Every output reloads successfully from its MAT file.
- The 21 rich/TFA models are semantically identical to the earlier validated Python-fixed set.
- Five distinct rich-medium and four distinct M9 GPR genes were reconnected in both FBA and TFA formulations.

KEGG metabolite snapshot SHA-256: `636aea94e2475e5819edaad781f244db3cc265c7472c12fd66cb46a124529809`

KEGG reaction snapshot SHA-256: `932f24189bf777e8b3a5b483c055811013c1b5be3ce3bc963f47771b3200469b`
