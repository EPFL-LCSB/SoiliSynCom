# Model curation

This directory retains the current, non-legacy MATLAB and R scripts used for
strain-specific, phenotype-guided, and shared model curation. MATLAB Live
Scripts are kept because many manual curation decisions were recorded in that
format. Models, workspaces, reports, solver logs, temporary files, MEMOTE
reruns, and `Old` trees are intentionally excluded.

- `strain_specific/` contains the staged curation notebooks for each final strain.
- `biolog/` contains the BioLog confusion-matrix and follow-up curation stages.
- `shared/` contains annotation, exchange, balance, and metabolic-task utilities.
- `shared_templates/` stores byte-identical scripts once; `duplicate_files.tsv`
  maps every omitted source alias to the retained copy.
- `repair_model_families.py` (when present) is the separate Python repair used
  to align annotations and GPR matrices in the publication model families.

| Original workflow number | Strain ID | Organism label |
|---:|---|---|
| `1` | `cuvk` | Curtobacterium |
| `2` | `cxy` | Cellulomonas |
| `3` | `bbet` | Bradyrhizobium |
| `4` | `lur` | Luteibacter |
| `5` | `pzu` | Phenylobacterium |
| `6` | `pmor` | Pseudomonas moraviensis |
| `7` | `dri` | Devosia |
| `8` | `cohm` | Cohnella |
| `9` | `tro` | Tardiphaga |
| `10` | `rfa` | Rhodococcus |
| `11` | `pscf` | Pseudomonas sp |
| `12` | `cpi` | Chitinophaga |
| `13` | `vpe` | Variovorax |
| `14` | `flk` | Flavobacterium |
| `16` | `cauj` | Caulobacter |
| `17` | `lgu` | Lysobacter |
| `18` | `mgk` | Mucilaginibacter |
| `19` | `racd` | Rahnella |
| `20` | `mlot` | Mesorhizobium |
| `21` | `buol` | Burkholderia |
| `22` | `mip` | Microbacterium |

These notebooks document the research curation process; they are not a
turnkey sequential workflow. They depend on MATLAB, RAVEN, COBRA Toolbox,
matTFA, compatible optimization solvers, and the external gap-filling tools
not redistributed here.
