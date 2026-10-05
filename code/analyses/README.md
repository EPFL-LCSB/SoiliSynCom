# Analysis code

The analysis modules retain the current top-level MATLAB scripts and Live
Scripts for in-silico minimal media (iMM), DiMEs/niche-overlap work, and NIS
biosynthetic comparisons. Generated enumerations, model copies, workspaces,
logs, and bundled dependency repositories are excluded.

Identical project-local helpers used by both iMM and DiMEs are stored once in
`shared/`; `duplicate_files.tsv` records their original module aliases.
Unmodified copied helpers from ReMIND/iMM and the MATLAB `issame` utility are
not redistributed. Install the upstream dependencies described in each module
README instead.
