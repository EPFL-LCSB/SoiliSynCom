# Draft reconstruction

`reconstruct_draft_models.m` collects the retained RAVEN reconstruction
commands for the 21 published strains. It preserves the original `1e-20`
e-value threshold and other scoring parameters while accepting portable
input and output roots. `AddRavenGeneAnnotations.m` transfers FASTA header
annotations into a RAVEN-compatible model structure.

The FASTA files, KEGG-search outputs, and RAVEN database snapshot are large
source inputs and are not copied into this GitHub repository. This directory
requires RAVEN and its supported solver/database dependencies.
