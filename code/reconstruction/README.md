# Draft reconstruction

We provide the reconstruction commands for the 21 SoilSynCom strains in
`reconstruct_draft_models.m`. The workflow uses RAVEN with the recorded
`1e-20` e-value threshold and accepts configurable input and output roots.
`AddRavenGeneAnnotations.m` transfers FASTA header annotations into a
RAVEN-compatible model structure.

Genome metadata are provided under
[`data/workflow_inputs/reconstruction/`](../../data/workflow_inputs/reconstruction/).
The companion Zenodo package contains the protein FASTAs and genome assemblies
under `workflow_inputs/reconstruction/`. The reconstruction entry point reads
that package through `SOILSYNCOM_ARCHIVE_ROOT` and accepts the extracted
KEGG-search output directory as its third argument or through
`SOILSYNCOM_KEGG_OUTPUT_ROOT`.

The workflow uses RAVEN, a KEGG-search snapshot, and a supported optimization
solver.
