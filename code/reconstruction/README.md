# Draft reconstruction

We provide the reconstruction commands for the 21 SoilSynCom strains in
`reconstruct_draft_models.m`. The workflow uses RAVEN with the recorded
`1e-20` e-value threshold and accepts configurable input and output roots.
`AddRavenGeneAnnotations.m` transfers FASTA header annotations into a
RAVEN-compatible model structure.

The workflow uses genome FASTA files, KEGG-search results, a RAVEN database
snapshot, and a supported optimization solver.
