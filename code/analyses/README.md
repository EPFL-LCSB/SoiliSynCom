# Analysis code

We provide analysis modules for:

- `imm/`: in-silico minimal-media calculations;
- `dimes_noi/`: DiMEs enumeration, conditioned-media analysis, and niche
  overlap;
- `nis/`: core, extracellular, and biosynthetic-subnetwork similarity; and
- `essentiality/`: links to the essentiality result data and plotting assets.

We collect project-local functions used across iMM and DiMEs in `shared/`. The
workflows use MATLAB, COBRA Toolbox, matTFA, ReMIND, NIS, redGEM, and compatible
optimization solvers. We provide study inputs and result tables under
[`data/`](../../data/).

`../soilSynComPaths.m` provides repository-relative locations for the input
data and final model families used by these modules.
