# Overview

- <b>1_data_prep</b>:
  Starting from imputed ancients, we merge these with present-day genomes, and generate all input files required to run Relate.

- <b>2_relate</b>:
  We use the input files generated in step 1 to run Relate.

- <b>3_f2s</b>:
  We use Twigstats to compute f2 statistics in blocks, using Relate output files. We do this in several ways, in particular grouping individuals according to different poplabels files for different analyses.

- <b>4_outgroup_f3</b>:
  We compute outgroup f3 statistics between all pairs of individuals, and compute an MDS of this variation.

- <b>5_qpAdm</b>:
  We model groups as mixtures of source groups as described in our paper.
