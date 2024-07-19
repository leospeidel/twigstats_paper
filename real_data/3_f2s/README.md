
1. Create a poplabels file. Second column determines population membership. Here we provide two poplabels files, one where every individual is in their own group (used to compute an MDS), and another where individuals are grouped (used for qpadm modelling)

2. Run submit.sh to compute f2 statistics in blocks for every chromosome. This will submit 88 jobs (22 chromosomes times 4 different parameter combinations: twigstats cutoff or no cutoff + Relate trees or Relate mutations )

3. Concatenate chromosomes together by running ```Rscript concat.R``` and ```Rscript concat_fine.R```


