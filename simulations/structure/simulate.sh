#!/bin/bash

#chr=${SLURM_ARRAY_TASK_ID}
chr=1

python3 ./simulate.py ${chr}
gzip -f msprime_chr${chr}.vcf

