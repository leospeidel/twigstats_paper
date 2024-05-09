#!/bin/bash

i=${SLURM_ARRAY_TASK_ID}

admix="0.00 0.01 0.02 0.03 0.04"
admix=($admix)
a=${admix[$(($i-1))]}
echo $a

file="msprime_ad${a}_${i}"

rm chr.txt
for chr in `seq 1 1 22`
do
  echo $chr >> chr.txt
  python3 ./simulate.py ${chr} ${a} ${file}
  gzip -f ${file}_chr${chr}.vcf
  plink2 --vcf ${file}_chr${chr}.vcf.gz --recode --make-bed --out ${file}
  mv ${file}_tmp.fam ${file}.fam
done

./run_relate.sh ${file}

Rscript mod_bim.R ${file}
Rscript calc.R ${file} 
Rscript calcf4.R ${file}

