#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=31G

i=${SLURM_ARRAY_TASK_ID}

chr=1
admix="0.0 0.2 0.4 0.6 0.8 1.0"
split="250 500"

run(){

  a=$1
  t=$2
  file="msprime_ad${a}_split${t}_${i}"

  python3 ./simulate.py ${chr} ${a} ${t} ${file}
  gzip -f ${file}_chr${chr}.vcf

  plink2 --vcf ${file}_chr${chr}.vcf.gz --recode --make-bed --out ${file}
  cp ${file}_tmp.fam ${file}.fam

  ./run_relate.sh ${file} #run relate
  Rscript mod_bim.R ${file} #add rec rates to bim

  Rscript calc.R ${file} 

  tar -cvf ${file}.tar ${file}_chr${chr}.* relate_${file}_chr${chr}.* ${file}.p* ${file}.b* ${file}.fam ${file}_*gz ${file}.coal ${file}.log
  rm ${file}_chr${chr}.* relate_${file}_chr${chr}.* ${file}.p* ${file}.b* ${file}.fam ${file}_*gz ${file}.coal ${file}.log

}

for a in ${admix}
do
  for t in $split
  do
    run $a $t &
  done
  wait
done


