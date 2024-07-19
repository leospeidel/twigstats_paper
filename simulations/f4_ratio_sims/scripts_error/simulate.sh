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

  bcftools view ${file}_error_chr${chr}.vcf -O z -o ${file}_error_chr${chr}.vcf.gz
  bcftools index ${file}_error_chr${chr}.vcf.gz
  shapeit4 --input ${file}_error_chr${chr}.vcf.gz --region 1 --output ${file}_error_phased_chr${chr}.vcf.gz &


  plink2 --vcf ${file}_error_phased_chr${chr}.vcf.gz --recode --make-bed --out ${file}
  cp ${file}_tmp.fam ${file}.fam

  rm ${file}_chr${chr}.vcf
  rm ${file}_error_chr${chr}.vcf

  ./run_relate.sh ${file} #run relate
  Rscript mod_bim.R ${file} #add rec rates to bim

  Rscript calc.R ${file} 

  tar -cvf ${file}.tar ${file}_chr${chr}.* relate_${file}_chr${chr}.* ${file}_error_chr${chr}.* relate_error_${file}_chr${chr}.* ${file}.p* ${file}.b* ${file}.fam ${file}_*gz ${file}.coal ${file}.log
  rm ${file}_chr${chr}.* relate_${file}_chr${chr}.* ${file}_error_chr${chr}.* relate_error_${file}_chr${chr}.* ${file}.p* ${file}.b* ${file}.fam ${file}_*gz ${file}.coal ${file}.log

}

for a in ${admix}
do
  for t in $split
  do
    run $a $t &
  done
  wait
done


