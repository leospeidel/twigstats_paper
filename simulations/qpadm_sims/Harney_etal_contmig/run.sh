#!/bin/bash

v=$SLURM_ARRAY_TASK_ID

admix="0.005"
chr=1

for ad in $admix
do

  pre=std
  file=${pre}_v${v}_mig${ad}

  python3 ./qpAdm_Harney_Supplementary_File_5.py ${v} ${pre} ${ad}
  echo $file

  pushd ./${v}
  gzip ${file}_chr${chr}.vcf
  plink2 --vcf ${file}_chr${chr}.vcf.gz --recode --make-bed --out ${file}
  Rscript ../make_fam.R ${file}
  Rscript ../mod_bim.R ${file}

  ../run_relate.sh ${file}
  gzip ${file}.coal

  popd

  Rscript calc.R ${ad} $v

done

