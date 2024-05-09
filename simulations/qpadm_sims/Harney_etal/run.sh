#!/bin/bash

v=$SLURM_ARRAY_TASK_ID
admix="0.0 0.1 0.2 0.3 0.4 0.5"

for ad in $admix
do

  pre=std
  file=${pre}_v${v}_a${ad}

  python3 ./qpAdm_Harney_Supplementary_File_1.py ${v} ${pre} ${ad}
  echo $file

  pushd ./${v}
    ../run_relate.sh ${file}
  popd

  gzip ./${v}/${file}.coal
  gzip ./${v}/${file}_chr1.vcf
  gzip ./${v}/${file}_chr1.geno
  gzip ./${v}/${file}_chr1.snp
  gzip ./${v}/${file}_chr1.ind

  Rscript mod_bim.R ./${v}/${file}_chr1.bim
  Rscript calc.R ${ad} $v

done

