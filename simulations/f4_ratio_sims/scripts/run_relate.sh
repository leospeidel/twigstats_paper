#!/bin/bash

path="/nemo/lab/skoglundp/working/leo/"
file=$1

chrs=$(cat ./chr.txt)

#convert true trees to Relate format
for chr in ${chrs}
do

  ${path}/software/relate_lib/bin/Convert --mode ConvertFromTreeSequence -i ${file}_chr${chr}.trees --compress --anc ${file}_chr${chr}.anc --mut ${file}_chr${chr}.mut
  
  gzip -f ${file}_chr${chr}.anc
  gzip -f ${file}_chr${chr}.mut

done
${path}/software/relate/bin/RelateCoalescentRate --mode CoalRateForTree -i ${file} -o ${file} --chr chr.txt --bins 3,7,0.2 


#run Relate
for chr in ${chrs}
do

  ${path}/software/relate/bin/RelateFileFormats --mode ConvertFromVcf -i ${file}_chr${chr} --haps ${file}_chr${chr}.haps --sample ${file}_chr${chr}.sample
  gzip -f ${file}_chr${chr}.haps
  gzip -f ${file}_chr${chr}.sample

  rm ${file}_masked_chr${chr}.*
  ${path}/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps ${file}_chr${chr}.haps.gz --sample ${file}_chr${chr}.sample.gz --ancestor ancestor.fa.gz -o ${file}_masked_chr${chr} --remove_ids remove.txt
  
  rm -rf relate_${file}_chr${chr}
  ${path}/software/relate/bin/Relate --mode All --haps ${file}_masked_chr${chr}.haps.gz --sample ${file}_masked_chr${chr}.sample.gz --map ../../recomb_rates/genetic_map_combined_b37_chr${chr}.txt.gz -m 1.25e-8 -N 20000 -o relate_${file}_chr${chr} --coal ${file}.coal  

  gzip -f relate_${file}_chr${chr}.anc
  gzip -f relate_${file}_chr${chr}.mut

done

