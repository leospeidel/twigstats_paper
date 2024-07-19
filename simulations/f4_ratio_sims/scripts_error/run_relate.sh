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

  ${path}/software/relate/bin/RelateFileFormats --mode ConvertFromVcf -i ${file}_error_phased_chr${chr} --haps ${file}_error_chr${chr}.haps --sample ${file}_error_chr${chr}.sample
  gzip -f ${file}_error_chr${chr}.haps
  gzip -f ${file}_error_chr${chr}.sample

  rm ${file}_masked_chr${chr}.*
  ${path}/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps ${file}_error_chr${chr}.haps.gz --sample ${file}_error_chr${chr}.sample.gz --ancestor ancestor.fa.gz -o ${file}_error_masked_chr${chr} --remove_ids remove.txt

  rm -rf relate_error_${file}_chr${chr}
  ${path}/software/relate/bin/Relate --mode All --haps ${file}_error_masked_chr${chr}.haps.gz --sample ${file}_error_masked_chr${chr}.sample.gz --map ${path}/datasets/human_genome/recomb_maps/genetic_map_chr${chr}_combined_b37.txt -m 4e-9 -N 20000 -o relate_error_${file}_chr${chr} --coal ${file}.coal  

  gzip -f relate_error_${file}_chr${chr}.anc
  gzip -f relate_error_${file}_chr${chr}.mut


done

