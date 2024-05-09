#!/bin/bash

path="~/Documents/genomics/"

for chr in `seq 1 1 1`
do

  ${path}/software/relate_lib/bin/Convert --mode ConvertFromTreeSequence -i msprime_chr${chr}.trees --anc msprime_chr${chr}.anc --mut msprime_chr${chr}.mut
  gzip -f msprime_chr${chr}.anc
  gzip -f msprime_chr${chr}.mut

done

RelateCoalescentRate --mode EstimatePopulationSize -i msprime -o msprime.pairwise --chr chr.txt --poplabels ./msprime.poplabels --bins 3,7,0.2
RelateCoalescentRate --mode CoalRateForTree -i msprime -o msprime --chr chr.txt --bins 3,7,0.2 

for chr in `seq 1 1 1`
do

	if false
	then
  #Run Relate
  RelateFileFormats --mode ConvertFromVcf -i msprime_chr${chr} --haps msprime_chr${chr}.haps --sample msprime_chr${chr}.sample
  gzip -f msprime_chr${chr}.haps
  gzip -f msprime_chr${chr}.sample

	${path}/software//relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps msprime_chr${chr}.haps.gz --sample msprime_chr${chr}.sample.gz --ancestor ../recomb_rates/ancestor.fa.gz -o msprime_masked_chr${chr} --remove_ids remove.txt
	fi

	if true
	then
    Relate --mode All --haps msprime_masked_chr${chr}.haps.gz --sample msprime_masked_chr${chr}.sample.gz --map ../recomb_rates/genetic_map_combined_b37_chr${chr}.txt.gz -m 1.25e-8 -N 20000 -o relate_homsap_chr${chr} --coal msprime.coal 

    gzip -f relate_homsap_chr${chr}.anc
    gzip -f relate_homsap_chr${chr}.mut
	fi

done

RelateCoalescentRate --mode EstimatePopulationSize -i relate_homsap -o relate_homsap --chr chr.txt --poplabels ./msprime.poplabels --bins 3,7,0.2

