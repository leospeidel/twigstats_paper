#!/bin/bash

#SBATCH -c 12
#SBATCH --mem 95G
#SBATCH --partition ncpu
#SBATCH --time 72:00:00

chr=${SLURM_ARRAY_TASK_ID}

relate_path="/nemo/lab/skoglundp/working/leo/software/"
recomb_rate="../../simulations/recomb_rates/genetic_map_combined_b37_chr${chr}.txt.gz"

####
#Use transversions only
${relate_path}/relate/scripts/RelateParallel/RelateParallel.sh --threads 12 --haps ../1_data_prep/SGDP_aDNA_mask_transv_chr${chr}.haps.gz --sample ../1_data_prep/SGDP_aDNA_mask_transv_chr${chr}.sample.gz --dist ../1_data_prep/SGDP_aDNA_mask_transv_chr${chr}.dist.gz --map ${recomb_rate} -o SGDP_aDNA_mask_transv_new_chr${chr} -m 4e-9 -N 20000 --coal SGDP_v1_annot_ne.coal --sample_ages ../1_data_prep/sample_ages.txt --memory 5 

gzip SGDP_aDNA_mask_transv_new_chr${chr}.anc
gzip SGDP_aDNA_mask_transv_new_chr${chr}.mut

#RelateCoalescentRate --mode EstimatePopulationSize -i SGDP_aDNA_mask_transv_new_chr${chr} -o SGDP_aDNA_mask_transv_new_chr${chr} --poplabels SGDPP_Phase3_sub_aDNA.poplabels --bins 3,7,0.2


####
#Use all SNPs
${relate_path}/relate/scripts/RelateParallel/RelateParallel.sh --threads 12 --haps ../1_data_prep/SGDP_aDNA_mask_chr${chr}.haps.gz --sample ../1_data_prep/SGDP_aDNA_mask_chr${chr}.sample.gz --dist ../1_data_prep/SGDP_aDNA_mask_chr${chr}.dist.gz --map ${recomb_rate} -o SGDP_aDNA_mask_new_chr${chr} -m 1.25e-8 -N 20000 --coal ./SGDP_v1_annot_ne.coal --sample_ages ../1_data_prep/sample_ages.txt --memory 5 --consistency 1

gzip SGDP_aDNA_mask_new_chr${chr}.anc
gzip SGDP_aDNA_mask_new_chr${chr}.mut

#RelateCoalescentRate --mode EstimatePopulationSize -i SGDP_aDNA_mask_new_chr${chr} -o SGDP_aDNA_mask_new_chr${chr} --poplabels ../SGDP_1kG_data/SGDP_sub_aDNA.poplabels --bins 3,7,0.2


