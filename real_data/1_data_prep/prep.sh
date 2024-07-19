#!/bin/bash

#SBATCH -c 2
#SBATCH --mem 31G
#SBATCH --partition ncpu
#SBATCH --time 72:00:00

chr=${SLURM_ARRAY_TASK_ID}

SGDP_vcf="/nemo/lab/skoglundp/working/leo/datasets/SGDP/vcfs_1000G/chr.sgdp.pub.${chr}.bcf"
aDNA_vcf="~/leo/datasets/ancients/GL/April2024/${outfile}.phased.chr${chr}.bcf"
aDNA_vcf_with_GP="~/leo/datasets/ancients/GL/April2024/${outfile}.chr${chr}.bcf" #note: same file as aDNA_vcf in GLIMPSE2, different file (unphased) in GLIMPSE1

outfile="aDNA_1000G.April2024.imputed"

#####
#1. merge SGDP with imputed ancients. Take the unfiltered phased vcf from GLIMPSE
bcftools merge ${SGDP_vcf} ${aDNA_vcf} -O z -o SGDP_aDNA_chr${chr}.vcf.gz
bcftools index SGDP_aDNA_chr${chr}.vcf.gz

#subset to samples of interest listed in samples_aDNA.txt and remove any SNPs with missingness just in case.
bcftools view -S samples_aDNA.txt --force-samples --max-alleles 2 --exclude-types indels -i 'F_MISSING==0.0' SGDP_aDNA_chr${chr}.vcf.gz -O z -o SGDP_aDNA_sub_chr${chr}.vcf.gz
bcftools index SGDP_aDNA_sub_chr${chr}.vcf.gz


#####
#2. Take the GLIMPSE vcf that has genotype imputation posteriors and set any genotype with posterior <= 0.8 to missing. Then extract BP positions where more than 2% of samples have a missing genotype
bcftools view -S samples_aDNA.txt --force-samples --max-alleles 2 --exclude-types indels ${aDNA_vcf_with_GP} | bcftools +setGT -- -t q -n . -i'SMPL_MAX(FORMAT/GP)<=0.8' | bcftools annotate -x FORMAT | bcftools view -i "F_MISSING > 0.02" | bcftools query -f "%CHROM\t%POS\n" > excl_0.02_chr${chr}.txt

#set sites in excl_0.02 to missing in the genomic mask. Genomic mask available from e.g., https://www.dropbox.com/scl/fi/ocvr7m00f4jhfcytew4fy/Relate_input_files.tar?rlkey=j1kuqbo2hgvgfwmudwao0v0wv&st=eea0niqm&dl=0
Rscript compile_mask.R ${chr}
gzip -f ./StrictMask_0.02_chr${chr}.fa


####
#3. Compile haps/sample (Relate input)
~/leo/software/relate/bin/RelateFileFormats --mode ConvertFromVcf -i SGDP_aDNA_sub_chr${chr} --haps SGDP_aDNA_chr${chr}.haps --sample SGDP_aDNA_chr${chr}.sample
gzip -f SGDP_aDNA_chr${chr}.haps
gzip -f SGDP_aDNA_chr${chr}.sample

#Filter by genomic mask and polarise using the ancestral genome. E.g. available from https://www.dropbox.com/scl/fi/ocvr7m00f4jhfcytew4fy/Relate_input_files.tar?rlkey=j1kuqbo2hgvgfwmudwao0v0wv&st=eea0niqm&dl=0
~/leo/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps SGDP_aDNA_chr${chr}.haps.gz --sample SGDP_aDNA_chr${chr}.sample.gz --ancestor /nemo/lab/skoglundp/working/leo/datasets/human_genome/human_ancestor_GRCh37_e59/human_ancestor_chr${chr}.fa -o SGDP_aDNA_mask_chr${chr} --remove_ids remove.txt --mask ./StrictMask_0.02_chr${chr}.fa.gz
rm SGDP_aDNA_chr${chr}.haps.gz
rm SGDP_aDNA_chr${chr}.sample.gz


####
#4. Do the same again but now just for transversions
bcftools view -i 'TYPE="snp" && ( (REF!="A" || ALT!="G") && (REF!="G" || ALT!="A") && (REF!="C" || ALT!="T") && (REF!="T" || ALT!="C"))' SGDP_aDNA_sub_chr${chr}.vcf.gz -O z -o SGDP_aDNA_sub_transv_chr${chr}.vcf.gz

~/leo/software/relate/bin/RelateFileFormats --mode ConvertFromVcf -i SGDP_aDNA_sub_transv_chr${chr} --haps SGDP_aDNA_transv_chr${chr}.haps --sample SGDP_aDNA_transv_chr${chr}.sample
gzip -f SGDP_aDNA_transv_chr${chr}.haps
gzip -f SGDP_aDNA_transv_chr${chr}.sample

~/leo/software/relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps SGDP_aDNA_transv_chr${chr}.haps.gz --sample SGDP_aDNA_transv_chr${chr}.sample.gz --ancestor /nemo/lab/skoglundp/working/leo/datasets/human_genome/human_ancestor_GRCh37_e59/human_ancestor_chr${chr}.fa -o SGDP_aDNA_mask_transv_chr${chr} --remove_ids remove.txt --mask ./StrictMask_0.02_chr${chr}.fa.gz
rm SGDP_aDNA_transv_chr${chr}.haps.gz
rm SGDP_aDNA_transv_chr${chr}.sample.gz


