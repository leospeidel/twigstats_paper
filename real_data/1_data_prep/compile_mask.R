library(seqinr)

chr <- commandArgs(trailingOnly = T)[1]

mask <- read.fasta(paste0("/camp/lab/skoglundp/working/leo/datasets/human_genome/genome_mask/StrictMask_anc/StrictMask_chr",chr,".fa.gz"))[[1]] 

pos <- as.matrix(read.table(paste0("excl_0.02_chr",chr,".txt")))
mask[as.matrix(pos[,2])] <- "n"
print(mean(mask == "p"))

write.fasta(mask, file = paste0("StrictMask_0.02_chr",chr,".fa"), names = paste0("chr",chr))

