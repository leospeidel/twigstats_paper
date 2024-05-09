options(scipen=999)

file <- commandArgs(trailingOnly = T)[1]

chrom <- 1
rec <- read.table(paste0("../../recomb_rates/HapmapII/genetic_map_GRCh37_chr", chrom, ".txt.gz"), header = T)
f <- Vectorize(approxfun(rec[,2], rec[,4]/100, rule = 2))

bim <- read.table(paste0(file,".bim"))
bim[,3] <- f(bim[,4])
bim[,5] <- 0
bim[,6] <- 1
bim[,1] <- 1
write.table(bim, file = paste0(file,".bim"), row.names = F, col.names =F, quote = F, sep = "\t")

