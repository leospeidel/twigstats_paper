options(scipen=999)

filename <- commandArgs(trailingOnly = T)[1]

chrom <- 1
rec <- read.table(paste0("../../recomb_rates/HapmapII/genetic_map_GRCh37_chr", chrom, ".txt.gz"), header = T)
f <- Vectorize(approxfun(rec[,2], rec[,4]/100, rule = 2))

bim <- read.table(filename)
bim[,3] <- f(bim[,4])
bim[,2] <- 1:nrow(bim)

write.table(bim, file = filename, row.names = F, col.names = F, quote = F, sep = "\t")



