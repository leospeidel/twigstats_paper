
filename <- commandArgs(trailingOnly = T)[1]
popl <- read.table(paste0(filename,"_chr1.poplabels"), header = T)
fam  <- read.table(paste0(filename, ".fam"))
fam[,1] <- popl[,2]
write.table(fam, file = paste0(filename, ".fam"), row.names = F, col.names = F, quote = F)
