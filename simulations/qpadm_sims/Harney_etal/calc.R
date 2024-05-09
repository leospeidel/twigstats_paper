library(admixtools)
library(twigstats)
options(scipen=999)

v <- commandArgs(trailingOnly = T)[2]
ad <- commandArgs(trailingOnly = T)[1]
chr <- 1

left = c('5', '9')
right = c('0', '7', '10', '12', '13')
target = '14'

pops <- sort(c(left, right, target))

df <- data.frame()

file_geno <- paste0("./",v,"/std_v",v,"_a",ad,"_chr",chr)

print(file_geno)
if(!file.exists(paste0(file_geno, ".RData"))){
  f2_blocks <- f2_from_geno(file_geno)
  save(f2_blocks, file = paste0(file_geno, ".RData"))
}
load(paste0(file_geno, ".RData"))
qpres <- qpadm(f2_blocks, left = left, right = right, target = target)
df <- rbind(df, cbind(qpres$weights, pval = qpres$popdrop$p[1], type = "genotype", cutoff = "Inf"))


filename <- paste0("./",v,"/relate_std_v",v,"_a",ad)
if(!file.exists(paste0(filename, "_Relate.RData"))){
  f2_blocks <- f2_blocks_from_RelateAges(file_geno, filename, blgsize = 0.05)
  save(f2_blocks, file = paste0(filename, "_Relate.RData"))
}
load(paste0(filename, "_Relate.RData"))
qpres <- qpadm(f2_blocks, left = left, right = right, target = target)
df <- rbind(df, cbind(qpres$weights, pval = qpres$popdrop$p[1], type = "Relate dated", cutoff = "Inf"))

for(x in seq(200,3000,200)){

  filename <- paste0("./",v,"/relate_std_v",v,"_a",ad)
  if(!file.exists(paste0(filename,"_t",x, "_Relate.RData"))){
    f2_blocks <- f2_blocks_from_RelateAges(file_geno, filename, blgsize = 0.05, t = x)
    save(f2_blocks, file = paste0(filename,"_t",x, "_Relate.RData"))
  }
  load(paste0(filename,"_t",x, "_Relate.RData"))
  qpres <- qpadm(f2_blocks, left = left, right = right, target = target)
  df <- rbind(df, cbind(qpres$weights, pval = qpres$popdrop$p[1], type = "Relate dated", cutoff = sprintf("%d",x)))

}

###########

filename <- paste0("./",v,"/relate_std_v",v,"_a",ad)
if(!file.exists(paste0(filename, "_Relate2.RData"))){
  f2_blocks <- f2_blocks_from_Relate(paste0(filename,"_chr1.anc.gz"), paste0(filename,"_chr1.mut.gz"), file_map = "../../recomb_rates/genetic_map_combined_b37_chr1.txt.gz", blgsize = 0.05, poplabels = paste0("std_v1_a0.0_chr1.poplabels"))
  save(f2_blocks, file = paste0(filename, "_Relate2.RData"))
}
load(paste0(filename, "_Relate2.RData"))
qpres <- qpadm(f2_blocks, left = left, right = right, target = target)
df <- rbind(df, cbind(qpres$weights, pval = qpres$popdrop$p[1], type = "Relate trees", cutoff = "Inf"))

for(x in seq(200,3000,200)){

  filename <- paste0("./",v,"/relate_std_v",v,"_a",ad)
  if(!file.exists(paste0(filename,"_t",x, "_Relate2.RData"))){
    f2_blocks <- f2_blocks_from_Relate(paste0(filename,"_chr1.anc.gz"), paste0(filename,"_chr1.mut.gz"), file_map = "../../recomb_rates/genetic_map_combined_b37_chr1.txt.gz", blgsize = 0.05, poplabels = paste0("std_v1_a0.0_chr1.poplabels"), t = x)
    save(f2_blocks, file = paste0(filename,"_t",x, "_Relate2.RData"))
  }
  load(paste0(filename,"_t",x, "_Relate2.RData"))
  qpres <- qpadm(f2_blocks, left = left, right = right, target = target)
  df <- rbind(df, cbind(qpres$weights, pval = qpres$popdrop$p[1], type = "Relate trees", cutoff = sprintf("%d",x)))

}

save(df, file = paste0("qpadm_v",v,"_a",ad,".RData"))

