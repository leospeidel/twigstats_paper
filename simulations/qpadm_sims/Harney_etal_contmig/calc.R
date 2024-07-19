library(admixtools)
library(twigstats)
options(scipen=999)

v <- commandArgs(trailingOnly = T)[2]
ad <- commandArgs(trailingOnly = T)[1]
chr <- 1

left = c('0','1','3','4','5')
right = c('6','7','8')
target = '2'

pops <- sort(c(left, right, target))

df <- data.frame()

run_qpadm <- function(f2_blocks, mtype, mcutoff){

  df <- data.frame()
  for(l1 in left){
    for(l2 in left){
      if(l1 > l2){
        right_focus <- c(right, left[!left %in% c(l1,l2)])
        qpres <- qpadm(f2_blocks, left = c(l1,l2), right = right_focus, target = target)
        df <- rbind(df, cbind(qpres$weights, pval = qpres$popdrop$p[1], source1 = l1, source2 = l2, type = mtype, cutoff = mcutoff))
      }
    }
  }
  return(df)

}

file_geno <- paste0("./",v,"/std_v",v,"_mig",ad,"_chr",chr)
file_bed <- paste0("./",v,"/std_v",v,"_mig",ad)


if(!file.exists(paste0(file_geno, ".RData"))){
  f2_blocks <- f2_from_geno(file_geno)
  save(f2_blocks, file = paste0(file_geno, ".RData"))
}
load(paste0(file_geno, ".RData"))
df <- rbind(df, run_qpadm(f2_blocks, "genotypes", "Inf"))

filename <- paste0("./",v,"/relate_std_v",v,"_mig",ad)
if(!file.exists(paste0(filename, "_Relate.RData"))){
  f2_blocks <- f2_blocks_from_RelateAges(file_bed, filename, blgsize = 0.05)
  save(f2_blocks, file = paste0(filename, "_Relate.RData"))
}
load(paste0(filename, "_Relate.RData"))
df <- rbind(df, run_qpadm(f2_blocks, "Relate dated", "Inf"))

for(x in seq(200,3000,200)){

  filename <- paste0("./",v,"/relate_std_v",v,"_mig",ad)
  if(!file.exists(paste0(filename,"_t",x, "_Relate.RData"))){
    f2_blocks <- f2_blocks_from_RelateAges(file_bed, filename, blgsize = 0.05, t = x)
    save(f2_blocks, file = paste0(filename,"_t",x, "_Relate.RData"))
  }
  load(paste0(filename,"_t",x, "_Relate.RData"))
  df <- rbind(df, run_qpadm(f2_blocks, "Relate dated", sprintf("%d",x)))

}

###########

filename <- paste0("./",v,"/relate_std_v",v,"_mig",ad)
if(!file.exists(paste0(filename, "_Relate2.RData"))){
  f2_blocks <- f2_blocks_from_Relate(paste0(filename,"_chr1.anc.gz"), paste0(filename,"_chr1.mut.gz"), file_map = "/nemo/lab/skoglundp/working/leo/datasets/human_genome/recomb_maps/genetic_map_combined_b37_chr1.txt", blgsize = 0.05, poplabels = paste0("1/std_v1_mig0.001_chr1.poplabels"))
  save(f2_blocks, file = paste0(filename, "_Relate2.RData"))
}
load(paste0(filename, "_Relate2.RData"))
df <- rbind(df, run_qpadm(f2_blocks, "Relate trees", "Inf"))

for(x in seq(200,3000,200)){

  filename <- paste0("./",v,"/relate_std_v",v,"_mig",ad)
  if(!file.exists(paste0(filename,"_t",x, "_Relate2.RData"))){
    f2_blocks <- f2_blocks_from_Relate(paste0(filename,"_chr1.anc.gz"), paste0(filename,"_chr1.mut.gz"), file_map = "/nemo/lab/skoglundp/working/leo/datasets/human_genome/recomb_maps/genetic_map_combined_b37_chr1.txt", blgsize = 0.05, poplabels = paste0("1/std_v1_mig0.001_chr1.poplabels"), t = x)
    save(f2_blocks, file = paste0(filename,"_t",x, "_Relate2.RData"))
  }
  load(paste0(filename,"_t",x, "_Relate2.RData"))
  df <- rbind(df, run_qpadm(f2_blocks, "Relate trees", sprintf("%d",x)))

}


save(df, file = paste0("qpadm_v",v,"_mig",ad,".RData"))

print(df)

