library(admixtools)
library(twigstats)
options(scipen=999)

v <- commandArgs(trailingOnly = T)[2]
ad <- commandArgs(trailingOnly = T)[1]
chr <- 1

left = c('p0','p1','p2','p3','p4','p5','p6','p7','p8')
right = c('p11','p10','p9')
targets = left

df <- data.frame()

run_qpadm <- function(f2_blocks, mtype, mcutoff){

  df <- data.frame()
  for(target in targets){
    left2 <- left[!left %in% target]
    for(l1 in left2){
      for(l2 in left2){
        if(l1 > l2){
          right_focus <- c(right, left2[!left2 %in% c(target,l1,l2)])
          print(c(target, l1, l2))
          print(right_focus)
          qpres <- qpadm(f2_blocks, left = c(l1,l2), right = right_focus, target = target)
          df <- rbind(df, cbind(qpres$weights, pval = qpres$popdrop$p[1], target = target, source1 = l1, source2 = l2, type = mtype, cutoff = mcutoff))
        }
      }
    }
  }
  return(df)

}

file_bed <- paste0("./",v,"/std_v",v,"_mig",ad)

if(!file.exists(paste0(file_bed, ".RData"))){
  f2_blocks <- f2_from_geno(file_bed)
  save(f2_blocks, file = paste0(file_bed, ".RData"))
}
load(paste0(file_bed, ".RData"))
df <- rbind(df, run_qpadm(f2_blocks, "genotypes", "Inf"))

filename <- paste0("./",v,"/relate_std_v",v,"_mig",ad)
if(!file.exists(paste0(filename, "_Relate.RData"))){
  f2_blocks <- f2_blocks_from_RelateAges(file_bed, filename, blgsize = 0.05)
  save(f2_blocks, file = paste0(filename, "_Relate.RData"))
}
load(paste0(filename, "_Relate.RData"))
df <- rbind(df, run_qpadm(f2_blocks, "Relate dated", "Inf"))

for(x in c(1000,2000,3000)){

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

for(x in c(1000,2000,3000)){

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

