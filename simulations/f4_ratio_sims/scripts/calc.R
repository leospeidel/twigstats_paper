library(admixtools)
library(twigstats)
library(lsei)

file <- commandArgs(trailingOnly = T)[1]

df <- data.frame()
df_ratio <- data.frame()

#regular f-stats
f2_blocks <- f2_from_geno(paste0("", file))
df <- rbind(df, cbind(f4(f2_blocks, "P1", "PX", "P3", "P4"), type = "Genotype", cut = Inf))
df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks, "P4", "P1", "P2", "P3", "PX"), type = "Genotype", cut = Inf))

#ascertain by MAF
for(m in seq(0.025,0.45,0.025)){
  f2_blocks <- f2_from_geno(paste0("", file), maxmaf = m)
  df <- rbind(df, cbind(f4(f2_blocks, "P1", "PX", "P3", "P4"), type = "Genotype", cut = m))
  df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks, "P4", "P1", "P2", "P3", "PX"), type = "Genotype", cut = m))
}

########################

#twigstats on Relate trees
f2_blocks_Relate <- f2_blocks_from_Relate(paste0("relate_",file,"_chr1.anc.gz"), paste0("relate_",file,"_chr1.mut.gz"), file_map = "../../recomb_rates/genetic_map_combined_b37_chr1.txt.gz", blgsize = 0.05, poplabels = paste0(file,".poplabels"))
df <- rbind(df, cbind(f4(f2_blocks_Relate, "P1", "PX", "P3", "P4"), type = "Relate trees", cut = Inf))
df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks_Relate, "P4", "P1", "P2", "P3", "PX"), type = "Relate trees", cut = Inf))

cutoffs <- seq(100,5000,100)
for(t in cutoffs){
  f2_blocks_Relate <- f2_blocks_from_Relate(paste0("relate_",file,"_chr1.anc.gz"), paste0("relate_",file,"_chr1.mut.gz"), file_map = "../../recomb_rates/genetic_map_combined_b37_chr1.txt.gz", blgsize = 0.05, poplabels = paste0(file,".poplabels"), t = t)
  df <- rbind(df, cbind(f4(f2_blocks_Relate, "P1", "PX", "P3", "P4"), type = "Relate trees", cut = t))
  df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks_Relate, "P4", "P1", "P2", "P3", "PX"), type = "Relate trees", cut = t))
}

########################

#twigstats on true trees
f2_blocks_Relate <- f2_blocks_from_Relate(paste0(file,"_chr1.anc.gz"), paste0(file,"_chr1.mut.gz"), file_map = "../../recomb_rates/genetic_map_combined_b37_chr1.txt.gz", blgsize = 0.05, poplabels = paste0(file,".poplabels"), use_muts = F)
df <- rbind(df, cbind(f4(f2_blocks_Relate, "P1", "PX", "P3", "P4"), type = "True trees", cut = Inf))
df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks_Relate, "P4", "P1", "P2", "P3", "PX"), type = "True trees", cut = Inf))

for(t in cutoffs){
  f2_blocks_Relate <- f2_blocks_from_Relate(paste0(file,"_chr1.anc.gz"), paste0(file,"_chr1.mut.gz"), file_map = "../../recomb_rates/genetic_map_combined_b37_chr1.txt.gz", blgsize = 0.05, poplabels = paste0(file,".poplabels"), t = t, use_muts = F)
  df <- rbind(df, cbind(f4(f2_blocks_Relate, "P1", "PX", "P3", "P4"), type = "True trees", cut = t))
  df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks_Relate, "P4", "P1", "P2", "P3", "PX"), type = "True trees", cut = t))
}

########################

#twigstats on Relate trees but only ascertaining mutations
f2_blocks_Relate <- f2_blocks_from_RelateAges(pref = paste0("", file), paste0("relate_",file), blgsize = 0.05)
df <- rbind(df, cbind(f4(f2_blocks_Relate, "P1", "PX", "P3", "P4"), type = "Relate dated mutations", cut = Inf))
df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks_Relate, "P4", "P1", "P2", "P3", "PX"), type = "Relate dated mutations", cut = Inf))

#cutoffs <- seq(50,4000,50)
for(t in cutoffs){
  f2_blocks_Relate <- f2_blocks_from_RelateAges(pref = paste0("", file), paste0("relate_",file), blgsize = 0.05, t = t)
  df <- rbind(df, cbind(f4(f2_blocks_Relate, "P1", "PX", "P3", "P4"), type = "Relate dated mutations", cut = t))
  df_ratio <- rbind(df_ratio, cbind(f4_ratio(f2_blocks_Relate, "P4", "P1", "P2", "P3", "PX"), type = "Relate dated mutations", cut = t))
}


#chromosome painting
mat <- ExpPaintingProfile(paste0(file,"_chr1.anc.gz"), paste0(file,"_chr1.mut.gz"), pops = c("P1","P2","P3","P4"), poplabels = paste0(file,".poplabels"))
res <- pnnls(t(mat[c("P2","P3"),c("P1","P2","P3","P4"),1]),(mat["PX",c("P1","P2","P3","P4"),1]), sum = 1)
df_ratio <- rbind(df_ratio, data.frame(popO = "P4", popI = "P1", pop1 = "P2", pop2 = "P3", popX = "PX", val = res$x[1], se = NA, type = "True trees, nnls", cut = Inf))

mat <- ExpPaintingProfile(paste0("relate_",file,"_chr1.anc.gz"), paste0("relate_",file,"_chr1.mut.gz"), pops = c("P1","P2","P3","P4"), poplabels = paste0(file,".poplabels"))
res <- pnnls(t(mat[c("P2","P3"),c("P1","P2","P3","P4"),1]),(mat["PX",c("P1","P2","P3","P4"),1]), sum = 1)
df_ratio <- rbind(df_ratio, data.frame(popO = "P4", popI = "P1", pop1 = "P2", pop2 = "P3", popX = "PX", val = res$x[1], se = NA, type = "Relate trees, nnls", cut = Inf))


save(df, df_ratio, file = paste0(file,".RData"))

