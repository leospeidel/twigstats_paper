library(admixtools)
library(twigstats)

file <- commandArgs(trailingOnly = T)[1]

df <- data.frame()
df2 <- data.frame()



f2_blocks_true <- f2_blocks_from_Relate(paste0("",file), paste0("",file), chrs = 1:22, file_map = "../../recomb_rates/genetic_map_combined_b37", blgsize = 0.05, use_muts = T, paste0(file,".poplabels"))
df <- rbind(df, cbind(f4(f2_blocks_true, "Root", "P2", "PX", "P3"), type = "True trees", cutoff = "Inf"))

f2_blocks_Relate <- f2_blocks_from_Relate(paste0("relate_",file), paste0("relate_",file), chrs = 1:22, file_map = "../../recomb_rates/genetic_map_combined_b37", blgsize = 0.05, use_muts = T, paste0(file,".poplabels"))
df <- rbind(df, cbind(f4(f2_blocks_Relate, "Root", "P2", "PX", "P3"), type = "Relate", cutoff = "Inf"))

for(x in c(4e3,6e3,8e3,1e4,2e4,3e4,5e4,7e4,1e5)){

  f2_blocks_true <- f2_blocks_from_Relate(paste0("./",file), paste0("",file), chrs = 1:22, file_map = "../../recomb_rates/genetic_map_combined_b37", blgsize = 0.05, use_muts = T, paste0(file,".poplabels"),t=x)
  df <- rbind(df, cbind(f4(f2_blocks_true, "Root", "P2", "PX", "P3"), type = "True trees", cutoff = sprintf("%d",x))) 

  f2_blocks_Relate <- f2_blocks_from_Relate(paste0("relate_",file), paste0("relate_",file), chrs = 1:22, file_map = "../../recomb_rates/genetic_map_combined_b37", blgsize = 0.05, use_muts = T, paste0(file,".poplabels"), t=x)
  df <- rbind(df, cbind(f4(f2_blocks_Relate, "Root", "P2", "PX", "P3"), type = "Relate", cutoff = sprintf("%d",x)))

}

print(df)

save(df, f2_blocks_true, file = paste0(file,".f4.RData"))

