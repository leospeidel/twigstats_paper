library(twigstats)
library(admixtools)

i <- as.numeric(commandArgs(trailingOnly = T)[1])
j <- ceiling(i/44)
mut <- (i %% 2 == 0)
chr <- (ceiling(i/2) %% 22)+1

ancmut_prefix <- "../../SGDP_aDNA_mask_transv_new"
recmap <- "../../simulations/recomb_rates/genetic_map_combined_b37"

poplabels <- "SGDP_aDNA.poplabels"

print(chr)

if(j == 1){

  cat("no cutoff\n")
  if(mut){
    cat("use dated mutations\n")
    f2_blocks <- f2_blocks_from_Relate(file_anc = ancmut_prefix, file_mut = ancmut_prefix, file_map = recmap, chrs = c(chr), poplabels = popl, use_muts = T, transitions = 0)
    save(f2_blocks, file = paste0("f2_Relate_transv_dated_chr",chr,".RData"))
  }else{
    cat("use trees\n")
    f2_blocks <- f2_blocks_from_Relate(file_anc = ancmut_prefix, file_mut = ancmut_prefix, file_map = recmap, chrs = c(chr), poplabels = popl, use_muts = F, transitions = 0)
    save(f2_blocks, file = paste0("f2_Relate_transv_chr",chr,".RData"))
  }

}else{

  thresh <- c(1000)
  cut <- thresh[j-1]

  cat(paste0("cutoff: ", cut,"\n"))
  if(mut){
    cat("use dated mutations\n")
    f2_blocks <- f2_blocks_from_Relate(file_anc = ancmut_prefix, file_mut = ancmut_prefix, file_map = recmap, chrs = c(chr), poplabels = popl, use_muts = T, transitions = 0, t = cut)
    save(f2_blocks, file = paste0("f2_Relate_transv_dated_",cut,"_chr",chr,".RData"))
  }else{
    cat("use trees\n")
    f2_blocks <- f2_blocks_from_Relate(file_anc = ancmut_prefix, file_mut = ancmut_prefix, file_map = recmap, chrs = c(chr), poplabels = popl, use_muts = F, transitions = 0, t = cut)
    save(f2_blocks, file = paste0("f2_Relate_transv_",cut,"_chr",chr,".RData"))
  }

}


