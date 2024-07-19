library(abind)

load("f2_Relate_transv_dated_chr1.RData")
d <- dim(f2_blocks)[1:2]
rm(f2_blocks)

concat <- function(file, label){

  print(file)
  f2_blocks_all <- array(,dim = c(d,0))

  if(file.exists(paste0(file,"_chr1.RData"))){
    for(chr in 1:22){
      load(paste0(file,"_chr",chr,".RData"))
      f2_blocks_all <- abind(f2_blocks_all, f2_blocks)
      rm(f2_blocks)
    }
    f2_blocks <- f2_blocks_all
    rm(f2_blocks_all)
  }else{
    load(paste0(file,".RData"))
  }

  save(f2_blocks, file = paste0(label, ".RData"))
  rm(f2_blocks)
  gc()
}

concat("./f2_Relate_transv_dated", "f2_genotypes_new")
concat("./f2_Relate_transv_dated_1000", "f2_genotypes_new_1000")

concat("./f2_Relate_transv", "f2_Relate_new")
concat("./f2_Relate_transv_1000", "f2_Relate_new_1000")

