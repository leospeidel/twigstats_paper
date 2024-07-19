library(twigstats)

run_qpadm <- function(models){

  df_all <- data.frame()

  for(n in c("auto")){
    df <- data.frame()
    df_p <- data.frame()

    if(file.exists(paste0(dir,"/f2_",save_name,".RData"))){
      load(paste0(dir,"/f2_",save_name,".RData"))
      f2_blocks <- f2_blocks[alls,alls,]
      results = qpadm_multi(f2_blocks, models)
      results %>% map('weights') %>% bind_rows(.id = 'model') -> df
      results %>% map('popdrop') %>% bind_rows(.id = 'model') %>% filter(wt == 0) -> df_p

      df %>% group_by(model) %>% mutate(max_se = max(se), min_weight = min(weight), max_weight = max(weight)) -> df
      df <- merge(df, df_p[,c("model", "p")], by = "model")

      df_all <- rbind(df_all, cbind(df, chrom = n, type = "all SNPs"))
    }
  }

  for(cut in c(1000)){
    print(cut)
    for(n in c("auto")){
      df <- data.frame()
      df_p <- data.frame()

      if(file.exists(paste0(dir,"/f2_",save_name,"_",cut,".RData"))){
        load(paste0(dir,"/f2_",save_name,"_",cut,".RData"))
        f2_blocks <- f2_blocks[alls,alls,]
        results = qpadm_multi(f2_blocks, models)
        results %>% map('weights') %>% bind_rows(.id = 'model') -> df
        results %>% map('popdrop') %>% bind_rows(.id = 'model') %>% filter(wt == 0) -> df_p

        df %>% group_by(model) %>% mutate(max_se = max(se), min_weight = min(weight), max_weight = max(weight)) -> df
        df <- merge(df, df_p[,c("model", "p")], by = "model")
        df_all <- rbind(df_all, cbind(df, chrom = n, type = paste0("twigstats ",cut)))
      }
    }
  }

  return(df_all)
}
