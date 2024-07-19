library(ggplot2)
library(cowplot)
library(admixtools)
library(twigstats)
library(ggh4x)

dir    <- "../3_f2s/"

for(save_name in c("Relate_fine_new", "Relate_dated_fine_new")){

  popl <- read.table("SGDP_aDNA_fine.poplabels", header = T)
  targets <- as.matrix(unique(popl[,2]))
  out_name <- save_name

  v <- 1
  #v1
  if(v == 1){
    refs <- c("Russia_Shamanka_Eneolithic", "Anatolia_EBA", "IrelandOrkney_BA", "Yamnaya")
    ia_sources <- c("Britain.IronRoman","CentralEurope.IronRoman(I)","CentralEurope.IronRoman(II)","CentralEurope.IronRoman(III)",
                    "HungarySlovakia.IronRoman", "Italy.Imperial(I)","Portugal.IronRoman",
                    "Lithuania.IronRoman","PolandUkraine_MLBA(I)","PolandUkraine_MLBA(II)",
                    "Scandinavian_Peninsula_EIA(I)","Scandinavian_Peninsula_EIA(II)",
                    "Russia_Sarmatian", "Kyrgyzstan_TianShanHun"
                    )
    refs <- c(refs, ia_sources)
  }

  cat(paste0("Number of source groups: ", length(ia_sources), "\n"))
  cat(paste0("Number of combinations: ", choose(length(ia_sources),3) + choose(length(ia_sources),2) + length(ia_sources), "\n"))

  load(paste0(dir,"/f2_",save_name,".RData"))

  targets <- targets[targets %in% colnames(f2_blocks)]
  refs <- refs[refs %in% colnames(f2_blocks)]

  print("targets")
  print(targets)
  print("refs")
  print(refs)

  df_sources <- data.frame()
  allrefs <- unique(c(refs))
  alls <- unique(c(targets,refs,ia_sources))

  source("run_qpadm.R")
  df_all <- data.frame()

  ###################
  #1source

  left <- list()
  right <- list()
  target <- list()
  i <- 1
  for(t in targets){
    for(s1 in ia_sources){
      if(t != s1){
        left[[i]] <- c(s1)
        right[[i]] <- as.vector(c(allrefs[!allrefs %in% c(s1,t)]))
        target[[i]] <- t
        df_sources <- rbind(df_sources, data.frame(model = i, source1 = s1, source2 = NA, source3 = NA, source4 = NA, nsources = 1))
        i <- i+1
      }
    }
  }

  models = tibble(left = left,
                  right = right,
                  target = target
                  )

  df_all <- rbind(df_all, cbind(run_qpadm(models), nsources = 1))

  ###################
  #2source
  if(length(targets) > 0){
    left <- list()
    right <- list()
    target <- list()
    i <- 1
    for(t in targets){
      for(s1 in ia_sources){
        for(s2 in ia_sources){
          if(t != s2 & t != s1 & s1 < s2){
            left[[i]] <- c(s1,s2)
            right[[i]] <- as.vector(c(allrefs[!allrefs %in% c(s1,s2,t)]))
            target[[i]] <- t
            df_sources <- rbind(df_sources, data.frame(model = i, source1 = s1, source2 = s2, source3 = NA, source4 = NA, nsources = 2))
            i <- i+1
          }
        }
      }
    }

    models = tibble(left = left,
                    right = right,
                    target = target
                    )

    df_all <- rbind(df_all, cbind(run_qpadm(models), nsources = 2))
  }

  ###################
  #3source

  if(length(targets) > 0){
    left <- list()
    right <- list()
    target <- list()
    i <- 1
    for(t in targets){
      for(s1 in ia_sources){
        for(s2 in ia_sources){
          for(s3 in ia_sources){
            if(t != s2 & t != s1 & t != s3 & s1 < s2 & s2 < s3){
              left[[i]] <- c(s1,s2,s3)
              right[[i]] <- as.vector(c(allrefs[!allrefs %in% c(s1,s2,s3,t)]))
              target[[i]] <- t
              df_sources <- rbind(df_sources, data.frame(model = i, source1 = s1, source2 = s2, source3 = s3, source4 = NA, nsources = 3))
              i <- i + 1
            }
          }
        }
      }
    }

    models = tibble(left = left,
                    right = right,
                    target = target
                    )

    df_all <- rbind(df_all, cbind(run_qpadm(models), nsources = 3))
  }

  df_all$target <- as.matrix(df_all$target)
  df_all <- merge(df_all, df_sources, by = c("model", "nsources"))

  df_all %>% filter(type == "twigstats 1000",min_weight >= 0 & max_weight <= 1 & max_se < 0.5) %>% group_by(target) %>% filter(p > 1e-2 | p == max(p)) -> foo
  print(as.data.frame(foo[order(-foo$p),]))

  file <- paste0("qpadm_2way_cycle_earlymedieval_v",v,"_",out_name,".RData")
  print(file)
  save(df_all, file = file)

}
