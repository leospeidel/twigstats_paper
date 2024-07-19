library(twigstats)
library(dplyr)

df_all <- data.frame()
for(v in 1:50){
  for(mig in c(0.001, 0.005)){
    load(paste0("./qpadm_v",v,"_mig",mig,".RData"))
    df_all <- rbind(df_all, cbind(df, mig = mig, ver = v))
  }
}

df_all$sources <- paste0(df_all$source1,"-",df_all$source2)
df <- data.frame(sources = c("1-0", "4-3", "5-3", "5-4", "3-1", "4-1", "5-1", "3-0", "4-0", "5-0"), label = c(rep("one-sided",4), rep("two-sided", 6)))
df_all <- merge(df_all, df, by = "sources")
df_all <- subset(df_all, label == "two-sided")

df_all %>% mutate( is_feasible = (weight >= 0 & weight <= 1) ) -> df_all
df_all$pval2 <- df_all$pval
df_all$pval2[df_all$is_feasible == 0] <- 0
df_all %>% filter(left == source2) %>% group_by(type, cutoff, ver, mig) %>% mutate( prank = rank(-pval2), minrank = min(prank) ) -> df_all


df_all %>% group_by(type, cutoff, mig, sources, source1, source2) %>% 
  summarize(median_pval = median(pval), median_pval2 = median(pval2), 
            mean_pval = mean(pval), mean_pval2 = mean(pval2), 
            meanprank = mean(prank), 
            mean_best = mean(prank == minrank),
            mean_rank1 = mean(prank == 1),
            mean_feasible = mean(is_feasible), 
            seprank = sd(prank), 
            se_feasible = sd(is_feasible), 
            median_weight = median(weight[is_feasible]), 
            median_se = median(se[is_feasible]),
            weight = mean(weight[is_feasible]), 
            se = mean(se[is_feasible]) ) -> df_sum

#print(head(subset(df_all, type == "Relate trees" & cutoff == 3000 & source1 == 4 & source2 == 0)))
print(as.data.frame(subset(df_sum, type == "Relate trees" & cutoff == 1000)))

print(as.data.frame(subset(df_sum, type == "Relate trees" & cutoff == 3000)))

print(as.data.frame(subset(df_sum, type == "Relate trees" & cutoff == Inf)))


save(df_sum, df_all, file = "sim_migration.RData")


