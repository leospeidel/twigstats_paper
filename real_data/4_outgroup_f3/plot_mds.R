library(plyr)
library(ggplot2)
library(dplyr)
library(tidyverse)

#poplabels
labs <- read.table("../3_f2s/SGDP_aDNA.poplabels", header = T)
labs <- labs[!duplicated(labs$ID),]

#######################
## Calculate PCA

source("./calc_mds.R")
df <- data.frame()
df <- rbind(df, calc_pca("../3_f2s/f2_genotypes_new.RData", "f3 (genotypes)", labs))
df <- rbind(df, calc_pca("../3_f2s/f2_Relate_new_1000.RData", "f3 (twigstats 1000)", labs))

save(df, file = "f3_pca_ia.RData")

########################
## Plot

load("./f3_pca_ia.RData")
df <- subset(df, method == "MDS")

find_hull <- function(df){ df[chull(df$PC1, df$PC2), ]}
hulls <- ddply(df, c("GROUP", "type"), find_hull)

k <- length(unique(df$GROUPS))

p <- ggplot(df) + 
	geom_polygon(data = hulls, aes(x = PC1, y = PC2, fill = GROUP, colour = GROUP), alpha = 0.2) +         
	geom_point(data = df, aes(x = PC1, y = PC2, colour = GROUP, shape = GROUP, size = coverage)) +
	scale_shape_manual(values = rep(11:19,k), drop = T, name = "") +
	scale_colour_manual(values = rainbow(k), drop = T, name = "") + 
	scale_fill_manual(values = rainbow(k), drop = T, name = "") +
	coord_cartesian(clip = "off") +
	theme_bw() + theme(legend.position = "bottom", legend.box="vertical", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
	facet_wrap(~type, scales = "free") +
	xlab("dim1") + ylab("dim2") 

ggsave(p, file = "plot_MDS.pdf", width = 20, height = 10)



