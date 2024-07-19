library(dplyr)
library(tidyverse)
library(tidyr)
library(umap)

calcf2 <- function(f2_blocks){

	N <- ncol(f2_blocks)
	L <- length(f2_blocks[1,1,])

	weights <- as.numeric(gsub(names(f2_blocks[1,1,]), pattern = "l", replace = ""))

	f2_blocks_2d <- matrix(f2_blocks, nrow = N * N, ncol = L)
	weight_2d <- matrix(weights,ncol = L, nrow = N*N, byrow = T)

	weighted_values <- f2_blocks_2d * weight_2d
	# Compute the sum along the L dimension
	sum_weights <- rowSums(weighted_values)
	# Compute the sum of weights
	total_weight <- sum(weights)
	# Calculate the final result for each pair i, j
	result <- sum_weights / total_weight
	# Reshape the result back to a matrix (N x N)
	result_matrix <- matrix(result, nrow = N, ncol = N)

	return(result_matrix)

}

calc_pca <- function(file, label, labs){

	load(file)
	d <- dim(f2_blocks)[1:2]
	subn <- colnames(f2_blocks)
	subn <- subn[subn != "Root"]

	labs <- subset(labs, ID %in% subn)
	subn <- unique(as.matrix(labs$ID))

	n <- colnames(f2_blocks)
	oname <- c("Han")
	oname <- oname[oname %in% n]
	subn <- intersect(subn,n)

	f2_blocks <- f2_blocks[c(oname,subn),c(oname,subn),]

	A <- calcf2(f2_blocks)
	colnames(A) <- colnames(f2_blocks)
	rownames(A) <- colnames(f2_blocks)
	
  #compute f3 from f2
	N <- ncol(A)
	Ak_i <- matrix(A[oname, ], nrow = N, ncol = N, byrow = F)
	Ak_j <- matrix(A[oname, ], nrow = N, ncol = N, byrow = T)
	A_i_j <- A
	result <- 0.5*(Ak_i + Ak_j - A_i_j)


	result <- cbind(result, pop2 = rownames(result))
	result <- as.data.frame(result)
	result %>% pivot_longer(cols = !pop2,names_to = "pop3", values_to = "est") -> f3_matrix
	f3_matrix     <- subset(f3_matrix, pop2 != oname & pop3 != oname)
	f3_matrix$est <- as.numeric(as.matrix(f3_matrix$est))

	f3_matrix %>% pivot_wider(id_cols = pop2, names_from = pop3, values_from = est) -> df_wider
	colnames(df_wider)[1] <- "pop1"

	names <- unique(df_wider$pop1)
	cols <- names(df_wider[1,])
	df_wider <- df_wider[order(df_wider$pop1),]
	n <- as.matrix(df_wider$pop1)
	df_wider <- df_wider[, c("pop1", n[n %in% subn])]

	df_ret <- data.frame()

	pca <- df_wider %>% 
		select(where(is.numeric)) %>% # retain only numeric columns
		prcomp(scale = T) 
	labs$POP <- factor(labs$POP, levels = n)
	labs <- labs[order(labs$POP),]

	df_ret <- rbind(df_ret, cbind(data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], 
																					 ID = as.matrix(labs$POP[!is.na(labs$POP)]), 
																					 type = label), labs[,-1], method = "PCA"))
	custom.settings = umap.defaults
	custom.settings$n_neighbors = 15
	custom.settings

	umap_fit <- df_wider %>% 
		select(where(is.numeric)) %>% # retain only numeric columns
		scale() %>% umap(config = custom.settings)

	labs$POP <- factor(labs$POP, levels = n)
	labs <- labs[order(labs$POP),]
	df_ret <- rbind(df_ret, cbind(data.frame(PC1 = umap_fit$layout[,1], PC2 = umap_fit$layout[,2], PC3 = NA, 
																					 ID = as.matrix(labs$POP), 
																					 type = label), labs[,-1], method = "UMAP"))

	cmd_fit <- cmdscale(1-df_wider[,-1], eig = T, k = 2)
	labs$POP <- factor(labs$POP, levels = n)
	labs <- labs[order(labs$POP),]
	df_ret <- rbind(df_ret, cbind(data.frame(PC1 = cmd_fit$points[,1], PC2 = cmd_fit$points[,2], PC3 = NA, 
																					 ID = as.matrix(labs$POP), 
																					 type = label), labs[,-1], method = "MDS"))

}



