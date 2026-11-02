rm(list=ls())

library(FactoMineR)
library(GMD)
library(gplots)
library(ecodist)
# install.packages("ClusterR")
library(ClusterR)
library(ggplot2)

load("output/2018-04-18_hclust_MPG_trajectories.Rdata")

# ------------------------------------------------------------------------------
# ALL SAMPLES
# ------------------------------------------------------------------------------

subj.table.all.mat
subj.table.all.mat3

# Get dataframe for all samples, columns ordered by MPG type, not timepoint
# (omit 1.5)
subj.table.all.df4 <- 
	subj.table.all.df3[, c(seq(1,73,12),
												 seq(2,74,12),
												 seq(3,75,12),
												 seq(4,76,12),
												 seq(5,77,12),
												 seq(6,78,12),
												 seq(7,79,12),
												 seq(8,80,12),
												 seq(9,81,12),
												 seq(10,82,12),
												 seq(11,83,12),
												 seq(12,84,12))]
subj.table.all.df4 <- subj.table.all.df4[,-seq(6,83,7)]
subj.table.all.mat4 <- as.matrix(subj.table.all.df4)

# Mean imputation for mat3 and mat4
# Mean imputation is okay because sum of mean MPG rel. ab. within each timepoint 
# should still equal 1 

subj.table.all.mat3.imp <- subj.table.all.mat3
for(i in 1:ncol(subj.table.all.mat3.imp)){
	subj.table.all.mat3.imp[is.na(subj.table.all.mat3.imp[,i]), i] <- 
		mean(subj.table.all.mat3.imp[,i], na.rm = TRUE)
}

for(i in 1:ncol(subj.table.all.mat4)){
 	subj.table.all.mat4[is.na(subj.table.all.mat4[,i]), i] <- 
 		mean(subj.table.all.mat4[,i], na.rm = TRUE)
}
 
colnames(subj.table.all.mat4) <- 
		gsub('([[:punct:]])|\\s+', ".", colnames(subj.table.all.mat4))
colnames(subj.table.all.mat4) <- c("2m.Haemo.2", "4m.Haemo.2", "6m.Haemo.2", "9m.Haemo.2", "1y.Haemo.2", "2y.Haemo.2",
																	 "2m.Strep.1", "4m.Strep.1", "6m.Strep.1", "9m.Strep.1", "1y.Strep.1", "2y.Strep.1",
																	 "2m.Haemo.1", "4m.Haemo.1", "6m.Haemo.1", "9m.Haemo.1", "1y.Haemo.1", "2y.Haemo.1",
																	 "2m.Morax", "4m.Morax", "6m.Morax", "9m.Morax", "1y.Morax", "2y.Morax",
																	 "2m.Strep.2", "4m.Strep.2", "6m.Strep.2", "9m.Strep.2", "1y.Strep.2", "2y.Strep.2",
																	 "2m.Coryn.1", "4m.Coryn.1", "6m.Coryn.1", "9m.Coryn.1", "1y.Coryn.1", "2y.Coryn.1",
																	 "2m.Alloi", "4m.Alloi", "6m.Alloi", "9m.Alloi", "1y.Alloi", "2y.Alloi",
																	 "2m.Staph", "4m.Staph",  "6m.Staph", "9m.Staph", "1y.Staph", "2y.Staph",
																	 "2m.Moraxe.1", "4m.Moraxe.1", "6m.Moraxe.1", "9m.Moraxe.1", "1y.Moraxe.1", "2y.Moraxe.1",
																	 "2m.Moraxe.2", "4m.Moraxe.2", "6m.Moraxe.2", "9m.Moraxe.2", "1y.Moraxe.2", "2y.Moraxe.2",
																	 "2m.Haemo.r", "4m.Haemo.r", "6m.Haemo.r", "9m.Haemo.r", "1y.Haemo.r", "2y.Haemo.r",
																	 "2m.Morax.r", "4m.Morax.r", "6m.Morax.r", "9m.Morax.r", "1y.Morax.r", "2y.Morax.r")

mfa.results <- MFA(subj.table.all.mat4,
									 group = rep(6, 12),
									 type = rep("c", 12),
									 ncp = 5,
									 name.group = c("Haemo2", "Strep1", "Haemo1", "Morax", 
									 							 "Strep2", "Coryn", "Alloi", "Staph", 
									 							 "Moraxe1", "Moraxe2", "Haemo.r", "Morax.r"))

mfa.results <- MFA(subj.table.all.mat4[,-55],
									 group = c(rep(6, 9), 5, rep(6,2)),
									 type = rep("c", 12),
									 ncp = 30,
									 name.group = c("Haemo2", "Strep1", "Haemo1", "Morax",
									 							 "Strep2", "Coryn", "Alloi", "Staph",
									 							 "Moraxe1", "Moraxe2",  "Haemo.r", "Morax.r"))

pdf("output/2018-05-03_mfa.graphs.all.pdf", width=16, height=10)
plot(mfa.results, choix="ind")
plot(mfa.results, choix="var")
plot(mfa.results, choix="group")
plot(mfa.results, choix="axes")
plot((100-(mfa.results$eig)[,3])/100,
		 ylab = "Proportion of variance unexplained",
		 xlab = "Number of dimensions")
dev.off()

# Coordinates of dimensions per individual
mfa.results$ind$coord
coord.mat <- as.matrix(mfa.results$ind$coord)

# Using hclust, choose number of clusters k for kmeans
# ---------------------------------------
dist.p <- dist(coord.mat)
dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
grp <- hclust(dist.p, method = "ward")

# Choose number of clusters based on silhouette
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21),
						 xlab = "Number of clusters", 
						 ylab = "Silhouette")
# Iterate for number of clusters = 2 to 21
for(i in 2:21){
	si <- cluster::silhouette(cutree(grp, k=i), dist.p)
	y <- si[,3]
	m.si[i-1] = mean(y)
	trim.m.si[i-1] = mean(y[y>0])
	med[i-1] = median(y)
	trim.med[i-1] = median(y[y>0])
	b[[i-1]] = boxplot(si[,3], add = TRUE, at = i-1)
}
# Mean silhouette
lines(1:20, m.si, lwd = 2)
# Mean silhouette (positive only)
lines(1:20, trim.m.si, col = 2, lwd = 2)
axis(1, at=1:20, labels = 2:21, cex.axis=1.2)

# Choose k (number of clusters) that maximises median/mean silhouette
which(med == max(med))+1 # 11
which(m.si == max(m.si))+1 # 18
k <- which(med == max(med))+1 
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Kmeans
# ---------------------------------------

# Get optimal number of clusters using ClusterR function
# Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
# 												criterion = 'variance_explained', 
# 												num_init = 100, max_iters = 100,
# 												initializer = 'optimal_init', tol_optimal_init = 0.2)
# Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
# 												criterion = 'WCSSE', 
# 												num_init = 100, max_iters = 100,
# 												initializer = 'optimal_init', tol_optimal_init = 0.2)
# Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
# 												criterion = 'dissimilarity', 
# 												num_init = 100, max_iters = 100,
# 												initializer = 'optimal_init', tol_optimal_init = 0.2)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'silhouette', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
# Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
# 												criterion = 'distortion_fK', fK_threshold = 0.85,
# 												num_init = 100, max_iters = 100,
# 												initializer = 'optimal_init', tol_optimal_init = 0.2)
# Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
# 												criterion = 'AIC', 
# 												num_init = 100, max_iters = 100,
# 												initializer = 'optimal_init', tol_optimal_init = 0.2)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'BIC', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
# Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
# 												criterion = 'Adjusted_Rsquared', 
# 												num_init = 100, max_iters = 100,
# 												initializer = 'optimal_init', tol_optimal_init = 0.2)

k <- 5
kmeans.results <- KMeans_rcpp(coord.mat, clusters=k, 
															num_init=100, max_iters=100,
															initializer="optimal_init")
str(kmeans.results)

pdf("output/2018-05-03_mfa.kmeans.clusters.all.pdf", width=16, height=10)
plot(coord.mat[,1], coord.mat[,2], col=kmeans.results$clusters, 
		 xlim=c(-10,10), ylim=c(-4,8),
		 xlab="Dimension 1", ylab="Dimension 2")
points(kmeans.results$centroids[,1], kmeans.results$centroids[,2], col=c(1:k), 
			 pch=4, cex=5)
dev.off()

# Heatmap for dimensions
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.dimensions.all.pdf", width = 16, height=10)
p1 <- heatmap.3(coord.mat[order(kmeans.results$cluster),], 
							 diss = FALSE, 
							 Colv = FALSE, 
							 Rowv = FALSE,
							 dendrogram = "none", 
							 trace="none", 
							 density.info="none",
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])), 
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# Heatmap with raw data
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.raw.all.pdf", width = 16, height=10)
p2 <- heatmap.3(subj.table.all.mat3[order(kmeans.results$cluster),], 
							 diss = FALSE, 
							 Rowv = FALSE,
							 Colv = FALSE, 
							 dendrogram = "none", 
							 # dist.FUN = dist, 
							 # hclust.FUN = hclust,
							 # hclust.row = hclust(dist.p, method = ("ward")),
							 # kr=k,
							 # cluster.by.row = TRUE, 
							 sepList = list(NULL, list(c(12,0,12.05,282),
							 													c(24,0,24.05,282),
							 													c(36,0,36.05,282),
							 													c(48,0,48.05,282),
							 													c(60,0,60.05,282),
							 													c(72,0,72.05,282),
							 													c(84,0,84.05,282))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(MPG.colour), 
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# Heatmap with mean imputed data
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.imp.all.pdf", width = 16, height=10)
p3 <- heatmap.3(subj.table.all.mat3.imp[order(kmeans.results$cluster),-c(61:72)], 
							 diss = FALSE, 
							 Rowv = FALSE,
							 Colv = FALSE, 
							 dendrogram = "none", 
							 # dist.FUN = dist, 
							 # hclust.FUN = hclust,
							 # hclust.row = hclust(dist.p, method = ("ward")),
							 # kr=k,
							 # cluster.by.row = TRUE, 
							 sepList = list(NULL, list(c(12,0,12.05,282),
							 													c(24,0,24.05,282),
							 													c(36,0,36.05,282),
							 													c(48,0,48.05,282),
							 													c(60,0,60.05,282),
							 													c(72,0,72.05,282),
							 													c(84,0,84.05,282))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(MPG.colour[-c(61:72)]), 
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

pdf("output/2018-05-03_mfa.kmeans.heatmap.imp2.all.pdf", width = 16, height=10)
p4 <- heatmap.3(subj.table.all.mat4[order(kmeans.results$cluster),-c(1,55)], 
							 diss = FALSE, 
							 Rowv = FALSE,
							 Colv = FALSE, 
							 dendrogram = "none", 
							 # dist.FUN = dist, 
							 # hclust.FUN = hclust,
							 # hclust.row = hclust(dist.p, method = ("ward")),
							 # kr=k,
							 # cluster.by.row = TRUE, 
							 sepList = list(NULL, list(c(5,0,5.05,282),
							 													c(11,0,11.05,282),
							 													c(17,0,17.05,282),
							 													c(23,0,23.05,282),
							 													c(29,0,29.05,282),
							 													c(35,0,35.05,282),
							 													c(41,0,41.05,282),
							 													c(47,0,47.05,282),
							 													c(53,0,53.05,282),
							 													c(58,0,58.05,282),
							 													c(64,0,64.05,282),
							 													c(70,0,70.05,282))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(c(
							 	rep(MPG.colour[1],5), rep(MPG.colour[2],6), rep(MPG.colour[3],6),
							 	rep(MPG.colour[4],6), rep(MPG.colour[5],6), rep(MPG.colour[6],6),
							 	rep(MPG.colour[7],6), rep(MPG.colour[8],6), rep(MPG.colour[9],6),
							 	rep(MPG.colour[10],5), rep(MPG.colour[11],6), rep(MPG.colour[12],6))),
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# ------------------------------------------------------------------------------
# ROUTINE SAMPLES ONLY
# ------------------------------------------------------------------------------

subj.table.routine.mat 
dim(subj.table.routine.mat3)

subj.table.routine.df4 <- 
	subj.table.routine.df3[, c(seq(1,73,12),
												 seq(2,74,12),
												 seq(3,75,12),
												 seq(4,76,12),
												 seq(5,77,12),
												 seq(6,78,12),
												 seq(7,79,12),
												 seq(8,80,12),
												 seq(9,81,12),
												 seq(10,82,12),
												 seq(11,83,12),
												 seq(12,84,12))]
subj.table.routine.df4 <- subj.table.routine.df4[,-seq(6,83,7)]
subj.table.routine.mat4 <- as.matrix(subj.table.routine.df4)

# Mean imputation for mat3 and mat4
subj.table.routine.mat3.imp <- subj.table.routine.mat3
for(i in 1:ncol(subj.table.routine.mat3.imp)){
	subj.table.routine.mat3.imp[is.na(subj.table.routine.mat3.imp[,i]), i] <- 
		mean(subj.table.routine.mat3.imp[,i], na.rm = TRUE)
}

for(i in 1:ncol(subj.table.routine.mat4)){
	subj.table.routine.mat4[is.na(subj.table.routine.mat4[,i]), i] <- 
		mean(subj.table.routine.mat4[,i], na.rm = TRUE)
}

colnames(subj.table.routine.mat4) <- 
	gsub('([[:punct:]])|\\s+', ".", colnames(subj.table.routine.mat4))
colnames(subj.table.routine.mat4) <- c("2m.Haemo.2", "4m.Haemo.2", "6m.Haemo.2", "9m.Haemo.2", "1y.Haemo.2", "2y.Haemo.2",
																	 "2m.Strep.1", "4m.Strep.1", "6m.Strep.1", "9m.Strep.1", "1y.Strep.1", "2y.Strep.1",
																	 "2m.Haemo.1", "4m.Haemo.1", "6m.Haemo.1", "9m.Haemo.1", "1y.Haemo.1", "2y.Haemo.1",
																	 "2m.Morax", "4m.Morax", "6m.Morax", "9m.Morax", "1y.Morax", "2y.Morax",
																	 "2m.Strep.2", "4m.Strep.2", "6m.Strep.2", "9m.Strep.2", "1y.Strep.2", "2y.Strep.2",
																	 "2m.Coryn.1", "4m.Coryn.1", "6m.Coryn.1", "9m.Coryn.1", "1y.Coryn.1", "2y.Coryn.1",
																	 "2m.Alloi", "4m.Alloi", "6m.Alloi", "9m.Alloi", "1y.Alloi", "2y.Alloi",
																	 "2m.Staph", "4m.Staph",  "6m.Staph", "9m.Staph", "1y.Staph", "2y.Staph",
																	 "2m.Moraxe.1", "4m.Moraxe.1", "6m.Moraxe.1", "9m.Moraxe.1", "1y.Moraxe.1", "2y.Moraxe.1",
																	 "2m.Moraxe.2", "4m.Moraxe.2", "6m.Moraxe.2", "9m.Moraxe.2", "1y.Moraxe.2", "2y.Moraxe.2",
																	 "2m.Haemo.r", "4m.Haemo.r", "6m.Haemo.r", "9m.Haemo.r", "1y.Haemo.r", "2y.Haemo.r",
																	 "2m.Morax.r", "4m.Morax.r", "6m.Morax.r", "9m.Morax.r", "1y.Morax.r", "2y.Morax.r")

mfa.results <- MFA(subj.table.routine.mat4,
									 group = rep(6, 12),
									 type = rep("c", 12),
									 ncp = 5,
									 name.group = c("Haemo2", "Strep1", "Haemo1", "Morax", 
									 							 "Strep2", "Coryn", "Alloi", "Staph", 
									 							 "Moraxe1", "Moraxe2", "Haemo.r", "Morax.r"))

mfa.results <- MFA(subj.table.routine.mat4[,-c(1,55)],
									 group = c(5,rep(6, 8), 5, rep(6,2)),
									 type = rep("c", 12),
									 ncp = 30,
									 name.group = c("Haemo2", "Strep1", "Haemo1", "Morax",
									 							 "Strep2", "Coryn", "Alloi", "Staph",
									 							 "Moraxe1", "Moraxe2", "Haemo.r", "Morax.r"))

pdf("output/2018-05-03_mfa.graphs.routine.pdf", width=16, height=10)
plot(mfa.results, choix="ind")
plot(mfa.results, choix="var")
plot(mfa.results, choix="group")
plot(mfa.results, choix="axes")
plot((100-(mfa.results$eig)[,3])/100,
		 ylab = "Proportion of variance unexplained",
		 xlab = "Number of dimensions")
dev.off()

# Coordinates of dimensions per individual
mfa.results$ind$coord
coord.mat <- as.matrix(mfa.results$ind$coord)

# Using hclust, choose number of clusters k for kmeans
# ---------------------------------------
dist.p <- dist(coord.mat)
dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
grp <- hclust(dist.p, method = "ward")

pdf("output/2018-05-03_mfa.kmeans.silhouette.routine.pdf", width=16, height=10)
# Choose number of clusters based on silhouette
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21),
						 xlab = "Number of clusters", 
						 ylab = "Silhouette")
# Iterate for number of clusters = 2 to 21
for(i in 2:21){
	si <- cluster::silhouette(cutree(grp, k=i), dist.p)
	y <- si[,3]
	m.si[i-1] = mean(y)
	trim.m.si[i-1] = mean(y[y>0])
	med[i-1] = median(y)
	trim.med[i-1] = median(y[y>0])
	b[[i-1]] = boxplot(si[,3], add = TRUE, at = i-1)
}
# Mean silhouette
lines(1:20, m.si, lwd = 2)
# Mean silhouette (positive only)
lines(1:20, trim.m.si, col = 2, lwd = 2)
axis(1, at=1:20, labels = 2:21, cex.axis=1.2)

# Choose k (number of clusters) that maximises median/mean silhouette
which(med == max(med))+1 # 2
which(m.si == max(m.si))+1 # 2
k <- which(m.si == max(m.si))+1
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Kmeans
# ---------------------------------------
# Get optimal number of clusters using ClusterR function
pdf("output/2018-05-03_mfa.kmeans.optimal.silhouette.routine.pdf", width=16, height=10)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'silhouette', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
dev.off()
pdf("output/2018-05-03_mfa.kmeans.optimal.bic.routine.pdf", width=16, height=10)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'BIC', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
dev.off()

# k <- 7
kmeans.results <- KMeans_rcpp(coord.mat, clusters=k, 
															num_init=100, max_iters=100,
															initializer="optimal_init")
str(kmeans.results)

pdf("output/2018-05-03_mfa.kmeans.clusters.routine.pdf", width=16, height=10)
plot(coord.mat[,1], coord.mat[,2], col=kmeans.results$clusters, 
		 xlim=c(-10,10), ylim=c(-4,8),
		 xlab="Dimension 1", ylab="Dimension 2")
points(kmeans.results$centroids[,1], kmeans.results$centroids[,2], col=c(1:k), 
			 pch=4, cex=5)
dev.off()

# Heatmap for dimensions
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.dimensions.routine.pdf", width = 16, height=10)
p1 <- heatmap.3(coord.mat[order(kmeans.results$cluster),], 
							 diss = FALSE, 
							 Colv = FALSE, 
							 Rowv = FALSE,
							 dendrogram = "none", 
							 trace="none", 
							 density.info="none",
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])), 
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# Heatmap with raw data
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.raw.routine.pdf", width = 16, height=10)
p2 <- heatmap.3(subj.table.routine.mat3[order(kmeans.results$cluster),], 
							 diss = FALSE, 
							 Rowv = FALSE,
							 Colv = FALSE, 
							 dendrogram = "none", 
							 # dist.FUN = dist, 
							 # hclust.FUN = hclust,
							 # hclust.row = hclust(dist.p, method = ("ward")),
							 # kr=k,
							 # cluster.by.row = TRUE, 
							 sepList = list(NULL, list(c(12,0,12.05,282),
							 													c(24,0,24.05,282),
							 													c(36,0,36.05,282),
							 													c(48,0,48.05,282),
							 													c(60,0,60.05,282),
							 													c(72,0,72.05,282),
							 													c(84,0,84.05,282))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(MPG.colour[-seq(60,71)]), 
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# Heatmap with mean imputed data
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.imp.routine.pdf", width = 16, height=10)
p3 <- heatmap.3(subj.table.routine.mat3.imp[order(kmeans.results$cluster),], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sepList = list(NULL, list(c(12,0,12.05,282),
																					c(24,0,24.05,282),
																					c(36,0,36.05,282),
																					c(48,0,48.05,282),
																					c(60,0,60.05,282),
																					c(72,0,72.05,282),
																					c(84,0,84.05,282))),
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(MPG.colour[-c(61:72)]), 
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

pdf("output/2018-05-03_mfa.kmeans.heatmap.imp2.routine.pdf", width = 16, height=10)
p4 <- heatmap.3(subj.table.routine.mat4[order(kmeans.results$cluster),-c(1,55)], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sepList = list(NULL, list(c(5,0,5.05,282),
																					c(11,0,11.05,282),
																					c(17,0,17.05,282),
																					c(23,0,23.05,282),
																					c(29,0,29.05,282),
																					c(35,0,35.05,282),
																					c(41,0,41.05,282),
																					c(47,0,47.05,282),
																					c(53,0,53.05,282),
																					c(58,0,58.05,282),
																					c(64,0,64.05,282),
																					c(70,0,70.05,282))),
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(c(
									rep(MPG.colour[1],5), rep(MPG.colour[2],6), rep(MPG.colour[3],6),
									rep(MPG.colour[4],6), rep(MPG.colour[5],6), rep(MPG.colour[6],6),
									rep(MPG.colour[7],6), rep(MPG.colour[8],6), rep(MPG.colour[9],6),
									rep(MPG.colour[10],5), rep(MPG.colour[11],6), rep(MPG.colour[12],6))),
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

trajectories.mpg.group <- kmeans.results$clusters

# 
# # Kmeans (ALT - KMEANS++)
# # ---------------------------------------
# kmeans.results <- KMeans_rcpp(coord.mat, clusters=k, 
# 															num_init=100, max_iters=100,
# 															initializer="kmeans++")
# str(kmeans.results)
# 
# plot(coord.mat[,1], coord.mat[,2], col=kmeans.results$clusters, 
# 		 xlim=c(-10,10), ylim=c(-4,8))
# points(kmeans.results$centroids[,1], kmeans.results$centroids[,2], col=c(1,2,3), 
# 			 pch=4, cex=5)
# 
# # Heatmap for dimensions
# # ---------------------------------------
# p <- heatmap.3(coord.mat[order(kmeans.results$cluster),], 
# 							 diss = FALSE, 
# 							 Colv = FALSE, 
# 							 Rowv = FALSE,
# 							 dendrogram = "none", 
# 							 trace="none", 
# 							 density.info="none",
# 							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])), 
# 							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
# 							 # color.FUN = seqPalette,
# 							 # labCol = NULL, 
# 							 keysize = 1, mapsize = 9, 
# 							 main = ""
# )
# 
# # Heatmap with raw data
# # ---------------------------------------
# p <- heatmap.3(subj.table.routine.mat3[order(kmeans.results$cluster),], 
# 							 diss = FALSE, 
# 							 Rowv = FALSE,
# 							 Colv = FALSE, 
# 							 dendrogram = "none", 
# 							 # dist.FUN = dist, 
# 							 # hclust.FUN = hclust,
# 							 # hclust.row = hclust(dist.p, method = ("ward")),
# 							 # kr=k,
# 							 # cluster.by.row = TRUE, 
# 							 sepList = list(NULL, list(c(12,0,12.05,282),
# 							 													c(24,0,24.05,282),
# 							 													c(36,0,36.05,282),
# 							 													c(48,0,48.05,282),
# 							 													c(60,0,60.05,282),
# 							 													c(72,0,72.05,282),
# 							 													c(84,0,84.05,282))),
# 							 sep.color="white",
# 							 trace="none", 
# 							 density.info="none",
# 							 ColIndividualColors = as.character(MPG.colour[-seq(60,71)]), 
# 							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
# 							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
# 							 # color.FUN = seqPalette,
# 							 # labCol = NULL, 
# 							 keysize = 1, mapsize = 9, 
# 							 main = ""
# )
# 
# # Heatmap with mean imputed data
# # ---------------------------------------
# p <- heatmap.3(subj.table.routine.mat4[order(kmeans.results$cluster),-c(1,55)], 
# 							 diss = FALSE, 
# 							 Rowv = FALSE,
# 							 Colv = FALSE, 
# 							 dendrogram = "none", 
# 							 # dist.FUN = dist, 
# 							 # hclust.FUN = hclust,
# 							 # hclust.row = hclust(dist.p, method = ("ward")),
# 							 # kr=k,
# 							 # cluster.by.row = TRUE, 
# 							 sepList = list(NULL, list(c(5,0,5.05,282),
# 							 													c(11,0,11.05,282),
# 							 													c(17,0,17.05,282),
# 							 													c(23,0,23.05,282),
# 							 													c(29,0,29.05,282),
# 							 													c(35,0,35.05,282),
# 							 													c(41,0,41.05,282),
# 							 													c(47,0,47.05,282),
# 							 													c(53,0,53.05,282),
# 							 													c(58,0,58.05,282),
# 							 													c(64,0,64.05,282),
# 							 													c(70,0,70.05,282))),
# 							 sep.color="white",
# 							 trace="none", 
# 							 density.info="none",
# 							 ColIndividualColors = as.character(c(
# 							 	rep(MPG.colour[1],5), rep(MPG.colour[2],6), rep(MPG.colour[3],6),
# 							 	rep(MPG.colour[4],6), rep(MPG.colour[5],6), rep(MPG.colour[6],6),
# 							 	rep(MPG.colour[7],6), rep(MPG.colour[8],6), rep(MPG.colour[9],6),
# 							 	rep(MPG.colour[10],5), rep(MPG.colour[11],6), rep(MPG.colour[12],6))),
# 							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
# 							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
# 							 # color.FUN = seqPalette,
# 							 # labCol = NULL, 
# 							 keysize = 1, mapsize = 9, 
# 							 main = ""
# )
# dev.off()

# ------------------------------------------------------------------------------
# HEALTHY ROUTINE SAMPLES ONLY
# ------------------------------------------------------------------------------

subj.table.healthy.mat 
dim(subj.table.healthy.mat3)

subj.table.healthy.df4 <- 
	subj.table.healthy.df3[, c(seq(1,73,12),
														 seq(2,74,12),
														 seq(3,75,12),
														 seq(4,76,12),
														 seq(5,77,12),
														 seq(6,78,12),
														 seq(7,79,12),
														 seq(8,80,12),
														 seq(9,81,12),
														 seq(10,82,12),
														 seq(11,83,12),
														 seq(12,84,12))]
subj.table.healthy.df4 <- subj.table.healthy.df4[,-seq(6,83,7)]
subj.table.healthy.mat4 <- as.matrix(subj.table.healthy.df4)

# Mean imputation for mat3 and mat4
subj.table.healthy.mat3.imp <- subj.table.healthy.mat3
for(i in 1:ncol(subj.table.healthy.mat3.imp)){
	subj.table.healthy.mat3.imp[is.na(subj.table.healthy.mat3.imp[,i]), i] <- 
		mean(subj.table.healthy.mat3.imp[,i], na.rm = TRUE)
}

for(i in 1:ncol(subj.table.healthy.mat4)){
	subj.table.healthy.mat4[is.na(subj.table.healthy.mat4[,i]), i] <- 
		mean(subj.table.healthy.mat4[,i], na.rm = TRUE)
}

colnames(subj.table.healthy.mat4) <- 
	gsub('([[:punct:]])|\\s+', ".", colnames(subj.table.healthy.mat4))
colnames(subj.table.healthy.mat4) <- c("2m.Haemo.2", "4m.Haemo.2", "6m.Haemo.2", "9m.Haemo.2", "1y.Haemo.2", "2y.Haemo.2",
																			 "2m.Strep.1", "4m.Strep.1", "6m.Strep.1", "9m.Strep.1", "1y.Strep.1", "2y.Strep.1",
																			 "2m.Haemo.1", "4m.Haemo.1", "6m.Haemo.1", "9m.Haemo.1", "1y.Haemo.1", "2y.Haemo.1",
																			 "2m.Morax", "4m.Morax", "6m.Morax", "9m.Morax", "1y.Morax", "2y.Morax",
																			 "2m.Strep.2", "4m.Strep.2", "6m.Strep.2", "9m.Strep.2", "1y.Strep.2", "2y.Strep.2",
																			 "2m.Coryn.1", "4m.Coryn.1", "6m.Coryn.1", "9m.Coryn.1", "1y.Coryn.1", "2y.Coryn.1",
																			 "2m.Alloi", "4m.Alloi", "6m.Alloi", "9m.Alloi", "1y.Alloi", "2y.Alloi",
																			 "2m.Staph", "4m.Staph",  "6m.Staph", "9m.Staph", "1y.Staph", "2y.Staph",
																			 "2m.Moraxe.1", "4m.Moraxe.1", "6m.Moraxe.1", "9m.Moraxe.1", "1y.Moraxe.1", "2y.Moraxe.1",
																			 "2m.Moraxe.2", "4m.Moraxe.2", "6m.Moraxe.2", "9m.Moraxe.2", "1y.Moraxe.2", "2y.Moraxe.2",
																			 "2m.Haemo.r", "4m.Haemo.r", "6m.Haemo.r", "9m.Haemo.r", "1y.Haemo.r", "2y.Haemo.r",
																			 "2m.Morax.r", "4m.Morax.r", "6m.Morax.r", "9m.Morax.r", "1y.Morax.r", "2y.Morax.r")

mfa.results <- MFA(subj.table.healthy.mat4,
									 group = rep(6, 12),
									 type = rep("c", 12),
									 ncp = 5,
									 name.group = c("Haemo2", "Strep1", "Haemo1", "Morax", 
									 							 "Strep2", "Coryn", "Alloi", "Staph", 
									 							 "Moraxe1", "Moraxe2", "Haemo.r", "Morax.r"))

mfa.results <- MFA(subj.table.healthy.mat4[,-c(1,2,5,55,65,66)],
									 group = c(3, rep(6, 8), 5, 4, 6),
									 type = rep("c", 12),
									 ncp = 30,
									 name.group = c("Haemo2", "Strep1", "Haemo1", "Morax", 
									 							 "Strep2", "Coryn", "Alloi", "Staph", 
									 							 "Moraxe1", "Moraxe2", "Haemo.r", "Morax.r"))

pdf("output/2018-05-03_mfa.graphs.healthy.pdf", width=16, height=10)
plot(mfa.results, choix="ind")
plot(mfa.results, choix="var")
plot(mfa.results, choix="group")
plot(mfa.results, choix="axes")
plot((100-(mfa.results$eig)[,3])/100,
		 ylab = "Proportion of variance unexplained",
		 xlab = "Number of dimensions")
dev.off()

# Coordinates of dimensions per individual
mfa.results$ind$coord
coord.mat <- as.matrix(mfa.results$ind$coord)

# Using hclust, choose number of clusters k for kmeans
# ---------------------------------------
dist.p <- dist(coord.mat)
dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
grp <- hclust(dist.p, method = "ward")

# Choose number of clusters based on silhouette
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21),
						 xlab = "Number of clusters", 
						 ylab = "Silhouette")
# Iterate for number of clusters = 2 to 21
for(i in 2:21){
	si <- cluster::silhouette(cutree(grp, k=i), dist.p)
	y <- si[,3]
	m.si[i-1] = mean(y)
	trim.m.si[i-1] = mean(y[y>0])
	med[i-1] = median(y)
	trim.med[i-1] = median(y[y>0])
	b[[i-1]] = boxplot(si[,3], add = TRUE, at = i-1)
}
# Mean silhouette
lines(1:20, m.si, lwd = 2)
# Mean silhouette (positive only)
lines(1:20, trim.m.si, col = 2, lwd = 2)
axis(1, at=1:20, labels = 2:21, cex.axis=1.2)

# Choose k (number of clusters) that maximises median/mean silhouette
which(med == max(med))+1 # 17
which(m.si == max(m.si))+1 # 17
k <- which(m.si == max(m.si))+1
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Kmeans
# ---------------------------------------
# Get optimal number of clusters using ClusterR function
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'silhouette', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'BIC', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
# k <- 10
k <- 2
kmeans.results <- KMeans_rcpp(coord.mat, clusters=k, 
															num_init=100, max_iters=100,
															initializer="optimal_init")
str(kmeans.results)

pdf("output/2018-05-03_mfa.kmeans.clusters.healthy.pdf", width=16, height=10)
plot(coord.mat[,1], coord.mat[,2], col=kmeans.results$clusters, 
		 xlim=c(-10,10), ylim=c(-4,8),
		 xlab="Dimension 1", ylab="Dimension 2")
points(kmeans.results$centroids[,1], kmeans.results$centroids[,2], col=c(1:k), 
			 pch=4, cex=5)
dev.off()

# Heatmap for dimensions
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.dimensions.healthy.pdf", width = 16, height=10)
p1 <- heatmap.3(coord.mat[order(kmeans.results$clusters),], 
							 diss = FALSE, 
							 Colv = FALSE, 
							 Rowv = FALSE,
							 dendrogram = "none", 
							 trace="none", 
							 density.info="none",
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])), 
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# Heatmap with raw data
# ---------------------------------------
# subj.table.healthy.mat3 missing row number 214 (row index 193).
healthy.kmeans.clusters <- kmeans.results$clusters[-193]

pdf("output/2018-05-03_mfa.kmeans.heatmap.raw.healthy.pdf", width = 16, height=10)
p2 <- heatmap.3(subj.table.healthy.mat3[order(healthy.kmeans.clusters),], 
							 diss = FALSE, 
							 Rowv = FALSE,
							 Colv = FALSE, 
							 dendrogram = "none", 
							 # dist.FUN = dist, 
							 # hclust.FUN = hclust,
							 # hclust.row = hclust(dist.p, method = ("ward")),
							 # kr=k,
							 # cluster.by.row = TRUE, 
							 sepList = list(NULL, list(c(12,0,12.05,282),
							 													c(24,0,24.05,282),
							 													c(36,0,36.05,282),
							 													c(48,0,48.05,282),
							 													c(60,0,60.05,282),
							 													c(72,0,72.05,282),
							 													c(84,0,84.05,282))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(MPG.colour), 
							 RowIndividualColors = as.character(col2hex(healthy.kmeans.clusters[order(healthy.kmeans.clusters)])),
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# Heatmap with mean imputed data
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.imp.healthy.pdf", width = 16, height=10)
p3 <- heatmap.3(subj.table.all.mat3.imp[order(kmeans.results$cluster),-c(61:72)], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sepList = list(NULL, list(c(12,0,12.05,282),
																					c(24,0,24.05,282),
																					c(36,0,36.05,282),
																					c(48,0,48.05,282),
																					c(60,0,60.05,282),
																					c(72,0,72.05,282),
																					c(84,0,84.05,282))),
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(MPG.colour[-c(61:72)]), 
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

pdf("output/2018-05-03_mfa.kmeans.heatmap.imp2.healthy.pdf", width = 16, height=10)
p4 <- heatmap.3(subj.table.healthy.mat4[order(healthy.kmeans.clusters),-c(1,2,5,55,65,66)], 
							 diss = FALSE, 
							 Rowv = FALSE,
							 Colv = FALSE, 
							 dendrogram = "none", 
							 # dist.FUN = dist, 
							 # hclust.FUN = hclust,
							 # hclust.row = hclust(dist.p, method = ("ward")),
							 # kr=k,
							 # cluster.by.row = TRUE, 
							 sepList = list(NULL, list(c(3,0,3.05,282),
							 													c(9,0,9.05,282),
							 													c(15,0,15.05,282),
							 													c(21,0,21.05,282),
							 													c(27,0,27.05,282),
							 													c(33,0,33.05,282),
							 													c(39,0,39.05,282),
							 													c(45,0,45.05,282),
							 													c(51,0,51.05,282),
							 													c(56,0,56.05,282),
							 													c(60,0,60.05,282),
							 													c(66,0,66.05,282))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(c(
							 	rep(MPG.colour[1],3), rep(MPG.colour[2],6), rep(MPG.colour[3],6),
							 	rep(MPG.colour[4],6), rep(MPG.colour[5],6), rep(MPG.colour[6],6),
							 	rep(MPG.colour[7],6), rep(MPG.colour[8],6), rep(MPG.colour[9],6),
							 	rep(MPG.colour[10],5), rep(MPG.colour[11],4), rep(MPG.colour[12],6))),
							 RowIndividualColors = as.character(col2hex(healthy.kmeans.clusters[order(healthy.kmeans.clusters)])),
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

# ------------------------------------------------------------------------------
# ROUTINE SAMPLES ONLY - group by timepoint
# ------------------------------------------------------------------------------

mfa.results <- MFA(subj.table.routine.mat3.imp,
									 group = rep(12, 6),
									 type = rep("c", 6),
									 ncp = 5,
									 name.group = c("2m", "4m", "6m", "9m", "1y", "2y"))

mfa.results <- MFA(subj.table.routine.mat3.imp[,-c(1,10)],
									 group = c(10, rep(12, 5)),
									 type = rep("c", 6),
									 ncp = 20,
									 name.group = c("2m", "4m", "6m", "9m", "1y", "2y"))

pdf("output/2018-05-03_mfa.graphs.routine.timepoint.pdf", width=16, height=10)
plot(mfa.results, choix="ind")
plot(mfa.results, choix="var")
plot(mfa.results, choix="group")
plot(mfa.results, choix="axes")
plot((100-(mfa.results$eig)[,3])/100,
		 ylab = "Proportion of variance unexplained",
		 xlab = "Number of dimensions")
dev.off()

# Coordinates of dimensions per individual
mfa.results$ind$coord
coord.mat <- as.matrix(mfa.results$ind$coord)

# Using hclust, choose number of clusters k for kmeans
# ---------------------------------------
dist.p <- dist(coord.mat)
dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
grp <- hclust(dist.p, method = "ward")

# Choose number of clusters based on silhouette
pdf("output/2018-05-03_mfa.kmeans.silhouette.routine.timepoint.pdf", width=16, height=10)
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21),
						 xlab = "Number of clusters", 
						 ylab = "Silhouette")
# Iterate for number of clusters = 2 to 21
for(i in 2:21){
	si <- cluster::silhouette(cutree(grp, k=i), dist.p)
	y <- si[,3]
	m.si[i-1] = mean(y)
	trim.m.si[i-1] = mean(y[y>0])
	med[i-1] = median(y)
	trim.med[i-1] = median(y[y>0])
	b[[i-1]] = boxplot(si[,3], add = TRUE, at = i-1)
}
# Mean silhouette
lines(1:20, m.si, lwd = 2)
# Mean silhouette (positive only)
lines(1:20, trim.m.si, col = 2, lwd = 2)
axis(1, at=1:20, labels = 2:21, cex.axis=1.2)

# Choose k (number of clusters) that maximises median/mean silhouette
which(med == max(med))+1 # 5
which(m.si == max(m.si))+1 # 21
k <- which(med == max(med))+1
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Kmeans
# ---------------------------------------
# Get optimal number of clusters using ClusterR function
pdf("output/2018-05-03_mfa.kmeans.optimal.silhouette.routine.timepoint.pdf", width=16, height=10)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'silhouette', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
dev.off()
pdf("output/2018-05-03_mfa.kmeans.optimal.bic.routine.timepoint.pdf", width=16, height=10)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'BIC', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
dev.off()

# k <- 5
kmeans.results <- KMeans_rcpp(coord.mat, clusters=k, 
															num_init=100, max_iters=100,
															initializer="optimal_init")
str(kmeans.results)

pdf("output/2018-05-03_mfa.kmeans.clusters.routine.timepoint.pdf", width=16, height=10)
plot(coord.mat[,1], coord.mat[,2], col=kmeans.results$clusters, 
		 xlim=c(-6,6), ylim=c(-4,3),
		 xlab="Dimension 1", ylab="Dimension 2")
points(kmeans.results$centroids[,1], kmeans.results$centroids[,2], col=c(1:k), 
			 pch=4, cex=5)
dev.off()

# Heatmap for dimensions
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.dimensions.routine.timepoint.pdf", width = 16, height=10)
p1 <- heatmap.3(coord.mat[order(kmeans.results$cluster),], 
								diss = FALSE, 
								Colv = FALSE, 
								Rowv = FALSE,
								dendrogram = "none", 
								trace="none", 
								density.info="none",
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])), 
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

# Heatmap with raw data
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.raw.routine.timepoint.pdf", width = 16, height=10)
p2 <- heatmap.3(subj.table.routine.mat3[order(kmeans.results$cluster),], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sepList = list(NULL, list(c(12,0,12.05,282),
																					c(24,0,24.05,282),
																					c(36,0,36.05,282),
																					c(48,0,48.05,282),
																					c(60,0,60.05,282),
																					c(72,0,72.05,282),
																					c(84,0,84.05,282))),
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(MPG.colour[-seq(60,71)]), 
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

# Heatmap with mean imputed data
# ---------------------------------------
pdf("output/2018-05-03_mfa.kmeans.heatmap.imp.routine.timepoint.pdf", width = 16, height=10)
p3 <- heatmap.3(subj.table.routine.mat3.imp[order(kmeans.results$cluster),], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sepList = list(NULL, list(c(12,0,12.05,282),
																					c(24,0,24.05,282),
																					c(36,0,36.05,282),
																					c(48,0,48.05,282),
																					c(60,0,60.05,282),
																					c(72,0,72.05,282),
																					c(84,0,84.05,282))),
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(MPG.colour[-c(61:72)]), 
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

pdf("output/2018-05-03_mfa.kmeans.heatmap.imp2.routine.timepoint.pdf", width = 16, height=10)
p4 <- heatmap.3(subj.table.routine.mat4[order(kmeans.results$cluster),-c(1,55)], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sepList = list(NULL, list(c(5,0,5.05,282),
																					c(11,0,11.05,282),
																					c(17,0,17.05,282),
																					c(23,0,23.05,282),
																					c(29,0,29.05,282),
																					c(35,0,35.05,282),
																					c(41,0,41.05,282),
																					c(47,0,47.05,282),
																					c(53,0,53.05,282),
																					c(58,0,58.05,282),
																					c(64,0,64.05,282),
																					c(70,0,70.05,282))),
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(c(
									rep(MPG.colour[1],5), rep(MPG.colour[2],6), rep(MPG.colour[3],6),
									rep(MPG.colour[4],6), rep(MPG.colour[5],6), rep(MPG.colour[6],6),
									rep(MPG.colour[7],6), rep(MPG.colour[8],6), rep(MPG.colour[9],6),
									rep(MPG.colour[10],5), rep(MPG.colour[11],6), rep(MPG.colour[12],6))),
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

trajectories.timepoint.group <- kmeans.results$clusters

# ------------------------------------------------------------------------------
# GET GRAPHS FROM CLUSTERS
# ------------------------------------------------------------------------------

df <- data.frame(subj = rownames(subj.table.routine.mat3),
								 subj.table.routine.mat3)
df2 <- melt(df, id.vars = "subj")
df2$age <- factor(c(rep("2m", 283*12), rep("4m", 283*12), rep("6m", 283*12),
						 rep("9m", 283*12), rep("1y", 283*12), rep("2y", 283*12)),
						 levels = c("2m", "4m", "6m", "9m", "1y", "2y"))
df2$mpg <- factor(rep(c(rep("Haemophilus.MPG2", 283), rep("Streptococcus.MPG1", 283), 
						 rep("Haemophilus.MPG1", 283), rep("Moraxella.MPG", 283),
						 rep("Streptococcus.MPG2", 283), rep("Corynebacterium.MPG", 283),
						 rep("Alloiococcus.MPG", 283), rep("Staphylococcus.MPG", 283),
						 rep("Moraxellaceae.MPG1", 283), rep("Moraxellaceae.MPG2", 283),
						 rep("Haemophilus.rare.MPG", 283), rep("Moraxella.rare.MPG", 283)), 6), 
						 levels = c("Haemophilus.MPG2", "Streptococcus.MPG1", "Haemophilus.MPG1", "Moraxella.MPG",
						 					 "Streptococcus.MPG2", "Corynebacterium.MPG", "Alloiococcus.MPG", "Staphylococcus.MPG",
						 					 "Moraxellaceae.MPG1", "Moraxellaceae.MPG2", "Haemophilus.rare.MPG", "Moraxella.rare.MPG"))
trajectories.mpg.group.df <- data.frame(subj = df$subj, 
																				traj.mpg = trajectories.mpg.group)
df3 <- merge(df2, trajectories.mpg.group.df, by = "subj", all = T)

pdf("output/2018-05-03_mfa.kmeans.traj.pdf", width=16, height=10)
df3 %>% filter(!(mpg %in% c("Moraxellaceae.MPG1", "Moraxellaceae.MPG2", 
														"Haemophilus.rare.MPG", "Moraxella.rare.MPG"))) %>%
	mutate(traj.mpg = factor(traj.mpg, labels=c("Traj.1", "Traj.2"))) %>%
	ggplot(aes(x = age, y = value, color = mpg)) + 
	geom_violin() +  
	stat_summary(fun.y="mean", geom="point") +
	stat_summary(aes(group = 1), fun.y=mean, geom="line", linetype=2) +
	# stat_summary(aes(group = 1), fun.data=mean_cl_normal, geom="errorbar", width=0.2) +
	scale_color_manual(values = MPG.colour[1:12]) + 
	facet_grid(mpg~traj.mpg) + 
	labs(y = "Proportion of samples in timepoint", x = "Timepoint", group = "MPG", color = "MPG") + 
	theme(strip.text.y = element_blank())
dev.off()
	

trajectories.timepoint.group.df <- data.frame(subj = df$subj, 
																				traj.timepoint = trajectories.timepoint.group)
df4 <- merge(df2, trajectories.timepoint.group.df, by = "subj", all = T)

pdf("output/2018-05-03_mfa.kmeans.traj.timepoint.pdf", width=16, height=10)
df4  %>% filter(!(mpg %in% c("Moraxellaceae.MPG1", "Moraxellaceae.MPG2", 
														 "Haemophilus.rare.MPG", "Moraxella.rare.MPG"))) %>%
	mutate(traj.timepoint = factor(traj.timepoint, 
																 labels=c("Traj.1", "Traj.2", "Traj.3", 
																 				 "Traj.4", "Traj.5"))) %>%
	ggplot(aes(x = age, y = value, color = mpg)) + 
	geom_violin() +  
	stat_summary(fun.y="mean", geom="point") +
	stat_summary(aes(group = 1), fun.y=mean, geom="line", linetype=2) +
	# stat_summary(aes(group = 1), fun.data=mean_cl_normal, geom="errorbar", width=0.2) +
	scale_color_manual(values = MPG.colour[1:12]) + 
	facet_grid(mpg~traj.timepoint) + 
	labs(y = "Proportion of samples in timepoint", x = "Timepoint", group = "MPG", color = "MPG") + 
	theme(strip.text.y = element_blank())
dev.off()

# ------------------------------------------------------------------------------
# TRAJECTORIES VS. URI/LRI
# ------------------------------------------------------------------------------

df <- merge(log.table.df, trajectories.mpg.group.df, by="subj", all=T)
df$traj.mpg <- factor(df$traj.mpg, 
											labels=c("Traj.1", "Traj.2"))

pdf("output/2018-05-03_mpg_mfa.traj_uri.lri.pdf", width=8, height=5)
df %>% 
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.mpg)) %>% 
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(df$MPG.colour.hex.top)) + 
	facet_grid(traj.mpg~uri.lri) +
	labs(fill = "MPG",
			 x = "Timepoint",
			 y = "Proportion of samples") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) 
dev.off()


# Trajectories vs. no. of LRI/URI samples
# -----------------------------------------

df2 <-df %>%
	filter(!is.na(age.custom), !is.na(well.sick), !is.na(traj.mpg)) %>% 
	group_by(subj, well.sick, traj.mpg) %>%
	summarize(n = n())
# Average number of well, sick samples per subject in Traj, across first 2 years

df3 <- df2 %>%
	ungroup() %>%
	group_by(traj.mpg, well.sick) %>%
	summarize(mean = mean(n),
						sd = sd(n),
						median = median(n, na.rm=T),
						Q1 = quantile(n, na.rm=T)[2],
						Q3 = quantile(n, na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.mpg,
											 well.sick,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), "(",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))

write.csv(df4, "output/2018-05-03_traj_sick_sample_count.csv", 
					quote=F, row.names=F)

# Subject-based analysis
with(df2[which(df2$well.sick == "Well"),], kruskal.test(x=n, 
																												g=traj.mpg)) # p = 0.006
with(df2[which(df2$well.sick == "Sick"),], kruskal.test(x=n, 
																												g=traj.mpg)) # p = 0.89


# Average number of well, uri, lri samples per subject in Traj, across first 2 years
df2 <-df %>%
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.mpg)) %>% 
	group_by(subj, uri.lri, traj.mpg) %>%
	summarize(n = n())
df3 <- df2 %>%
	ungroup() %>%
	group_by(traj.mpg, uri.lri) %>%
	summarize(mean = mean(n),
						sd = sd(n),
						median = median(n, na.rm=T),
						Q1 = quantile(n, na.rm=T)[2],
						Q3 = quantile(n, na.rm=T)[4])

write.csv(df3, "output/2018-05-03_traj_uri.lri_sample_count.csv", 
					quote=F, row.names=F)

# Subject-based analysis
with(df2[which(df2$uri.lri == "Well"),], kruskal.test(x=n, 
																											g=traj.mpg)) # p = 0.006
with(df2[which(df2$uri.lri == "URI"),], kruskal.test(x=n, 
																										 g=traj.mpg)) # p = 0.95
with(df2[which(df2$uri.lri == "LRI"),], kruskal.test(x=n, 
																										 g=traj.mpg)) # p = 0.26

# Average number of well, uri, lri samples per subject in Traj, across first 2 years

# Trajectories vs. severity
# -----------------------------------------
df2 <-df %>%
	filter(!is.na(age.custom), !is.na(score), !is.na(traj.mpg)) %>% 
	group_by(traj.mpg) %>%
	summarize(mean = mean(score),
						sd = sd(score),
						median = median(score, na.rm=T),
						Q1 = quantile(score, na.rm=T)[2],
						Q3 = quantile(score, na.rm=T)[4])

with(df, kruskal.test(x=score, g=traj.mpg)) # p = 0.42

# Trajectories vs. pet ownership
# -----------------------------------------

subject.table.df2 <- merge(subject.table.df, trajectories.mpg.group.df, by="subj", all=T)
subject.table.df2$traj.mpg <- factor(subject.table.df2$traj.mpg, 
											labels=c("Traj.1", "Traj.2"))


# Subject-based analysis
with(subject.table.df2, 
		 table(dog0, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(dog0, traj.mpg, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(dog0, traj.mpg))) # p=0.70

with(subject.table.df2, 
		 table(cat0, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(cat0, traj.mpg, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(cat0, traj.mpg))) # p=0.42

# Trajectories vs. day care
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(dayc, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(dayc, traj.mpg,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.mpg, dayc))) # p=1


# Trajectories vs. number of older siblings
# -----------------------------------------
with(subject.table.df2, hist(livsib))

subject.table.df2 %>%
	group_by(traj.mpg) %>%
	summarize(mean = mean(livsib, na.rm=T),
						sd = sd(livsib, na.rm=T),
						median = median(livsib, na.rm=T),
						Q1 = quantile(livsib, na.rm=T)[2],
						Q3 = quantile(livsib, na.rm=T)[4])

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=livsib, g=traj.routine)) # p = 0.26


# Trajectories vs. breastfeeding
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(brfeed, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(brfeed, traj.mpg,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.mpg, brfeed))) # p=1

# Trajectories vs. aeroallergen sensitisation
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(aero.pattern, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(aero.pattern, traj.mpg,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(aero.pattern, traj.mpg))) # p=0.29

with(subject.table.df2, 
		 kruskal.test(
		 	x=as.numeric(factor(aero.pattern, levels=c("C.0by5", "B.0by1.1+by5", "A.1+by1"))), 
		 	g=traj.mpg)) # p = 0.36

# Trajectories vs. total IgE sensitisation at age 3
# -----------------------------------------
with(subject.table.df2, hist(totige))
with(subject.table.df2, hist(log10(totige)))

df3 <- subject.table.df2 %>%
	group_by(traj.mpg) %>%
	summarize(mean = mean(log10(totige), na.rm=T),
						sd = sd(log10(totige), na.rm=T),
						median = median(log10(totige), na.rm=T),
						Q1 = quantile(log10(totige), na.rm=T)[2],
						Q3 = quantile(log10(totige), na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.mpg,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), " (",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))

df3
df4

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=log10(totige), g=traj.mpg)) # p = 0.009
ggplot(subject.table.df2, aes(x=traj.mpg, y=log10(totige))) + geom_boxplot()

# Trajectories vs. asthma
# -----------------------------------------
# Subject-based analysis
with(subject.table.df2, 
		 table(asthma.yr6, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr6, traj.mpg, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.mpg, asthma.yr6))) # p=0.81

with(subject.table.df2, 
		 table(asthma.yr8, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr8, traj.mpg, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.mpg, asthma.yr8))) # p=0.36

with(subject.table.df2, 
		 table(asthma.yr11, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr11, traj.mpg, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.mpg, asthma.yr11))) # p=0.30

with(subject.table.df2, 
		 table(asthma.yr13, traj.mpg, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr13, traj.mpg, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.mpg, asthma.yr13))) # p=0.29

# Trajectories vs. lung function at age 6
# -----------------------------------------
lung.function.prebd.df <- read.csv("data/melbourne_spirometry_preBD.csv")
lung.function.prebd.df.6y <- lung.function.prebd.df %>%
	filter(whatyear==6)
subject.table.df3 <- merge(subject.table.df2, lung.function.prebd.df.6y, 
													 by="subj", all=T)

with(subject.table.df3, hist(FEV1FVC_percent_predicted/100))

df3 <- subject.table.df3 %>%
	group_by(traj.mpg) %>%
	summarize(mean = mean(FEV1FVC_percent_predicted/100, na.rm=T),
						sd = sd(FEV1FVC_percent_predicted/100, na.rm=T),
						median = median(FEV1FVC_percent_predicted/100, na.rm=T),
						Q1 = quantile(FEV1FVC_percent_predicted/100, na.rm=T)[2],
						Q3 = quantile(FEV1FVC_percent_predicted/100, na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.mpg,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), " (",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))

# Subject-based analysis
with(subject.table.df3, kruskal.test(x=FEV1FVC_percent_predicted/100, g=traj.mpg)) # p = 0.047
ggplot(subject.table.df3, aes(x = traj.mpg, y = FEV1FVC_percent_predicted/100)) + geom_boxplot()
# ------------------------------------------------------------------------------
# VIRUS
# ------------------------------------------------------------------------------

# Create new variables
# -------------------------------------
log.table.df$RV <- as.numeric(grepl("Rhino", log.table.df$virus))
with(log.table.df, table(RV, !is.na(RV.type), useNA="ifany"))
log.table.df[log.table.df$RV == 0 & !is.na(log.table.df$RV.type),]
log.table.df$RV.new <- as.numeric(log.table.df$RV == 1 | !is.na(log.table.df$RV.type))

log.table.df$virus.any <- as.numeric(log.table.df$virus != "neg virus")
log.table.df$Adenovirus <-as.numeric(grepl("Adeno", log.table.df$virus))
log.table.df$Bocavirus <-as.numeric(grepl("Boca|BoV", log.table.df$virus))
log.table.df$CMV <- as.numeric(grepl("CMV", log.table.df$virus))
log.table.df$Coronavirus <- as.numeric(grepl("Corona|CoV", log.table.df$virus))
log.table.df$Coxsackie <- as.numeric(grepl("Cx|Cox", log.table.df$virus))
log.table.df$Echovirus <- as.numeric(grepl("Echo", log.table.df$virus))
log.table.df$Enterovirus <- as.numeric(grepl("Entero|EV", log.table.df$virus))
log.table.df$HPV <- as.numeric(grepl("HPV", log.table.df$virus))
log.table.df$HSV <- as.numeric(grepl("HSV", log.table.df$virus))
log.table.df$Influenza.A <- as.numeric(grepl("Flu A", log.table.df$virus))
log.table.df$Influenza.B <- as.numeric(grepl("Flu B", log.table.df$virus))
log.table.df$Influenza.C <- as.numeric(grepl("Flu C", log.table.df$virus))
log.table.df$Influenza <- as.numeric(log.table.df$Influenza.A == 1 | 
																		 	log.table.df$Influenza.B == 1 |
																		 	log.table.df$Influenza.C == 1 )
log.table.df$MPV <- as.numeric(grepl("MPV", log.table.df$virus))
log.table.df$Parainfluenza <- as.numeric(grepl("Para|PIV", log.table.df$virus))
log.table.df$RSV.A <- as.numeric(grepl("RSV A|RSVA", log.table.df$virus))
log.table.df$RSV.B <- as.numeric(grepl("RSV B|RSVB", log.table.df$virus))
log.table.df$RSV <- as.numeric(grepl("RSV", log.table.df$virus))
log.table.df$TTV <- as.numeric(grepl("TTV", log.table.df$virus))
log.table.df$WU <- as.numeric(grepl("WU Polyomavirus", log.table.df$virus))
log.table.df$RV.A <- as.numeric(grepl("RV-A|RVA", log.table.df$RV.type))
log.table.df$RV.B <- as.numeric(grepl("RV-B|RVB", log.table.df$RV.type))
log.table.df$RV.C <- as.numeric(grepl("RV-C|RVC", log.table.df$RV.type))

log.table.df$Other.vir <- as.numeric(with(log.table.df, virus.any == 1 &
																						CMV + Coxsackie + Echovirus + 
																						Enterovirus + HPV + HSV + TTV + WU > 0))
log.table.df$No.vir <- as.numeric(with(log.table.df, virus.any == 0))

# Get top viruses
# ---------------------------------------
vir.count.table.df <- with(log.table.df, data.frame(LogNo,  
																										RV.new, RV.A, RV.B, RV.C, RSV, Influenza, Adenovirus, 
																										Bocavirus, CMV, Coronavirus, Coxsackie, Echovirus,
																										Enterovirus, HPV, HSV, MPV, Parainfluenza, TTV,WU))

vir.count.table.df.m <- melt(vir.count.table.df, id.var = c("LogNo"),
														 variable.name="virus.name",
														 value.name="count")
vir.count.table.df.m2 <- with(vir.count.table.df.m, 
															as.data.frame.matrix(table(virus.name, count == 1)))
colnames(vir.count.table.df.m2) <- c("Absent", "Present")

vir.count.table.df.m2 <- 
	vir.count.table.df.m2[rev(order(vir.count.table.df.m2$Present)),]

write.csv(vir.count.table.df.m2, "output/2018-02-12_virus_count.csv", quote=F,
					row.names=T)

# Top viruses:
# RV.new, RSV, Parainfluenza, Coronavirus, Adenovirus, Bocavirus, Influenza, MPV, Other

virus.real.name <- c("Adenovirus", "Bocavirus", "RV", "RV-A", "RV-B", "RV-C", 
										 "Coronavirus", "Influenza", "RSV", "Parainfluenza", 
										 "MPV", "Other", "None")
virus.var.name <- c("Adenovirus", "Bocavirus", "RV.new", "RV.A", "RV.B", "RV.C", 
										"Coronavirus", "Influenza",  "RSV", "Parainfluenza",
										"MPV", "Other.vir", "No.vir")
virus.simp.name <- c("adeno", "boca", "RV", "RVA", "RVB", "RVC", 
										 "corona", "influ", "RSV", "para",
										 "MPV", "other", "novir")
virus.colors <- c("#33a02c", "#b2df8a", "#c70002", "#800000", "#fb9a99", "#e31a1c", 
									"#ff7f00", "#1f78b4", "#a6cee3", "#00008B", 
									"#6a3d9a", "#000000", "#808080")

df <- merge(log.table.df, trajectories.mpg.group.df, by="subj", all=T)
df$traj.mpg <- factor(df$traj.mpg, 
											labels=c("Traj.1", "Traj.2"))

df %>% 
	filter(!is.na(age.custom), (visittype == "Well" | visittype == "Well/Sick"), 
				 !is.na(virus.any), !is.na(traj.mpg)) %>% 
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(df$MPG.colour.hex.top)) + 
	facet_grid(traj.mpg~virus.any) +
	labs(fill = "MPG",
			 x = "Timepoint",
			 y = "Proportion of samples") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) 

# Subject-based analysis
with(df, table(traj.mpg,virus.any))
with(df, prop.table(table(traj.mpg,virus.any), 1))
with(df, chisq.test(table(traj.mpg,virus.any))) # p = 0.78
with(df, chisq.test(table(traj.mpg,RV))) # p = 0.78
with(df, chisq.test(table(traj.mpg,RV.A))) # p = 0.99

with(df, prop.table(table(traj.mpg,RV.B), 1))
with(df, chisq.test(table(traj.mpg,RV.B))) # p = 0.029
with(df, chisq.test(table(traj.mpg,RV.C))) # p = 0.13
with(df, chisq.test(table(traj.mpg,RSV))) # p = 0.32
with(df, chisq.test(table(traj.mpg,Influenza))) # p = 0.096

ggplot(df, aes(x = traj.mpg, y = RV.B)) + geom_boxplot()

# ------------------------------------------------------------------------------
# SEASON OF BIRTH
# ------------------------------------------------------------------------------

subject.table.df2$date.birth <- as.Date(subject.table.df2$DOB, "%m/%d/%y")

# Define seasons
 	
# Spring runs from March 1 to May 31;
# Summer runs from June 1 to August 31;
# Fall (autumn) runs from September 1 to November 30; and
# Winter runs from December 1 to February 28 (February 29 in a leap year).

subject.table.df2$month.birth <- as.numeric(format(
	as.POSIXct(subject.table.df2$date.birth), "%m"))
subject.table.df2$season.birth <- factor(
	ifelse(subject.table.df2$month.birth %in% c(3,4,5), "Spring",
				 ifelse(subject.table.df2$month.birth %in% c(6,7,8), "Summer",
				 			 ifelse(subject.table.df2$month.birth %in% c(9,10,11), "Autumn",
				 			 			 ifelse(subject.table.df2$month.birth %in% c(12,1,2), "Winter", NA))))
)
subject.table.df2$season.birth <- factor(subject.table.df2$season.birth, 
															levels = c("Winter","Spring", "Summer", "Autumn"))
subject.table.df2$season.birth2 <- subject.table.df2$season.birth
levels(subject.table.df2$season.birth2) <- c("Winter.Spring", "Winter.Spring",
																	"Summer.Fall", "Summer.Fall")


with(subject.table.df2, table(traj.mpg, season.birth)) 
with(subject.table.df2, prop.table(table(traj.mpg, season.birth), 1)) 
with(subject.table.df2, table(traj.mpg, season.birth2)) 
with(subject.table.df2, prop.table(table(traj.mpg, season.birth2), 1)) 

with(subject.table.df2, chisq.test(table(traj.mpg, season.birth))) # p = 1e-12
with(subject.table.df2, chisq.test(table(traj.mpg, season.birth2))) # p = 0.22	

