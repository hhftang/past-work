rm(list=ls())

# Bioconductor-specific packages
# ---------------------------------------
desired.biocon.packages <- c("biomformat",
														 "microbiome", 
														 "metagenomeSeq")

# install new packages
new.biocon.packages <- 	desired.biocon.packages[
	!(desired.biocon.packages %in% installed.packages()[,"Package"])
	]
if(length(new.biocon.packages)>0){
	biocLite(new.biocon.packages)
}

# load new packages
for (i in 1:length(desired.biocon.packages)){
	library(desired.biocon.packages[i], character.only=T)
}

# Standard packages
# ---------------------------------------
desired.packages <- c("dplyr", 
											"ggplot2", 
											"GMD", 
											"ecodist",
											"gtools",
											"ecodist",
											"ClusterR",
											"reshape2",
											"FactoMineR",
											"gee")

# install new packages
new.packages <- desired.packages[
	!(desired.packages %in% installed.packages()[,"Package"])
	]
if(length(new.packages)>0){
	install.packages(new.packages)
}

# load new packages
for (i in 1:length(desired.packages)){
	library(desired.packages[i], character.only=T)
}

load("output/2018-04-18_hclust_MPG_trajectories.Rdata")

theme_set(theme_bw())
# ------------------------------------------------------------------------------
# ROUTINE SAMPLES ONLY
# ------------------------------------------------------------------------------


# Otu table - 2926 samples, 18 common taxa
coast.biom.method1.css.common
# Otu table - 2926 samples, 12464 taxa
coast.biom.method1.css

# Will probably group by timepoint in MFA

# Get otu table for routine samples only
coast.biom.method1.css.common.routine <- 
	subset_samples(coast.biom.method1.css.common,
								 visittype %in% c("Well", "Well/Sick"))

# Otu table - 1488 samples, 18 common taxa
coast.biom.method1.css.common.routine
dim(log.table.df[which(log.table.df$visittype %in% c("Well", "Well/Sick")),])

otu.table.common.routine.df <- 
	data.frame(t(otu_table(coast.biom.method1.css.common.routine)))
x <- gsub(" ", ".", gsub(", ", ".", tax.df$taxa))
colnames(otu.table.common.routine.df) <-substr(x, 1, nchar(x)-28)

otu.table.common.routine.df$LogNo <- row.names(otu.table.common.routine.df)

otu.table.common.routine.df.m <- 
	melt(data = otu.table.common.routine.df,
			 id.vars = "LogNo", variable.name = "otu",
			 value.name = "rel.ab")
	
otu.log.routine.df <- merge(
	otu.table.common.routine.df.m, 
	log.table.df[which(log.table.df$visittype %in% c("Well", "Well/Sick")),],
	by="LogNo", all=T)

subj.table.routine.df <- 
	otu.log.routine.df %>% 
	group_by(age.custom, subj, otu) %>%
	summarize(mean.rel.ab = mean(rel.ab))

subj.table.routine.df$age.otu <-
	with(subj.table.routine.df,
			 interaction(otu, age.custom))

subj.table.routine.df2 <- 
	with(subj.table.routine.df,
			 as.data.frame.matrix(xtabs(mean.rel.ab ~ subj + age.otu)))

# Fill NAs
subj.table.routine.df3 <- subj.table.routine.df2
for (i in 1:7) {
	colstart <- 18*(i-1)+1
	colend <- 18*i
	subj.table.routine.df3[which(rowSums(subj.table.routine.df3[,colstart:colend]) == 0), 
												 colstart:colend] <- NA
}

subj.table.routine.mat <- as.matrix(subj.table.routine.df2)
subj.table.routine.mat3 <- as.matrix(subj.table.routine.df3)

# Missing data stats
# Per column
colSums(is.na(subj.table.routine.df3))[seq(1,126,18)]
(colSums(is.na(subj.table.routine.df3))/nrow(subj.table.routine.df3))[seq(1,126,18)]

# Per row
rowSums(is.na(subj.table.routine.df3))
table(rowSums(is.na(subj.table.routine.df3))/ncol(subj.table.routine.df3)*7)

# Get new colours for OTUs
with(MPG.colour.scheme.df, data.frame(MPG.names, MPG.name2, hex2, descript2,
																			stringsAsFactors = F))
otu.colour.scheme.df <- data.frame(otu = factor(levels(subj.table.routine.df$otu),
																								levels = levels(subj.table.routine.df$otu)),
																	 hex = c("#ff7f00", "#e31a1c", "#b2df8a",
																	 				"#fdbf6f", "#ff8c1a", "#a6cee3", 
																	 				"#cab2db", "#808080", "#ff9933",
																	 				"#ffa54d", "#ffd700", "#b29600",
																	 				"#ff00ff", "#999999", "#ffccff", 
																	 				"#b2b2b2", "#cccccc", "#e5e5e5"), 
																	 stringAsFactors = F)
otu.colour.scheme.df$hex <- as.character(otu.colour.scheme.df$hex)
ggplot(otu.colour.scheme.df, aes(x = otu, fill = otu)) + geom_bar() + 
	scale_fill_manual(values = otu.colour.scheme.df$hex) + 
	theme(axis.text.x = element_text(angle=45, hjust=1))

dist.p <- bcdist(subj.table.routine.mat)

dist.mat <- as.matrix(dist.p)

dim(subj.table.routine.mat) 
dim(subj.table.routine.mat3)

# Exclude 1.5 timepoint
subj.table.routine.mat4 <-  subj.table.routine.mat3[,-c(91:108)]

# Mean imputation for mat3 and mat4
subj.table.routine.mat3.imp <- subj.table.routine.mat3
for(i in 1:ncol(subj.table.routine.mat3.imp)){
	subj.table.routine.mat3.imp[is.na(subj.table.routine.mat3.imp[,i]), i] <- 
		mean(subj.table.routine.mat3.imp[,i], na.rm = TRUE)
}

subj.table.routine.mat4.imp <- subj.table.routine.mat4
for(i in 1:ncol(subj.table.routine.mat4.imp)){
	subj.table.routine.mat4.imp[is.na(subj.table.routine.mat4.imp[,i]), i] <- 
		mean(subj.table.routine.mat4.imp[,i], na.rm = TRUE)
}

mfa.results <- MFA(subj.table.routine.mat4.imp,
									 group = rep(18, 6),
									 type = rep("c", 6),
									 ncp = 24,
									 name.group = c("2mth", "4mth", "6mth", "9mth", 
									 							 "1y", "2y"))

# Rearrange columns
subj.table.routine.mat4.2.imp <- 
	subj.table.routine.mat4.imp[,c(seq(1,108,18),
														 seq(2,108,18),
														 seq(3,108,18),
														 seq(4,108,18),
														 seq(5,108,18),
														 seq(6,108,18),
														 seq(7,108,18),
														 seq(8,108,18),
														 seq(9,108,18),
														 seq(10,108,18),
														 seq(11,108,18),	
														 seq(12,108,18),
														 seq(13,108,18),
														 seq(14,108,18),
														 seq(15,108,18),
														 seq(16,108,18),
														 seq(17,108,18),
														 seq(18,108,18))]

pdf("output/2018-05-18_mfa.graphs.routine.pdf", width=16, height=10)
plot(mfa.results, choix="ind")
plot(mfa.results, choix="var")
plot(mfa.results, choix="group")
plot(mfa.results, choix="axes")
plot((100-(mfa.results$eig)[,3])/100,
		 ylab = "Proportion of variance unexplained",
		 xlab = "Number of dimensions")
dev.off()

plot((100-(mfa.results$eig)[,3])/100,
		 ylab = "Proportion of variance unexplained",
		 xlab = "Number of dimensions")
which((100-(mfa.results$eig)[,3])/100 < 0.2)

# Coordinates of dimensions per individual
mfa.results$ind$coord
coord.mat <- as.matrix(mfa.results$ind$coord)

# Using hclust, choose number of clusters k for kmeans
# ---------------------------------------
dist.p <- dist(coord.mat)
dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
grp <- hclust(dist.p, method = "ward")

pdf("output/2018-05-18_mfa.kmeans.silhouette.routine.pdf", width=16, height=10)
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
which(med == max(med))+1 # 5
which(m.si == max(m.si))+1 # 4
k <- which(m.si == max(m.si))+1 
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Kmeans
# ---------------------------------------
# Get optimal number of clusters using ClusterR function
pdf("output/2018-05-18_mfa.kmeans.optimal.silhouette.routine.pdf", width=16, height=10)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'silhouette', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
dev.off()
pdf("output/2018-05-18_mfa.kmeans.optimal.bic.routine.pdf", width=16, height=10)
Optimal_Clusters_KMeans(coord.mat, max_clusters = 21, plot_clusters = T,
												criterion = 'BIC', 
												num_init = 100, max_iters = 100,
												initializer = 'optimal_init', tol_optimal_init = 0.2)
dev.off()

# k <- 4 
kmeans.results <- KMeans_rcpp(coord.mat, clusters=k, 
															num_init=100, max_iters=100,
															initializer="optimal_init")
str(kmeans.results)

pdf("output/2018-05-18_mfa.kmeans.clusters.routine.pdf", width=16, height=10)
plot(coord.mat[,1], coord.mat[,2], col=kmeans.results$clusters, 
		 xlim=c(-10,10), ylim=c(-4,8),
		 xlab="Dimension 1", ylab="Dimension 2")
points(kmeans.results$centroids[,1], kmeans.results$centroids[,2], col=c(1:k), 
			 pch=4, cex=5)
dev.off()

# Heatmap for dimensions
# ---------------------------------------
pdf("output/2018-05-18_mfa.kmeans.heatmap.dimensions.routine.pdf", width = 16, height=10)
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
pdf("output/2018-05-18_mfa.kmeans.heatmap.raw.routine.pdf", width = 16, height=10)
p2 <- heatmap.3(subj.table.routine.mat4[order(kmeans.results$cluster),], 
							 diss = FALSE, 
							 Rowv = FALSE,
							 Colv = FALSE, 
							 dendrogram = "none", 
							 # dist.FUN = dist, 
							 # hclust.FUN = hclust,
							 # hclust.row = hclust(dist.p, method = ("ward")),
							 # kr=k,
							 # cluster.by.row = TRUE, 
							 # sepList = list(NULL, list(c(12,0,12.05,282),
							 # 													c(24,0,24.05,282),
							 # 													c(36,0,36.05,282),
							 # 													c(48,0,48.05,282),
							 # 													c(60,0,60.05,282),
							 # 													c(72,0,72.05,282),
							 # 													c(84,0,84.05,282))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(rep(otu.colour.scheme.df$hex,6)), 
							 RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()

subj.table.routine.mat4.2 <- 
	subj.table.routine.mat4[,c(seq(1,108,18),
																 seq(2,108,18),
																 seq(3,108,18),
																 seq(4,108,18),
																 seq(5,108,18),
																 seq(6,108,18),
																 seq(7,108,18),
																 seq(8,108,18),
																 seq(9,108,18),
																 seq(10,108,18),
																 seq(11,108,18),	
																 seq(12,108,18),
																 seq(13,108,18),
																 seq(14,108,18),
																 seq(15,108,18),
																 seq(16,108,18),
																 seq(17,108,18),
																 seq(18,108,18))]
pdf("output/2018-05-18_mfa.kmeans.heatmap.raw2.routine.pdf", width = 16, height=10)
p2.1 <- heatmap.3(subj.table.routine.mat4.2[order(kmeans.results$cluster),], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								# sepList = list(NULL, list(c(12,0,12.05,282),
								# 													c(24,0,24.05,282),
								# 													c(36,0,36.05,282),
								# 													c(48,0,48.05,282),
								# 													c(60,0,60.05,282),
								# 													c(72,0,72.05,282),
								# 													c(84,0,84.05,282))),
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(rep(otu.colour.scheme.df$hex,each=6)), 
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
pdf("output/2018-05-18_mfa.kmeans.heatmap.imp.routine.pdf", width = 16, height=10)
p3 <- heatmap.3(subj.table.routine.mat4.imp[order(kmeans.results$cluster),], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(rep(otu.colour.scheme.df$hex,6)), 
								RowIndividualColors = as.character(col2hex(kmeans.results$cluster[order(kmeans.results$cluster)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

pdf("output/2018-05-18_mfa.kmeans.heatmap.imp2.routine.pdf", width = 16, height=10)
p4 <- heatmap.3(subj.table.routine.mat4.2.imp[order(kmeans.results$cluster),], 
								diss = FALSE, 
								Rowv = FALSE,
								Colv = FALSE, 
								dendrogram = "none", 
								# dist.FUN = dist, 
								# hclust.FUN = hclust,
								# hclust.row = hclust(dist.p, method = ("ward")),
								# kr=k,
								# cluster.by.row = TRUE, 
								sep.color="white",
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(rep(otu.colour.scheme.df$hex,each=6)), 
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
df2$age <- factor(c(rep("2m", 283*18), rep("4m", 283*18), rep("6m", 283*18),
						 rep("9m", 283*18), rep("1y", 283*18), rep("1.5y", 283*18), rep("2y", 283*18)),
						 levels = c("2m", "4m", "6m", "9m", "1y", "1.5y", "2y"))
df2$otu <- factor(rep(rep(levels(otu.colour.scheme.df$otu), each=283), 7), 
									levels=levels(otu.colour.scheme.df$otu))

trajectories.timepoint.group.df <- data.frame(subj = df$subj, 
																				traj.timepoint = trajectories.timepoint.group)
df4 <- merge(df2, trajectories.timepoint.group.df, by = "subj", all = T)

pdf("output/2018-05-18_mfa.kmeans.traj.timepoint.pdf", width=16, height=10)
df4  %>% 
	mutate(traj.timepoint = factor(traj.timepoint, 
																 labels=c("Traj.1", "Traj.2", "Traj.3",
																 				 "Traj.4"))) %>%
	ggplot(aes(x = age, y = value, color = otu)) + 
	geom_violin() +  
	stat_summary(fun.y="mean", geom="point") +
	stat_summary(aes(group = 1), fun.y=mean, geom="line", linetype=2) +
	scale_color_manual(values = as.character(otu.colour.scheme.df$hex)) + 
	facet_grid(otu~traj.timepoint) + 
	labs(y = "Relative abundance", x = "Timepoint", group = "OTU", color = "OTU") + 
	theme(strip.text.y = element_blank())
dev.off()

# ------------------------------------------------------------------------------
# TRAJECTORIES VS. URI/LRI
# ------------------------------------------------------------------------------

df <- merge(log.table.df, trajectories.timepoint.group.df, by="subj", all=T)
df$traj.timepoint <- factor(df$traj.timepoint, 
														labels=c("Traj.1", "Traj.2", "Traj.3", "Traj.4"))
	
pdf("output/2018-05-18_mfa.otu.traj_uri.lri.pdf", width=8, height=5)
df %>% 
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.timepoint)) %>% 
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(df$MPG.colour.hex.top)) + 
	facet_grid(uri.lri~traj.timepoint) +
	labs(fill = "MPG",
			 x = "Timepoint",
			 y = "Proportion of samples") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) 
dev.off()


# Trajectories vs. no. of LRI/URI samples
# -----------------------------------------

df2 <-df %>%
	filter(!is.na(age.custom), !is.na(well.sick), !is.na(traj.timepoint)) %>% 
	group_by(subj, well.sick, traj.timepoint) %>%
	summarize(n = n())
# Average number of well, sick samples per subject in Traj, across first 2 years

df3 <- df2 %>%
	ungroup() %>%
	group_by(traj.timepoint, well.sick) %>%
	summarize(mean = mean(n),
						sd = sd(n),
						median = median(n, na.rm=T),
						Q1 = quantile(n, na.rm=T)[2],
						Q3 = quantile(n, na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.timepoint,
											 well.sick,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), "(",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))
df4
write.csv(df4, "output/2018-05-18_traj_sick_sample_count.csv", 
					quote=F, row.names=F)

# Subject-based analysis
with(df2[which(df2$well.sick == "Well"),], kruskal.test(x=n, 
																												g=traj.timepoint)) # p = 0.39
with(df2[which(df2$well.sick == "Sick"),], kruskal.test(x=n, 
																												g=traj.timepoint)) # p = 0.53
with(df2[which(df2$well.sick == "Sick"),], kruskal.test(x=n, 
																												g=traj.timepoint== "Traj.3")) # p = 0.34

# Average number of well, uri, lri samples per subject in Traj, across first 2 years
df2 <-df %>%
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.timepoint)) %>% 
	group_by(subj, uri.lri, traj.timepoint) %>%
	summarize(n = n())
df3 <- df2 %>%
	ungroup() %>%
	group_by(traj.timepoint, uri.lri) %>%
	summarize(mean = mean(n),
						sd = sd(n),
						median = median(n, na.rm=T),
						Q1 = quantile(n, na.rm=T)[2],
						Q3 = quantile(n, na.rm=T)[4])
df3
write.csv(df3, "output/2018-05-18_traj_uri.lri_sample_count.csv", 
					quote=F, row.names=F)

# Subject-based analysis
with(df2[which(df2$uri.lri == "Well"),], kruskal.test(x=n, 
																											g=traj.timepoint)) # p = 0.39
with(df2[which(df2$uri.lri == "URI"),], kruskal.test(x=n, 
																										 g=traj.timepoint)) # p = 0.31
with(df2[which(df2$uri.lri == "LRI"),], kruskal.test(x=n, 
																										 g=traj.timepoint)) # p = 0.60

with(df2[which(df2$uri.lri == "Well"),], kruskal.test(x=n, 
																											g=traj.timepoint == "Traj.3")) # p = 0.80
with(df2[which(df2$uri.lri == "URI"),], kruskal.test(x=n, 
																										 g=traj.timepoint == "Traj.3")) # p = 0.060
with(df2[which(df2$uri.lri == "LRI"),], kruskal.test(x=n, 
																										 g=traj.timepoint == "Traj.3")) # p = 0.23

# Average number of well, uri, lri samples per subject in Traj, across first 2 years

# Trajectories vs. severity
# -----------------------------------------
df2 <-df %>%
	filter(!is.na(age.custom), !is.na(score), !is.na(traj.timepoint)) %>% 
	group_by(traj.timepoint) %>%
	summarize(mean = mean(score),
						sd = sd(score),
						median = median(score, na.rm=T),
						Q1 = quantile(score, na.rm=T)[2],
						Q3 = quantile(score, na.rm=T)[4])
df2
with(df, kruskal.test(x=score, g=traj.timepoint)) # p = 0.002
with(df, kruskal.test(x=score, g=traj.timepoint %in% c("Traj.1"))) # p = 0.017
with(df, kruskal.test(x=score, g=traj.timepoint %in% c("Traj.3"))) # p = 0.057
with(df, kruskal.test(x=score, g=traj.timepoint %in% c("Traj.2"))) # p = 0.0055

# Trajectories vs. pet ownership
# -----------------------------------------

subject.table.df2 <- merge(subject.table.df, trajectories.timepoint.group.df, by="subj", all=T)
subject.table.df2$traj.timepoint <- factor(subject.table.df2$traj.timepoint, 
											labels=c("Traj.1", "Traj.2","Traj.3","Traj.4"))


# Subject-based analysis
with(subject.table.df2, 
		 table(dog0, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(dog0, traj.timepoint, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(dog0, traj.timepoint))) # p=0.25
with(subject.table.df2, 
		 chisq.test(table(dog0, traj.timepoint == "Traj.3"))) # p=0.23


with(subject.table.df2, 
		 table(cat0, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(cat0, traj.timepoint, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(cat0, traj.timepoint))) # p=0.35
with(subject.table.df2, 
		 chisq.test(table(cat0, traj.timepoint == "Traj.3"))) # p=0.72


# Trajectories vs. day care
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(dayc, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(dayc, traj.timepoint,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint, dayc))) # p=0.59
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint == "Traj.3", dayc))) # p=0.99

# Trajectories vs. number of older siblings
# -----------------------------------------
with(subject.table.df2, hist(livsib))

subject.table.df2 %>%
	group_by(traj.timepoint) %>%
	summarize(mean = mean(livsib, na.rm=T),
						sd = sd(livsib, na.rm=T),
						median = median(livsib, na.rm=T),
						Q1 = quantile(livsib, na.rm=T)[2],
						Q3 = quantile(livsib, na.rm=T)[4])

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=livsib, g=traj.timepoint)) # p = 0.16
with(subject.table.df2, kruskal.test(x=livsib, g=traj.timepoint == "Traj.3")) # p-0.42
with(subject.table.df2, kruskal.test(x=livsib, g=traj.timepoint == "Traj.2")) # p-0.07

# Trajectories vs. breastfeeding
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(brfeed, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(brfeed, traj.timepoint,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint, brfeed))) # p=0.22
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint == "Traj.3", brfeed))) # p=0.062

# Trajectories vs. aeroallergen sensitisation
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(aero.pattern, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(aero.pattern, traj.timepoint,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(aero.pattern, traj.timepoint))) # p=0.59
with(subject.table.df2, 
		 chisq.test(table(aero.pattern, traj.timepoint == "Traj.3"))) # p=0.20

with(subject.table.df2, 
		 kruskal.test(
		 	x=as.numeric(factor(aero.pattern, levels=c("C.0by5", "B.0by1.1+by5", "A.1+by1"))), 
		 	g=traj.timepoint)) # p = 0.40
with(subject.table.df2, 
		 kruskal.test(
		 	x=as.numeric(factor(aero.pattern, levels=c("C.0by5", "B.0by1.1+by5", "A.1+by1"))), 
		 	g=traj.timepoint == "Traj.3")) # p = 0.14

# Trajectories vs. total IgE sensitisation at age 3
# -----------------------------------------
with(subject.table.df2, hist(totige))
with(subject.table.df2, hist(log10(totige)))

df3 <- subject.table.df2 %>%
	group_by(traj.timepoint) %>%
	summarize(mean = mean(log10(totige), na.rm=T),
						sd = sd(log10(totige), na.rm=T),
						median = median(log10(totige), na.rm=T),
						Q1 = quantile(log10(totige), na.rm=T)[2],
						Q3 = quantile(log10(totige), na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.timepoint,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), " (",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))

df4

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=log10(totige), g=traj.timepoint)) # p = 0.13
with(subject.table.df2, kruskal.test(x=log10(totige), 
																		 g=traj.timepoint == "Traj.1")) # p = 0.02
with(subject.table.df2, kruskal.test(x=log10(totige), 
																		 g=traj.timepoint == "Traj.3")) # p = 0.97
ggplot(subject.table.df2, aes(x=traj.timepoint, y=log10(totige))) + geom_boxplot()

# Trajectories vs. gender
# -----------------------------------------
with(subject.table.df2, 
		 table(gender, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(gender, traj.timepoint),2))
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint, gender))) # p=0.45
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint == "Traj.3", gender))) # p=0.45

# Trajectories vs. asthma
# -----------------------------------------
# Subject-based analysis
with(subject.table.df2, 
		 table(asthma.yr6, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr6, traj.timepoint),2))
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint, asthma.yr6))) # p=0.21
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint == "Traj.3", asthma.yr6))) # p=0.11

with(subject.table.df2, 
		 table(asthma.yr8, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr8, traj.timepoint),2))
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint, asthma.yr8))) # p=0.093
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint == "Traj.3", asthma.yr8))) # p=0.026

with(subject.table.df2, 
		 table(asthma.yr11, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr11, traj.timepoint),2))
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint, asthma.yr11))) # p=0.042
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint == "Traj.3", asthma.yr11))) # p=0.014

with(subject.table.df2, 
		 table(asthma.yr13, traj.timepoint, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr13, traj.timepoint),2))
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint, asthma.yr13))) # p=0.041
with(subject.table.df2, 
		 chisq.test(table(traj.timepoint == "Traj.3", asthma.yr13))) # p=0.028

# Trajectories vs. lung function at age 6
# -----------------------------------------
lung.function.prebd.df <- read.csv("data/melbourne_spirometry_preBD.csv")
lung.function.prebd.df.6y <- lung.function.prebd.df %>%
	filter(whatyear==6)
subject.table.df3 <- merge(subject.table.df2, lung.function.prebd.df.6y, 
													 by="subj", all=T)

with(subject.table.df3, hist(FEV1FVC_percent_predicted/100))

df3 <- subject.table.df3 %>%
	group_by(traj.timepoint) %>%
	summarize(mean = mean(FEV1FVC_percent_predicted/100, na.rm=T),
						sd = sd(FEV1FVC_percent_predicted/100, na.rm=T),
						median = median(FEV1FVC_percent_predicted/100, na.rm=T),
						Q1 = quantile(FEV1FVC_percent_predicted/100, na.rm=T)[2],
						Q3 = quantile(FEV1FVC_percent_predicted/100, na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.timepoint,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), " (",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))
df4
# Subject-based analysis
with(subject.table.df3, kruskal.test(x=FEV1FVC_percent_predicted/100, g=traj.timepoint)) # p = 0.36
with(subject.table.df3, kruskal.test(x=FEV1FVC_percent_predicted/100, 
																		 g=traj.timepoint == "Traj.3")) # p = 0.13
ggplot(subject.table.df3, aes(x = traj.timepoint, y = FEV1FVC_percent_predicted/100)) + geom_boxplot()

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

df <- merge(log.table.df, trajectories.timepoint.group.df, by="subj", all=T)
df$traj.timepoint <- factor(df$traj.timepoint, 
											labels=c("Traj.1", "Traj.2", "Traj.3", "Traj.4"))

df %>% 
	filter(!is.na(age.custom), (visittype == "Well" | visittype == "Well/Sick"), 
				 !is.na(virus.any), !is.na(traj.timepoint)) %>% 
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(df$MPG.colour.hex.top)) + 
	facet_grid(virus.any~traj.timepoint) +
	labs(fill = "MPG",
			 x = "Timepoint",
			 y = "Proportion of samples") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) 

# Subject-based analysis
with(df, table(traj.timepoint,virus.any))
with(df, prop.table(table(traj.timepoint,virus.any), 1))
with(df, chisq.test(table(traj.timepoint,virus.any))) # p = 0.25
with(df, chisq.test(table(traj.timepoint,RV))) # p = 0.51
with(df, chisq.test(table(traj.timepoint,RV.A))) # p = 0.79

with(df, prop.table(table(traj.timepoint,RV.B), 1))
with(df, chisq.test(table(traj.timepoint,RV.B))) # p = 0.043
with(df, chisq.test(table(traj.timepoint == "Traj.3",RV.B))) # p = 0.024

with(df, prop.table(table(traj.timepoint,RV.C), 1))
with(df, chisq.test(table(traj.timepoint,RV.C))) # p = 0.12
with(df, chisq.test(table(traj.timepoint == "Traj.3",RV.C))) # p = 0.02

with(df, chisq.test(table(traj.timepoint,RSV))) # p = 0.21
with(df, chisq.test(table(traj.timepoint,Influenza))) # p = 0.26

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


with(subject.table.df2, table(traj.timepoint, season.birth)) 
with(subject.table.df2, prop.table(table(traj.timepoint, season.birth), 1)) 
with(subject.table.df2, table(traj.timepoint, season.birth2)) 
with(subject.table.df2, prop.table(table(traj.timepoint, season.birth2), 1)) 

with(subject.table.df2, chisq.test(table(traj.timepoint, season.birth))) # p = 5.1e-08
with(subject.table.df2, chisq.test(table(traj.timepoint, season.birth2))) # p = 0.019

with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.3", 
																				 season.birth))) # p = 0.012

with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.1", 
																				 season.birth == "Spring"))) # p = 0.038
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.2", 
																				 season.birth == "Autumn"))) # p = 0.046
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.3", 
																				 season.birth == "Winter"))) # p = 0.036
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.4", 
																				 season.birth == "Winter"))) # p = 0.014

# ------------------------------------------------------------------------------
# DELIVERY MODE
# ------------------------------------------------------------------------------

with(subject.table.df2, table(traj.timepoint, birthby, useNA="ifany")) 
with(subject.table.df2, prop.table(table(traj.timepoint, birthby),1)) 
with(subject.table.df2, chisq.test(table(traj.timepoint, birthby))) # p = 0.80
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.1", birthby))) # p = 0.84
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.3", birthby))) # p = 0.95
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.4", birthby))) # p = 0.44

# ------------------------------------------------------------------------------
# CLUSTERS
# ------------------------------------------------------------------------------

with(subject.table.df2, table(traj.timepoint, group, useNA="ifany")) 
with(subject.table.df2, prop.table(table(traj.timepoint, group),1)) 
with(subject.table.df2, chisq.test(table(traj.timepoint, group))) # p = 0.85
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.3", group))) # p = 0.34
with(subject.table.df2, chisq.test(table(traj.timepoint == "Traj.3", group == 2))) # p = 0.29

with(subject.table.df2, prop.table(table(group, asthma.yr13),1))

glm.results <- glm(I(asthma.yr13 == "Yes") ~ 
									 	(group == 3) * (traj.timepoint == "Traj.3"),
		data=subject.table.df2,
		family=binomial)
s <- summary(glm.results)
OR <- exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*s$coefficients[,2])
uci <- exp(s$coefficients[,1] + 1.96*s$coefficients[,2])

data.frame(OR, lci, uci) 


glm.results <- glm(I(asthma.yr13 == "Yes") ~ 
									 	group * traj.timepoint,
									 data=subject.table.df2,
									 family=binomial)
s <- summary(glm.results)
OR <- exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*s$coefficients[,2])
uci <- exp(s$coefficients[,1] + 1.96*s$coefficients[,2])

data.frame(OR, lci, uci) 

glm.results <- glm(I(asthma.yr13 == "Yes") ~ 
									 	factor(group) + traj.timepoint,
									 data=subject.table.df2,
									 family=binomial)
s <- summary(glm.results)
OR <- exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*s$coefficients[,2])
uci <- exp(s$coefficients[,1] + 1.96*s$coefficients[,2])

data.frame(OR, lci, uci) 

save.image("output/2018-05-18_MFA_trajectories_otu_table.Rdata")
