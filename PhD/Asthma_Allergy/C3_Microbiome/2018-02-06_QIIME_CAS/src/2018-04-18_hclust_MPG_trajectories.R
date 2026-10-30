# ==============================================================================
# 2018-04-18 - Hclust MPG trajectories
# ==============================================================================

rm(list=ls())

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------

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
											"cluster",
											"gee",
											"reshape2")

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

# Install metagMisc
# ---------------------------------------
if(!("metagMisc" %in% installed.packages()[,"Package"])){
	library(devtools)
	install_github("vmikk/metagMisc")
}
library(metagMisc)

# Ggplot and graphical settings
# ---------------------------------------
theme_set(theme_bw())

seqPalette <- colorRampPalette(c("#fef0d9","#fdd49e","#fdbb84",
																 "#fc8d59","#e34a33","#b30000"))
quaPalette <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
								'#ff7f00','#ffff33','#a65628','#f781bf',
								'#999999')
# quaPalette <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4',
#								'#fed9a6','#ffffcc','#e5d8bd','#fddaec',
#								'#f2f2f2')
divPalette <- colorRampPalette(rev(c("#d73027","#fc8d59","#fee090",
																		 "#e0f3f8","#91bfdb","#4575b4")))
# Functions
# ---------------------------------------
# A better signif
good.signif <- function(x, n) {
	trimws(formatC(signif(x,digits=n), digits=n,format="fg"))
}
# Note that this converts number to character string

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("output/2018-04-03_micro_assoc_analysis.Rdata")

# ------------------------------------------------------------------------------
# All samples
# ------------------------------------------------------------------------------
subj.table.all.df <- 
	log.table.df %>% 
	filter(age.custom %in% c("2mth", "6mth", "1y", "1.5y", "2y")) %>%
	mutate(age.custom = factor(age.custom, labels=c("2mth", "6mth", "1y", "1.5y", "2y", NULL, NULL, NULL, NULL))) %>%
	group_by(age.custom, SUBJ_ID, MPG.name) %>%
	summarize(n = n())

subj.table.all.df$age.MPG <-
	with(subj.table.all.df,
			 interaction(age.custom, MPG.name))

# x <- levels(subj.table.all.df$age.MPG)[
# 	order(levels(subj.table.all.df$age.MPG))]
# 
# cat(paste(shQuote(x, type="cmd"), collapse=", "))

age.MPG.levels <- c("2mth.Haemophilus MPG2", "2mth.Streptococcus MPG1","2mth.Haemophilus MPG1",
			 "2mth.Moraxella MPG", "2mth.Mixed MPG", "2mth.Corynebacterium MPG", 
			 "2mth.Alloiococcus MPG", "2mth.Staphylococcus MPG", "2mth.Moraxellaceae MPG1", 
			 "2mth.Haemophilus rare MPG", "2mth.Moraxella rare MPG", "2mth.Moraxellaceae rare MPG",
			 "2mth.Streptococcus rare MPG", "2mth.Corynebacterium rare MPG",
			 
			 "6mth.Haemophilus MPG2", "6mth.Streptococcus MPG1","6mth.Haemophilus MPG1",
			 "6mth.Moraxella MPG", "6mth.Mixed MPG", "6mth.Corynebacterium MPG", 
			 "6mth.Alloiococcus MPG", "6mth.Staphylococcus MPG", "6mth.Moraxellaceae MPG1", 
			 "6mth.Haemophilus rare MPG", "6mth.Moraxella rare MPG", "6mth.Moraxellaceae rare MPG",
			 "6mth.Streptococcus rare MPG", "6mth.Corynebacterium rare MPG",

			 "1y.Haemophilus MPG2", "1y.Streptococcus MPG1","1y.Haemophilus MPG1",
			 "1y.Moraxella MPG", "1y.Mixed MPG", "1y.Corynebacterium MPG", 
			 "1y.Alloiococcus MPG", "1y.Staphylococcus MPG", "1y.Moraxellaceae MPG1", 
			 "1y.Haemophilus rare MPG", "1y.Moraxella rare MPG", "1y.Moraxellaceae rare MPG",
			 "1y.Streptococcus rare MPG", "1y.Corynebacterium rare MPG",
			 
			 "1.5y.Haemophilus MPG2", "1.5y.Streptococcus MPG1","1.5y.Haemophilus MPG1",
			 "1.5y.Moraxella MPG", "1.5y.Mixed MPG", "1.5y.Corynebacterium MPG", 
			 "1.5y.Alloiococcus MPG", "1.5y.Staphylococcus MPG", "1.5y.Moraxellaceae MPG1", 
			 "1.5y.Haemophilus rare MPG", "1.5y.Moraxella rare MPG", "1.5y.Moraxellaceae rare MPG",
			 "1.5y.Streptococcus rare MPG", "1.5y.Corynebacterium rare MPG",
			 
			 "2y.Haemophilus MPG2", "2y.Streptococcus MPG1","2y.Haemophilus MPG1",
			 "2y.Moraxella MPG", "2y.Mixed MPG", "2y.Corynebacterium MPG", 
			 "2y.Alloiococcus MPG", "2y.Staphylococcus MPG", "2y.Moraxellaceae MPG1", 
			 "2y.Haemophilus rare MPG", "2y.Moraxella rare MPG", "2y.Moraxellaceae rare MPG",
			 "2y.Streptococcus rare MPG", "2y.Corynebacterium rare MPG")

subj.table.all.df$age.MPG <- 
	factor(subj.table.all.df$age.MPG,
				 levels = age.MPG.levels)

subj.table.all.df$subj.age <- 
	with(subj.table.all.df,
			 interaction(SUBJ_ID, age.custom))

subj.table.all.df <- 
	subj.table.all.df %>% ungroup() %>%
	group_by(age.custom, SUBJ_ID) %>%
	mutate(prop = n/sum(n))

subj.table.all.df2 <- 
	with(subj.table.all.df,
			 as.data.frame.matrix(xtabs(prop ~ SUBJ_ID + age.MPG)))

# Fill NAs
subj.table.all.df3 <- subj.table.all.df2
for (i in 1:5) {
	colstart <- 14*(i-1)+1
	colend <- 14*i
	subj.table.all.df3[which(rowSums(subj.table.all.df3[,colstart:colend]) == 0), 
										 colstart:colend] <- NA
}

subj.table.all.mat <- as.matrix(subj.table.all.df2)
subj.table.all.mat3 <- as.matrix(subj.table.all.df3)

age.MPG.colors <- c("#b29600", "#fdbf6f", "#ffd700",
			 "#e31a1c", "#ff7f00", "#a6cee3",
			 "#b2df8a", "#cab2db", "#ff00ff",
			 "#e5c100", "#fb9a99", "#ffccff", 
			 "#b15928", "#1f78b4")

MPG.colour <- rep(age.MPG.colors, 5)

# Bray Curtis preferred to Euclidean, because still dealing with "compositional" data
# (this time, MPG composition of samples within time period, rather than
# ESV composition to determine MPGs)
dist.p <- bcdist(subj.table.all.mat)

dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
# ---------------------------------------
grp <- hclust(dist.p, method = "ward")

# Get optimal number of clusters
# ---------------------------------------
# Choose number of clusters based on silhouette
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21), ylim = c(-1,1), 
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
which(med == max(med))+1 # 8
which(m.si == max(m.si))+1 # 21
k <- 4
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Cut tree to number of clusters
# ---------------------------------------
traj <- cutree(grp, k=k)
length(traj) # 244, for each subject
traj.silhouette <- si[,3]
length(traj.silhouette) # 244, for each subject

# Get heatmap
# ---------------------------------------
pdf("output/2018-04-16_hclust_mpg_traj_all.pdf", height=10, width=16)
p <- heatmap.3(subj.table.all.mat, 
					diss = FALSE, 
					Rowv = TRUE,
					Colv = FALSE, 
					dendrogram = "row", 
					dist.FUN = dist, 
					hclust.FUN = hclust,
					hclust.row = hclust(dist.p, method = ("ward")),
					kr=k,
					cluster.by.row = TRUE, 
					sepList = list(NULL, list(c(14,0,14.05,245),
																		c(28,0,28.05,245),
																		c(42,0,42.05,245),
																		c(56,0,56.05,245),
																		c(70,0,70.05,245))),
					sep.color="black",
					trace="none", 
					density.info="none",
					ColIndividualColors = as.character(MPG.colour), 
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					# color.FUN = seqPalette,
					# labCol = NULL, 
					keysize = 1, mapsize = 9, 
					main = ""
)
dev.off()
mat.row.order <- rev(row.names(p$x))

pdf("output/2018-04-16_hclust_mpg_traj_all_NA.pdf", height=10, width=16)
heatmap.3(subj.table.all.mat3[match(mat.row.order,row.names(subj.table.all.mat3)),], 
					diss = FALSE, 
					Rowv = FALSE,
					Colv = FALSE, 
					dendrogram = "none", 
					sepList = list(NULL, list(c(14,0,14.05,245),
																		c(28,0,28.05,245),
																		c(42,0,42.05,245),
																		c(56,0,56.05,245),
																		c(70,0,70.05,245))),
					sep.color="black",
					trace="none", 
					density.info="none",
					ColIndividualColors = as.character(MPG.colour), 
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					# color.FUN = seqPalette,
					# labCol = NULL, 
					keysize = 1, mapsize = 9, 
					main = ""
					)
dev.off()

# ------------------------------------------------------------------------------
# Routine healthy samples
# ------------------------------------------------------------------------------
subj.table.healthy.df <- 
	log.table.df %>% 
	filter(age.custom %in% c("2mth", "6mth", "1y", "1.5y", "2y")) %>%
	mutate(age.custom = factor(age.custom, labels=c("2mth", "6mth", "1y", "1.5y", "2y", NULL, NULL, NULL, NULL))) %>%
	filter(Category == "control",
				 uri.lri == "Well") %>%
	group_by(age.custom, SUBJ_ID, MPG.name) %>%
	summarize(n = n())

subj.table.healthy.df$age.MPG <-
	with(subj.table.healthy.df,
			 interaction(age.custom, MPG.name))

subj.table.healthy.df$age.MPG <- 
	factor(subj.table.healthy.df$age.MPG,
				 levels = age.MPG.levels)

subj.table.healthy.df <- 
	subj.table.healthy.df %>% ungroup() %>%
	group_by(age.custom, SUBJ_ID) %>%
	mutate(prop = n/sum(n))

subj.table.healthy.df2 <- 
	with(subj.table.healthy.df,
			 as.data.frame.matrix(xtabs(prop ~ SUBJ_ID + age.MPG)))

subj.table.healthy.df2 <- subj.table.healthy.df2[
	which(rowSums(subj.table.healthy.df2 == 0) != 70),]

# Fill NAs
subj.table.healthy.df3 <- subj.table.healthy.df2
for (i in 1:5) {
	colstart <- 14*(i-1)+1
	colend <- 14*i
	subj.table.healthy.df3[which(rowSums(subj.table.healthy.df3[,colstart:colend]) == 0), 
												 colstart:colend] <- NA
}

subj.table.healthy.mat <- as.matrix(subj.table.healthy.df2)
subj.table.healthy.mat3 <- as.matrix(subj.table.healthy.df3)

# Missing data stats
# Per column
colSums(is.na(subj.table.healthy.df3))[c(1,15,29,43,57)]
(colSums(is.na(subj.table.healthy.df3))/nrow(subj.table.healthy.df3))[c(1,15,29,43,57)]

# Per row
rowSums(is.na(subj.table.healthy.df3))
table(rowSums(is.na(subj.table.healthy.df3))/ncol(subj.table.healthy.df3)*5)
table(rowSums(is.na(subj.table.healthy.df3))/ncol(subj.table.healthy.df3)*5)/244

MPG.colour <- rep(age.MPG.colors, 5)

dist.p <- bcdist(subj.table.healthy.mat)

dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
# ---------------------------------------
grp <- hclust(dist.p, method = "ward")

# Get optimal number of clusters
# ---------------------------------------
# Choose number of clusters based on silhouette
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21), ylim = c(-1,1), 
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
which(m.si == max(m.si))+1 # 21
k <- 6
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Cut tree to number of clusters
# ---------------------------------------
traj <- cutree(grp, k=k)
length(traj) # 219, for each subject
traj.silhouette <- si[,3]
length(traj.silhouette) # 219, for each sample

traj.healthy <- traj
grp.healthy <- grp
dist.p.healthy <- dist.p
dist.mat.healthy <- dist.mat

# Get heatmap
# ---------------------------------------
pdf("output/2018-04-16_hclust_mpg_traj_healthy.pdf", height=10, width=16)
p <- heatmap.3(subj.table.healthy.mat, 
							 diss = FALSE, 
							 Rowv = TRUE,
							 Colv = FALSE, 
							 dendrogram = "row", 
							 dist.FUN = dist, 
							 hclust.FUN = hclust,
							 hclust.row = hclust(dist.p, method = ("ward")),
							 kr=k,
							 cluster.by.row = TRUE, 
							 sepList = list(NULL, list(c(14,0,14.05,220),
							 													c(28,0,28.05,220),
							 													c(42,0,42.05,220),
							 													c(56,0,56.05,220),
							 													c(70,0,70.05,220))),
							 sep.color="white",
							 trace="none", 
							 density.info="none",
							 ColIndividualColors = as.character(MPG.colour), 
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()
mat.row.order <- rev(row.names(p$x))

pdf("output/2018-04-16_hclust_mpg_traj_healthy_NA.pdf", height=10, width=16)
heatmap.3(subj.table.healthy.mat3[match(mat.row.order,
																				row.names(subj.table.healthy.mat3)),], 
					diss = FALSE, 
					Rowv = FALSE,
					Colv = FALSE, 
					dendrogram = "none", 
					sepList = list(NULL, list(c(14,0,14.05,220),
																		c(28,0,28.05,220),
																		c(42,0,42.05,220),
																		c(56,0,56.05,220),
																		c(70,0,70.05,220))),
					sep.color="black",
					trace="none", 
					density.info="none",
					ColIndividualColors = as.character(MPG.colour), 
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					# color.FUN = seqPalette,
					# labCol = NULL, 
					keysize = 1, mapsize = 9, 
					main = ""
)
dev.off()

# ------------------------------------------------------------------------------
# Surveillance samples (healthy or not)
# ------------------------------------------------------------------------------

subj.table.routine.df <- 
	log.table.df %>% 
	filter(Category == "control") %>%
	group_by(age.custom, SUBJ_ID, MPG.name) %>%
	summarize(n = n())

subj.table.routine.df$age.MPG <-
	with(subj.table.routine.df,
		 interaction(age.custom, MPG.name))

subj.table.routine.df$age.MPG <- 
	factor(subj.table.routine.df$age.MPG,
						levels = age.MPG.levels)

subj.table.routine.df <- 
	subj.table.routine.df %>% ungroup() %>%
	group_by(age.custom, SUBJ_ID) %>%
	mutate(prop = n/sum(n))

subj.table.routine.df2 <- 
	with(subj.table.routine.df,
			 as.data.frame.matrix(xtabs(prop ~ SUBJ_ID + age.MPG)))

# Fill NAs
subj.table.routine.df3 <- subj.table.routine.df2
for (i in 1:7) {
	colstart <- 12*(i-1)+1
	colend <- 12*i
	subj.table.routine.df3[which(rowSums(subj.table.routine.df3[,colstart:colend]) == 0), 
										 colstart:colend] <- NA
}

subj.table.routine.mat <- as.matrix(subj.table.routine.df2)
subj.table.routine.mat3 <- as.matrix(subj.table.routine.df3)

# Missing data stats
# Per column
colSums(is.na(subj.table.routine.df3))[c(1,13,25,37,49,61,73)]
(colSums(is.na(subj.table.routine.df3))/nrow(subj.table.routine.df3))[
	c(1,13,25,37,49,61,73)]

# Per row
rowSums(is.na(subj.table.routine.df3))
table(rowSums(is.na(subj.table.routine.df3))/ncol(subj.table.routine.df3)*7)

MPG.colour <- rep(age.MPG.colors, 7)

dist.p <- bcdist(subj.table.routine.mat)

dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
# ---------------------------------------
grp <- hclust(dist.p, method = "ward")

# Get optimal number of clusters
# ---------------------------------------
# Choose number of clusters based on silhouette
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21), ylim = c(-1,1), 
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
which(med == max(med))+1 # 3
which(m.si == max(m.si))+1 # 3
k <- which(med == max(med))+1
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Cut tree to number of clusters
# ---------------------------------------
traj <- cutree(grp, k=k)
length(traj) # 282, for each subject
traj.silhouette <- si[,3]
length(traj.silhouette) # 282, for each sample

traj.routine <- traj
grp.routine <- grp
dist.p.routine <- dist.p
dist.mat.routine <- dist.mat

# Get heatmap
# ---------------------------------------
pdf("output/2018-04-16_hclust_mpg_traj_routine.pdf", height=10, width=16)
p <- heatmap.3(subj.table.routine.mat, 
					diss = FALSE, 
					Rowv = TRUE,
					Colv = FALSE, 
					dendrogram = "row", 
					dist.FUN = dist, 
					hclust.FUN = hclust,
					hclust.row = hclust(dist.p, method = ("ward")),
					kr=k,
					cluster.by.row = TRUE, 
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
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					# color.FUN = seqPalette,
					# labCol = NULL, 
					keysize = 1, mapsize = 9, 
					main = ""
)
dev.off()
mat.row.order <- rev(row.names(p$x))

pdf("output/2018-04-16_hclust_mpg_traj_routine_NA.pdf", height=10, width=16)
heatmap.3(subj.table.routine.mat3[match(mat.row.order,
																				row.names(subj.table.routine.mat3)),], 
					diss = FALSE, 
					Rowv = FALSE,
					Colv = FALSE, 
					dendrogram = "none", 
					sepList = list(NULL, list(c(12,0,12.05,282),
																		c(24,0,24.05,282),
																		c(36,0,36.05,282),
																		c(48,0,48.05,282),
																		c(60,0,60.05,282),
																		c(72,0,72.05,282),
																		c(84,0,84.05,282))),
					sep.color="black",
					trace="none", 
					density.info="none",
					ColIndividualColors = as.character(MPG.colour), 
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					# color.FUN = seqPalette,
					# labCol = NULL, 
					keysize = 1, mapsize = 9, 
					main = ""
)
dev.off()

# ------------------------------------------------------------------------------
# Surveillance samples, but remove 1.5 years 
# ------------------------------------------------------------------------------

subj.table.routine.no.1.5.df2 <- subj.table.routine.df2[,-c(61:72)]
subj.table.routine.no.1.5.df3 <- subj.table.routine.df3[,-c(61:72)]

subj.table.routine.mat <- as.matrix(subj.table.routine.no.1.5.df2)
subj.table.routine.mat3 <- as.matrix(subj.table.routine.no.1.5.df3)

MPG.colour <- rep(age.MPG.colors, 7)

dist.p <- bcdist(subj.table.routine.mat)

dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
# ---------------------------------------
grp <- hclust(dist.p, method = "ward")

# Get optimal number of clusters
# ---------------------------------------
# Choose number of clusters based on silhouette
b <- list()
# Define empty vectors for filling by number of clusters
med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# Cut tree by number of clusters
# Begin with number of clusters = 2
b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
						 at = 1, xlim = c(0,21), ylim = c(-1,1), 
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
which(med == max(med))+1 # 3
which(m.si == max(m.si))+1 # 3
k <- which(med == max(med))+1
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
dev.off()

# Cut tree to number of clusters
# ---------------------------------------
traj <- cutree(grp, k=k)
length(traj) # 282, for each subject
traj.silhouette <- si[,3]
length(traj.silhouette) # 282, for each sample

traj.routine <- traj
grp.routine <- grp
dist.p.routine <- dist.p
dist.mat.routine <- dist.mat

# Get heatmap
# ---------------------------------------
pdf("output/2018-04-16_hclust_mpg_traj_routine_no_1.5.pdf", height=10, width=16)
p2 <- heatmap.3(subj.table.routine.mat, 
							 diss = FALSE, 
							 Rowv = TRUE,
							 Colv = FALSE, 
							 dendrogram = "row", 
							 dist.FUN = dist, 
							 hclust.FUN = hclust,
							 hclust.row = hclust(dist.p, method = ("ward")),
							 kr=k,
							 cluster.by.row = TRUE, 
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
							 color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
							 # color.FUN = seqPalette,
							 # labCol = NULL, 
							 keysize = 1, mapsize = 9, 
							 main = ""
)
dev.off()
mat.row.order <- rev(row.names(p2$x))

pdf("output/2018-04-16_hclust_mpg_traj_routine_no_1.5_NA.pdf", height=10, width=16)
heatmap.3(subj.table.routine.mat3[match(mat.row.order,
																				row.names(subj.table.routine.mat3)),], 
					diss = FALSE, 
					Rowv = FALSE,
					Colv = FALSE, 
					dendrogram = "none", 
					sepList = list(NULL, list(c(12,0,12.05,282),
																		c(24,0,24.05,282),
																		c(36,0,36.05,282),
																		c(48,0,48.05,282),
																		c(60,0,60.05,282),
																		c(72,0,72.05,282),
																		c(84,0,84.05,282))),
					sep.color="black",
					trace="none", 
					density.info="none",
					ColIndividualColors = as.character(MPG.colour[-c(61:72)]), 
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					# color.FUN = seqPalette,
					# labCol = NULL, 
					keysize = 1, mapsize = 9, 
					main = ""
)
dev.off()

# Compare with and without 1.5 years

df1 <- data.frame(subj = names(p$row.clusters), 
					 traj = p$row.clusters)
df2 <- data.frame(subj = names(p2$row.clusters), 
									traj.no.1.5 = p2$row.clusters)
df3 <- merge(df1, df2, by="subj", all=T)
with(df3, table(traj, traj.no.1.5, useNA="ifany"))

# Accuracy = 83%	
(103+82+50)/sum(with(df3, table(traj, traj.no.1.5, useNA="ifany")))

# ------------------------------------------------------------------------------
# Plots and GEE for healthy trajectory groups
# ------------------------------------------------------------------------------

traj.routine.df <- data.frame(subj=names(traj.routine), 
															traj.routine=factor(traj.routine,
																									levels=c(1,2,3),
																									labels=c("Traj.1", "Traj.2", "Traj.3")), 
															stringsAsFactors=F)

subject.table.df2 <- merge(subject.table.df, traj.routine.df, by="subj", all=T)
with(subject.table.df2, 
		 table(traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 table(traj.routine, group, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(traj.routine, group, useNA="ifany")))

df <- merge(log.table.df, traj.routine.df, by="subj", all=T)

pdf("output/2018-04-16_mpg_traj_uri.lri.pdf", width=8, height=5)
df %>% 
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)) %>% 
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(df$MPG.colour.hex.top)) + 
	facet_grid(traj.routine~uri.lri) +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Proportion") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) 
dev.off()

# Healthy trajectory groups vs. no. of LRI/URI samples
# -----------------------------------------

df2 <-df %>%
	filter(!is.na(age.custom), !is.na(well.sick), !is.na(traj.routine)) %>% 
	group_by(SUBJ_ID, well.sick, traj.routine) %>%
	summarize(n = n())
# Average number of well, sick samples per subject in Traj, across first 2 years
df3 <- df2 %>%
	ungroup() %>%
	group_by(traj.routine, well.sick) %>%
	summarize(mean = mean(n),
						sd = sd(n),
						median = median(n, na.rm=T),
						Q1 = quantile(n, na.rm=T)[2],
						Q3 = quantile(n, na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.routine,
											 well.sick,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), "(",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))
														
write.csv(df4, "output/2018-04-16_traj_sick_sample_count.csv", 
					quote=F, row.names=F)

# Subject-based analysis
with(df2[which(df2$well.sick == "Sick"),], kruskal.test(x=n, 
																												g=traj.routine)) # p = 0.96
with(df2[which(df2$well.sick == "Sick"),], kruskal.test(x=n, 
																												g=traj.routine == "Traj.1")) # p = 0.79

# Average number of well, uri, lri samples per subject in Traj, across first 2 years
df2 <-df %>%
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)) %>% 
	group_by(SUBJ_ID, uri.lri, traj.routine) %>%
	summarize(n = n())
df3 <- df2 %>%
	ungroup() %>%
	group_by(traj.routine, uri.lri) %>%
	summarize(mean = mean(n),
						sd = sd(n),
						median = median(n, na.rm=T),
						Q1 = quantile(n, na.rm=T)[2],
						Q3 = quantile(n, na.rm=T)[4])

write.csv(df3, "output/2018-04-16_traj_uri.lri_sample_count.csv", 
						quote=F, row.names=F)

# Subject-based analysis
with(df2[which(df2$uri.lri == "Well"),], kruskal.test(x=n, 
																											g=traj.routine)) # p = 0.007
with(df2[which(df2$uri.lri == "URI"),], kruskal.test(x=n, 
																											g=traj.routine)) # p = 0.59
with(df2[which(df2$uri.lri == "LRI"),], kruskal.test(x=n, 
																										 g=traj.routine)) # p = 0.048

with(df %>%
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)),
	table(uri.lri, traj.routine))
with(df %>%
		 	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)),
		 prop.table(table(uri.lri, traj.routine), 2))
with(df %>%
		 	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)),
		 chisq.test(table(uri.lri, traj.routine))) # p = 0.011
with(df %>%
		 	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)),
		 chisq.test(table(uri.lri, traj.routine == "Traj.1"))) # p = 0.0015

with(df %>%
		 	filter(!is.na(age.custom), !is.na(well.sick), !is.na(traj.routine)),
		 table(well.sick, traj.routine))
with(df %>%
		 	filter(!is.na(age.custom), !is.na(well.sick), !is.na(traj.routine)),
		 prop.table(table(well.sick, traj.routine), 2))
with(df %>%
		 	filter(!is.na(age.custom), !is.na(well.sick), !is.na(traj.routine)),
		 chisq.test(table(well.sick, traj.routine))) # p = 0.40
with(df %>%
		 	filter(!is.na(age.custom), !is.na(well.sick), !is.na(traj.routine)),
		 chisq.test(table(well.sick, traj.routine == "Traj.1"))) # p = 0.19

with(df, table(well.sick, traj.routine, useNA="ifany"))
with(df, prop.table(table(well.sick, traj.routine, useNA="ifany"),2))

pdf("output/2018-04-16_mpg_traj_uri.lri.freq.pdf", width=8, height=5)
df %>% 
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)) %>% 
	ggplot(aes(x=age.custom, group=uri.lri, fill=uri.lri)) + 
	geom_bar() +
	facet_grid(traj.routine ~ .) + 
	labs(x = "Age or Timepoint", y = "Count", fill = "Sample type") + 
	scale_fill_manual(values=c("Blue", "Gold", "Red"))
dev.off()

df %>% 
	filter(!is.na(age.custom), !is.na(uri.lri), !is.na(traj.routine)) %>% 
	ggplot(aes(x=age.custom, group=uri.lri, fill=uri.lri)) + 
	geom_bar(position="fill") +
	facet_grid(traj.routine ~ .) + 
	labs(x = "Age or Timepoint", y = "Count", fill = "Sample type") + 
	scale_fill_manual(values=c("Blue", "Gold", "Red"))

# Healthy trajectory groups vs. pet ownership
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(dog0, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(dog0, traj.routine, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(dog0, traj.routine))) # p=0.90
with(subject.table.df2, 
		 chisq.test(table(dog0, traj.routine == "Traj.1"))) # p=0.87

with(subject.table.df2, 
		 table(cat0, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(cat0, traj.routine, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(cat0, traj.routine))) # p=0.11
with(subject.table.df2, 
		 chisq.test(table(cat0, traj.routine == "Traj.1"))) # p=0.11

# Sample-based GEE analysis - dog ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
								 	I(dog0) + gender + age.day + uri.lri,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
								 corstr = "unstructured"), 
								 error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
dog.results.df <- results.df

# Sample-based GEE analysis - cat ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														I(cat0) + gender + age.day + uri.lri,
													id = SUBJ_ID, 
													data = log.table.df, 
													family= binomial,
													corstr = "unstructured"), 
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
cat.results.df <- results.df

# Healthy trajectory groups vs. day care
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(dayc, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(dayc, traj.routine,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.routine, dayc))) # p=0.63
with(subject.table.df2, 
		 chisq.test(table(traj.routine == "Traj.1", dayc))) # p=0.48

# Sample-based GEE analysis - dayc ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														I(dayc) + gender + age.day + uri.lri,
													id = SUBJ_ID, 
													data = log.table.df, 
													family= binomial,
													corstr = "unstructured"), 
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
dayc.results.df <- results.df

# Healthy trajectory groups vs. number of older siblings
# -----------------------------------------
with(subject.table.df2, hist(livsib))

subject.table.df2 %>%
	group_by(traj.routine) %>%
	summarize(mean = mean(livsib, na.rm=T),
						sd = sd(livsib, na.rm=T),
						median = median(livsib, na.rm=T),
						Q1 = quantile(livsib, na.rm=T)[2],
						Q3 = quantile(livsib, na.rm=T)[4])
		 
# Subject-based analysis
with(subject.table.df2, kruskal.test(x=livsib, g=traj.routine)) # p = 0.45

# Sample-based GEE analysis - livsib ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <-  tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
													 	I(livsib) + gender + age.day + uri.lri,
													 id = SUBJ_ID, 
													 data = log.table.df, 
													 family= binomial,
													 corstr = "unstructured"), 
											 error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
livsib.results.df <- results.df

# Healthy trajectory groups vs. breastfeeding
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(brfeed, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(brfeed, traj.routine,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.routine, brfeed))) # p=0.63
with(subject.table.df2, 
		 chisq.test(table(traj.routine == "Traj.1", brfeed))) # p=0.48

# Sample-based GEE analysis - brfeed ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														I(brfeed) + gender + age.day + uri.lri,
													id = SUBJ_ID, 
													data = log.table.df, 
													family= binomial,
													corstr = "unstructured"), 
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
brfeed.results.df <- results.df

# Healthy trajectory groups vs. aeroallergen sensitisation
# -----------------------------------------

# Subject-based analysis
with(subject.table.df2, 
		 table(aero.pattern, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(aero.pattern, traj.routine,  useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(aero.pattern, traj.routine))) # p=0.58
with(subject.table.df2, 
		 chisq.test(table(aero.pattern, traj.routine == "Traj.1"))) # p=0.42

with(subject.table.df2, 
		 kruskal.test(
		 	x=as.numeric(factor(aero.pattern, levels=c("C.0by5", "B.0by1.1+by5", "A.1+by1"))), 
		 	g=traj.routine))
with(subject.table.df2, 
		 kruskal.test(
		 	x=as.numeric(factor(aero.pattern, levels=c("C.0by5", "B.0by1.1+by5", "A.1+by1"))), 
		 	g=traj.routine == "Traj.1"))

# Healthy trajectory groups vs. total IgE sensitisation at age 3
# -----------------------------------------
with(subject.table.df2, hist(totige))
with(subject.table.df2, hist(log10(totige)))

df3 <- subject.table.df2 %>%
	group_by(traj.routine) %>%
	summarize(mean = mean(log10(totige), na.rm=T),
						sd = sd(log10(totige), na.rm=T),
						median = median(log10(totige), na.rm=T),
						Q1 = quantile(log10(totige), na.rm=T)[2],
						Q3 = quantile(log10(totige), na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.routine,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), " (",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=log10(totige), g=traj.routine)) # p = 0.12
with(subject.table.df2, kruskal.test(x=log10(totige), g=traj.routine == "Traj.1")) # p = 0.05

pdf("output/2018-04-16_mpg_traj_totige.pdf", width=8, height=5)
subject.table.df2 %>%
	filter(!is.na(traj.routine)) %>%
	ggplot(aes(x = traj.routine, y=log10(totige))) + geom_boxplot() + 
	labs(x="Trajectory", y="Log10(total IgE)")
dev.off()


# Sample-based GEE analysis - log10totige ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														log10(totige) + gender + age.day + uri.lri,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= gaussian,
								 corstr = "unstructured"),
								 error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
totige.results.df <- results.df


# Healthy trajectory groups vs. dpter IgE sensitisation at age 3
# -----------------------------------------
with(subject.table.df2, hist(dpter))
with(subject.table.df2, hist(log10(dpter)))

subject.table.df2 %>%
	group_by(traj.routine) %>%
	summarize(mean = mean(log10(dpter), na.rm=T),
						sd = sd(log10(dpter), na.rm=T),
						median = median(log10(dpter), na.rm=T),
						Q1 = quantile(log10(dpter), na.rm=T)[2],
						Q3 = quantile(log10(dpter), na.rm=T)[4])

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=log10(dpter), g=traj.routine)) # p = 0.23

# Healthy trajectory groups vs. cat IgE sensitisation at age 3
# -----------------------------------------
with(subject.table.df2, hist(cat))
with(subject.table.df2, hist(log10(cat)))

subject.table.df2 %>%
	group_by(traj.routine) %>%
	summarize(mean = mean(log10(cat), na.rm=T),
						sd = sd(log10(cat), na.rm=T),
						median = median(log10(cat), na.rm=T),
						Q1 = quantile(log10(cat), na.rm=T)[2],
						Q3 = quantile(log10(cat), na.rm=T)[4])

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=log10(cat), g=traj.routine)) # p = 0.60

# Healthy trajectory groups vs. dog IgE sensitisation at age 3
# -----------------------------------------
with(subject.table.df2, hist(dog))
with(subject.table.df2, hist(log10(dog)))

subject.table.df2 %>%
	group_by(traj.routine) %>%
	summarize(mean = mean(log10(dog), na.rm=T),
						sd = sd(log10(dog), na.rm=T),
						median = median(log10(dog), na.rm=T),
						Q1 = quantile(log10(dog), na.rm=T)[2],
						Q3 = quantile(log10(dog), na.rm=T)[4])

# Subject-based analysis
with(subject.table.df2, kruskal.test(x=log10(dog), g=traj.routine)) # p = 0.39

# Healthy trajectory groups vs. asthma
# -----------------------------------------
# Subject-based analysis
with(subject.table.df2, 
		 table(asthma.yr6, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr6, traj.routine, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(traj.routine, asthma.yr6))) # p=0.56

# Sample-based GEE analysis - asthma.yr6 ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~
								 	I(asthma.yr6 == "Yes")+ gender + age.day + uri.lri,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
								 corstr = "unstructured"),
								 error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}								 
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
asthma.yr6.results.df <- results.df

# Subject-based analysis
with(subject.table.df2, 
		 table(asthma.yr8, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr8, traj.routine, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(asthma.yr8, traj.routine))) # p=0.072

# Sample-based GEE analysis - asthma.yr8 ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~
														I(asthma.yr8 == "Yes")+ gender + age.day + uri.lri,
													id = SUBJ_ID, 
													data = log.table.df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}								 
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
asthma.yr8.results.df <- results.df

# Subject-based analysis
with(subject.table.df2, 
		 table(asthma.yr11, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr11, traj.routine, useNA="ifany"),2))
with(subject.table.df2, 
		 chisq.test(table(asthma.yr11, traj.routine))) # p=0.108

# Sample-based GEE analysis - asthma.yr11 ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~
														I(asthma.yr11 == "Yes")+ gender + age.day + uri.lri,
													id = SUBJ_ID, 
													data = log.table.df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}								 
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
asthma.yr11.results.df <- results.df


# Subject-based analysis
with(subject.table.df2, 
		 table(asthma.yr13, traj.routine, useNA="ifany"))
with(subject.table.df2, 
		 prop.table(table(asthma.yr13, traj.routine,  useNA="ifany"), 2))
with(subject.table.df2, 
		 chisq.test(table( asthma.yr13, traj.routine))) # p=0.04
with(subject.table.df2, 
		 chisq.test(table( asthma.yr13, traj.routine == "Traj.3"))) # p=0.022

# Sample-based GEE analysis - asthma.yr13 ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														I(asthma.yr13 == "Yes") + gender + age.day + uri.lri,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
								 corstr = "unstructured"),
								 error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
asthma.yr13.results.df <- results.df

# Healthy trajectory groups vs. lung function at age 6
# -----------------------------------------
lung.function.prebd.df <- read.csv("data/melbourne_spirometry_preBD.csv")
lung.function.prebd.df.6y <- lung.function.prebd.df %>%
	filter(whatyear==6)
subject.table.df3 <- merge(subject.table.df2, lung.function.prebd.df.6y, 
													 by="subj", all=T)

with(subject.table.df3, hist(FEV1FVC_percent_predicted/100))

df3 <- subject.table.df3 %>%
	group_by(traj.routine) %>%
	summarize(mean = mean(FEV1FVC_percent_predicted/100, na.rm=T),
						sd = sd(FEV1FVC_percent_predicted/100, na.rm=T),
						median = median(FEV1FVC_percent_predicted/100, na.rm=T),
						Q1 = quantile(FEV1FVC_percent_predicted/100, na.rm=T)[2],
						Q3 = quantile(FEV1FVC_percent_predicted/100, na.rm=T)[4])
df4 <- with(df3, 
						data.frame(traj.routine,
											 mean.sd = paste0(good.signif(mean, 2), " (", good.signif(sd, 2), ")"),
											 median.iqr = paste0(good.signif(median, 2), " (",
											 										good.signif(Q1, 2), "-",
											 										good.signif(Q3, 2), ")")))

# Subject-based analysis
with(subject.table.df3, kruskal.test(x=FEV1FVC_percent_predicted/100, g=traj.routine)) # p = 0.14
with(subject.table.df3, kruskal.test(x=FEV1FVC_percent_predicted/100, g=traj.routine == "Traj.1")) # p = 0.048

pdf("output/2018-04-16_mpg_traj_fev1fvcpred.6y.pdf", width=8, height=5)
subject.table.df3 %>%
	filter(!is.na(traj.routine)) %>%
	ggplot(aes(x = traj.routine, y=FEV1FVC_percent_predicted/100)) + geom_boxplot() + 
	labs(x="Trajectory", y="FEV1:FVC % predicted, at age 6")
dev.off()


# Sample-based GEE analysis - log10totige ~ MPG
results.df <- data.frame(MPG = MPG.list, 
												 OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														FEV1FVC_percent_predicted/100 + gender + age.day + uri.lri,
													id = SUBJ_ID, 
													data = log.table.df, 
													family= gaussian,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

results.df <- results.df[match(MPG.sick.sig.order, results.df$MPG),]
totige.results.df <- results.df

save.image("output/2018-04-18_hclust_MPG_trajectories.Rdata")