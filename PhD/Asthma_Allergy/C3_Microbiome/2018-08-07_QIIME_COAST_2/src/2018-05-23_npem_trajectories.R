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
											"mixtools",
											"cluster")

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

otu.table.npem.list <- list()

j <- seq(1:20) # Number of components
k <- numeric(20) # Log likelihood measurements

subj.table.routine.mat4.imp.s <- scale(subj.table.routine.mat4.imp)

# Starting values strongly influence clustering outcomes
# i = number of initial centres
for (i in 1:20){
	set.seed(1)
	otu.table.npem.list[[i]] <- npEM(subj.table.routine.mat4.imp.s, 
																	 mu0=i, 
																	 stochastic=F)
	k[i] <- otu.table.npem.list[[i]]$loglik[length(otu.table.npem.list[[i]]$loglik)]
}
# BIC = -2*loglik + df * log(n)
# df approximately = number of components * number of variables (for mixtures)
# + (components - 1) (for mixing proportions)?
BIC <- -2*k + (j*ncol(subj.table.routine.mat4.imp.s)+j-1)*log(nrow(subj.table.routine.mat4.imp.s))

pdf("output/2018-05-23_BIC_npem.pdf")
plot(j, BIC, type="b",
		 xlab="Number of clusters",
		 ylab="BIC",
		 main="BIC vs. number of npEM clusters")
dev.off()

which(BIC==min(BIC))

# Group assignment and silhouette widths
# -----------------------------------
groups.otu.table.npem.list <- list()
otu.table.npem.groups.list <- list()
silhouette.npem.groups <- list()
for (i in 1:20){
	groups.otu.table.npem.list[[i]] <- data.frame(rownames(subj.table.routine.mat4.imp.s),
																								otu.table.npem.list[[i]]$posteriors>0.5)
	colnames(groups.otu.table.npem.list[[i]])[1] <- c("subj")
	
	groups.otu.table.npem.list[[i]] <- melt(groups.otu.table.npem.list[[i]], 
																					id.vars="subj",
																				 variable.name = "group")
	groups.otu.table.npem.list[[i]] <- groups.otu.table.npem.list[[i]][which(
		groups.otu.table.npem.list[[i]]$value==T),]
	groups.otu.table.npem.list[[i]]$value <- NULL
	groups.otu.table.npem.list[[i]]$group <- factor(groups.otu.table.npem.list[[i]]$group,
																								 labels=seq(1:i))
	df <- data.frame(subj.table.routine.mat4.imp.s)
	df$subj <- rownames(df)
	otu.table.npem.groups.list[[i]] <- merge(groups.otu.table.npem.list[[i]], df, 
																					by="subj")
	rownames(otu.table.npem.groups.list[[i]]) <- otu.table.npem.groups.list[[i]]$subj
	otu.table.npem.groups.mat <- otu.table.npem.groups.list[[i]][,which(!(
		colnames(otu.table.npem.groups.list[[i]]) %in% c("subj", "group")))]
	silhouette.npem.groups[[i]] <- silhouette(x = as.numeric(otu.table.npem.groups.list[[i]]$group), 
																						dist = dist(otu.table.npem.groups.mat, method="euclidean"))
}

# Mean silhouette width
mean.silhouette <- numeric(20)
for (i in 2:20){
	mean.silhouette[i] <- (summary(silhouette.npem.groups[[i]])$si.summary)[4]
}

pdf("output/2018-05-23_mean_silhouette_npem.pdf")
plot(j, mean.silhouette, type="b",
		 xlab="Number of clusters",
		 ylab="Mean silhouette width",
		 main="Mean silhouette width vs. number of npEM clusters")
dev.off()

which(mean.silhouette==max(mean.silhouette))

# 11 clusters?

#' Non-parametric EM clustering of scaled dataset (k=3)
#' -----------------------------------------------------------------------------
k <- 3
otu.table.npem <- otu.table.npem.list[[k]]

otu.table.npem.group.df <- data.frame(subj = rownames(subj.table.routine.mat4.imp),
																		 npem.group = otu.table.npem$posteriors>0.5)

colnames(otu.table.npem.group.df)[1] <- c("subj")

otu.table.npem.group.df <- melt(otu.table.npem.group.df, id.vars="subj",
																	 variable.name = "group")
otu.table.npem.group.df <- otu.table.npem.group.df[which(otu.table.npem.group.df$value==T),]
otu.table.npem.group.df$value <- NULL
otu.table.npem.group.df$group <- factor(otu.table.npem.group.df$group,
																					 labels=c(1:k))

# Heatmap with raw data
# ---------------------------------------

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
subj.table.routine.mat4.2.imp <- subj.table.routine.mat4.2
for(i in 1:ncol(subj.table.routine.mat4.2.imp)){
	subj.table.routine.mat4.2.imp[is.na(subj.table.routine.mat4.2.imp[,i]), i] <- 
		mean(subj.table.routine.mat4.2.imp[,i], na.rm = TRUE)
}

p4 <- heatmap.3(subj.table.routine.mat4.2.imp[match(otu.table.npem.group.df$subj, 
																										row.names(subj.table.routine.mat4.2.imp)),], 
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
								RowIndividualColors = as.character(col2hex(otu.table.npem.group.df$group[
									order(otu.table.npem.group.df$group)])),
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								# labCol = NULL, 
								keysize = 1, mapsize = 9, 
								main = ""
)
dev.off()

trajectories.timepoint.group <- otu.table.npem.group.df$group 

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

df4  %>% 
	mutate(traj.timepoint = factor(traj.timepoint)) %>%
	ggplot(aes(x = age, y = value, color = otu)) + 
	geom_violin() +  
	stat_summary(fun.y="mean", geom="point") +
	stat_summary(aes(group = 1), fun.y=mean, geom="line", linetype=2) +
	scale_color_manual(values = as.character(otu.colour.scheme.df$hex)) + 
	facet_grid(otu~traj.timepoint) + 
	labs(y = "Relative abundance", x = "Timepoint", group = "OTU", color = "OTU") + 
	theme(strip.text.y = element_blank())

