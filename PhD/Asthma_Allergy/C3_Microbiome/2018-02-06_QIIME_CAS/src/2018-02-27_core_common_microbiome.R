# ==============================================================================
# 2018-02-27 - Core microbiome analysis
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

load("output/2018-02-06_process_phyloseq_1.Rdata")

# ------------------------------------------------------------------------------
# Get common ESVs
# ------------------------------------------------------------------------------

# # Method 1: recommended by microbiome package, with or without css
# # ---------------------------------------
# cas.biom.method1.css.core # 13 taxa, 2926 samples
# # cas.biom.method1.css.core.rel <- microbiome::transform(cas.biom.method1.css.core, "compositional")
# # cas.biom.method1.css.core.rel <- subset_samples(cas.biom.method1.css.core.rel, 
# # 																								sample_sums(cas.biom.method1.css.core.rel) != 0)
# # cas.biom.method1.css.core.rel 
# otu1.mat <- otu_table(cas.biom.method1.css.core)@.Data
# tax1.df <- data.frame(tax_table(cas.biom.method1.css.core)@.Data)
# tax1.df$taxa <- paste0(tax1.df$Genus, " ", tax1.df$Species)
# tax1.df$taxa <- gsub(pattern="[a-z]__",replacement="", tax1.df$taxa)
# tax1.df$taxa <- gsub(pattern="NA",replacement="", tax1.df$taxa)
# tax1.df$taxa[11] <- "Gemellaceae"
# tax1.df$taxa <- trimws(tax1.df$taxa, which = c("right"))
# tax1.df$taxa <- paste0(tax1.df$taxa, ", ", row.names(tax1.df))
# tax1.df
# 
# rownames(otu1.mat) <- tax1.df$taxa
# 
# dist1.p <- bcdist(t(otu1.mat))
# 
# pdf("output/2018-02-06_core_hclust.pdf", width=16, height=10)
# h1 <- heatmap.3(otu1.mat, 
# 								diss = FALSE, 
# 								hclust.FUN = hclust,
# 								dendrogram = "col", 
# 								cluster.by.col = TRUE, 
# 								Rowv =FALSE, 
# 								trace="none", 
# 								density.info="none",
# 								# ColIndividualColors = as.character(otu.mat$MPG.colour), 
# 								# color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
# 								color.FUN = seqPalette,
# 								labCol = NULL, keysize = 2, mapsize = 5, 
# 								dist.FUN = dist, 
# 								hclust.col = hclust(dist1.p, method = "complete"), 
# 								main = ""
# 								#breaks = c(0,0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
# )
# dev.off()

# Method 2: Shumei's method, with css
# ---------------------------------------
# average abundance > 0.1% & present in >20% of samples and dominate (>50%) at least 1 sample
cas.biom.method1.css # 23441 taxa, 3120 samples
mean_per_OTU <- apply(otu_table(cas.biom.method1.css)@.Data,1,mean)
presence_per_OTU <- apply(otu_table(cas.biom.method1.css)@.Data,1,function(x){sum(x!=0)})
max_per_OTU <- apply(otu_table(cas.biom.method1.css)@.Data,1,max)
indiv.name <- row.names(otu_table(cas.biom.method1.css))

common.otus <- indiv.name[which(mean_per_OTU>0.001 & 
																	presence_per_OTU > 0.20*ncol(otu_table(cas.biom.method1.css)) & 
																	max_per_OTU > 0.5)]
length(common.otus) # 13 otus

cas.biom.method1.css.common <- prune_taxa(common.otus, cas.biom.method1.css)
cas.biom.method1.css.common # 13 taxa, 3120 samples

common.otus.across.cohorts <- c("dd2e15140d20c84712456e4c12ece08b",
																"cb50e79d177c169ec5b9d5e7d308874e",
																"d2a4add6029e32135562e2da9d4ed2e1",
																"d8009096d67cd7b4ce2853718bb17013", 
																"bc0d904a0469d05c65b270998488ee91",
																"f579ddd8d7c6a3d812ca33d72110f0d6",
																"d253ca966efcd811d057dc2aa6430774", 
																"6028872a168fbd6bc3381494da3b116f", 
																"a5a08d96bf11fdeed28fc530d7f5d8c7", 
																"c3d94028bb967baac75d5d1fabcb49c0",
																"03f429f55feb0e805ba48c8c00867cf3",
																"092523a4aed84953633a5f48bc2f87b9",
																"29eb8e14ce61e9c4549de0f9cf200c19",
																"357551f644064bb2b49614b1caab2fa5",
																"4060107ffdf52e44f9d72343b7332609",
																"a3a3b519e04bf11edbb1d82497c2b34c",
																"b069c8413c9eb478d4faa117c4372bd3",
																"be1bac8c38e726214223a5ac13184ab7",
																"fb81574d887fa68041a84122487a0e47")
cas.biom.method1.css.common.across.cohorts <- prune_taxa(common.otus.across.cohorts, 
																												 cas.biom.method1.css)
df <- data.frame(rowMeans(otu_table(cas.biom.method1.css.common.across.cohorts)@.Data))
df$esv <- row.names(df)
colnames(df) <- c("mean.rel.abund.cas", "esv")
write.table(df, "output/common.otus.rel.abund.txt", row.names=F, quote=F, col.names=T)

# Extract rel. abundance for each major ESV
# ---------------------------------------
tax_table(cas.biom.method1.css.common)

# Get moraxella spp.
cas.biom.moraxella <- subset_taxa(cas.biom.method1.css, Genus == "g__Moraxella")
cas.biom.moraxella # 502 otus
moraxella.rare <- row.names(otu_table(cas.biom.moraxella))[-which(
	row.names(otu_table(cas.biom.moraxella)) %in% 
		c("d253ca966efcd811d057dc2aa6430774"))] 
cas.biom.moraxella.rare <- prune_taxa(moraxella.rare, cas.biom.moraxella)
cas.biom.moraxella.rare # 501 otus
read.sum.moraxella.rare <- sample_sums(cas.biom.moraxella.rare)

# Moraxellaceae (excl. moraxella spp.)
cas.biom.moraxellaceae <- subset_taxa(cas.biom.method1.css, Family == "f__Moraxellaceae")
cas.biom.moraxellaceae # 988 otus
moraxellaceae.rare <- row.names(otu_table(cas.biom.moraxellaceae))[-which(
	row.names(otu_table(cas.biom.moraxellaceae)) %in% 
		c(row.names(otu_table(cas.biom.moraxella)),
			"6028872a168fbd6bc3381494da3b116f"))] 
cas.biom.moraxellaceae.rare <- prune_taxa(moraxellaceae.rare, cas.biom.moraxellaceae)
cas.biom.moraxellaceae.rare # 485 otus (988-501-1-1)
read.sum.moraxellaceae.rare <- sample_sums(cas.biom.moraxellaceae.rare)

# Streptococcus
cas.biom.streptococcus <- subset_taxa(cas.biom.method1.css, Genus == "g__Streptococcus")
cas.biom.streptococcus # 431 otus
streptococcus.rare <- row.names(otu_table(cas.biom.streptococcus))[-which(
	row.names(otu_table(cas.biom.streptococcus)) %in% 
		c("357551f644064bb2b49614b1caab2fa5", "4060107ffdf52e44f9d72343b7332609",
			"a3a3b519e04bf11edbb1d82497c2b34c"))] 
cas.biom.streptococcus.rare <- prune_taxa(streptococcus.rare, cas.biom.streptococcus)
cas.biom.streptococcus.rare # 428 otus
read.sum.streptococcus.rare <- sample_sums(cas.biom.streptococcus.rare)

# Alloiococcus
cas.biom.alloiococcus <- subset_taxa(cas.biom.method1.css, Genus == "g__Alloiococcus")
cas.biom.alloiococcus # 259 otus
alloiococcus.rare <- row.names(otu_table(cas.biom.alloiococcus))[-which(
	row.names(otu_table(cas.biom.alloiococcus)) %in% 
		c("dd2e15140d20c84712456e4c12ece08b"))] 
cas.biom.alloiococcus.rare <- prune_taxa(alloiococcus.rare, cas.biom.alloiococcus)
cas.biom.alloiococcus.rare # 258 otus
read.sum.alloiococcus.rare <- sample_sums(cas.biom.alloiococcus.rare)

# Corynebacterium
cas.biom.corynebacterium <- subset_taxa(cas.biom.method1.css, Genus == "g__Corynebacterium")
cas.biom.corynebacterium # 527 otus
corynebacterium.rare <- row.names(otu_table(cas.biom.corynebacterium))[-which(
	row.names(otu_table(cas.biom.corynebacterium)) %in% 
		c("cb50e79d177c169ec5b9d5e7d308874e"))] 
cas.biom.corynebacterium.rare <- prune_taxa(corynebacterium.rare, cas.biom.corynebacterium)
cas.biom.corynebacterium.rare # 526 otus
read.sum.corynebacterium.rare <- sample_sums(cas.biom.corynebacterium.rare)

# Staphylococcus
cas.biom.staphylococcus <- subset_taxa(cas.biom.method1.css, Genus == "g__Staphylococcus")
cas.biom.staphylococcus # 321 otus
staphylococcus.rare <- row.names(otu_table(cas.biom.staphylococcus))[-which(
	row.names(otu_table(cas.biom.staphylococcus)) %in% 
		c("29eb8e14ce61e9c4549de0f9cf200c19"))] 
cas.biom.staphylococcus.rare <- prune_taxa(staphylococcus.rare, cas.biom.staphylococcus)
cas.biom.staphylococcus.rare # 320 otus
read.sum.staphylococcus.rare <- sample_sums(cas.biom.staphylococcus.rare)

# Haemophilus
cas.biom.haemophilus <- subset_taxa(cas.biom.method1.css, Genus == "g__Haemophilus")
cas.biom.haemophilus # 350 otus
haemophilus.rare <- row.names(otu_table(cas.biom.haemophilus))[-which(
	row.names(otu_table(cas.biom.haemophilus)) %in% 
		c("f579ddd8d7c6a3d812ca33d72110f0d6", "bc0d904a0469d05c65b270998488ee91"))] 
cas.biom.haemophilus.rare <- prune_taxa(haemophilus.rare, cas.biom.haemophilus)
cas.biom.haemophilus.rare # 348 otus
read.sum.haemophilus.rare <- sample_sums(cas.biom.haemophilus.rare)

# Others
read.sum.others.rare <- 1 - sample_sums(cas.biom.method1.css.common) -
	read.sum.moraxella.rare - read.sum.moraxellaceae.rare - read.sum.streptococcus.rare -
	read.sum.alloiococcus.rare - read.sum.corynebacterium.rare - read.sum.staphylococcus.rare -
	read.sum.haemophilus.rare
read.sum.others.rare[which(read.sum.others.rare < 0)] <- 0

other.rare <- row.names(otu_table(cas.biom.method1.css))[-which(
	row.names(otu_table(cas.biom.method1.css)) %in% 
		c(row.names(otu_table(cas.biom.method1.css.common)),
				 moraxella.rare, moraxellaceae.rare,
				 streptococcus.rare, alloiococcus.rare,
				 corynebacterium.rare, staphylococcus.rare,
				 haemophilus.rare))]
cas.biom.other <- prune_taxa(other.rare, cas.biom.method1.css)
cas.biom.other # 20562 taxa, 3120 samples
top.other.rare <- top_taxa(cas.biom.other, n=10)
tax_table(cas.biom.other)[
	which(row.names(tax_table(cas.biom.other)) %in% top.other.rare),]

# Neisseria spp. (cinerea), Veillonella spp. (dispar, parvula), 
# Prevotella melaninogenica, Rothia mucilaginosa

# Create OTU/ESV table
# ---------------------------------------
otu.mat.common <- otu_table(cas.biom.method1.css.common)@.Data
otu.mat <- otu_table(cas.biom.method1.css.common)@.Data
write.csv(rownames(otu.mat.common), "output/common.otus.csv", quote=F, row.names=F, col.names=F)

tax.df <- data.frame(tax_table(cas.biom.method1.css.common)@.Data)
tax.df$taxa <- paste0(tax.df$Genus, " ", tax.df$Species)
tax.df$taxa <- gsub(pattern="[a-z]__",replacement="", tax.df$taxa)
tax.df$taxa <- gsub(pattern="NA",replacement="", tax.df$taxa)
tax.df$taxa <- trimws(tax.df$taxa, which = c("right"))
tax.df$taxa[9] <- "Moraxellaceae"
tax.df$taxa[11] <- "Gemellaceae"
tax.df$taxa <- paste0(tax.df$taxa, ", ", row.names(tax.df))
rownames(otu.mat.common) <- tax.df$taxa
otu.mat.common <- otu.mat.common[order(rownames(otu.mat.common)),]
otu.mat2 <- rbind(otu.mat.common, 
								 read.sum.alloiococcus.rare,
								 read.sum.corynebacterium.rare,
								 read.sum.haemophilus.rare,
								 read.sum.moraxella.rare,
								 read.sum.moraxellaceae.rare,
								 read.sum.staphylococcus.rare,
								 read.sum.streptococcus.rare,
								 read.sum.others.rare)
row.names(otu.mat2)[(nrow(otu.mat2)-7):nrow(otu.mat2)] <- 
	c("Alloiococcus, rare", "Corynebacterium, rare", "Haemophilus, rare", 
		"Moraxella, rare", "Moraxellaceae, rare", "Staphylococcus, rare", 
		"Streptococcus, rare", "Others, rare")

# Get distance matrix
# ---------------------------------------
dist.p <- bcdist(t(otu.mat2))

# ------------------------------------------------------------------------------
# Define MPGs
# ------------------------------------------------------------------------------

# Do Hclust
# ---------------------------------------
# MPG is defined using hclust with bray curtis distance and complete linkage
dist.mat <- as.matrix(dist.p)

# Run hclust on distance matrix from OTUs
grp <- hclust(dist.p, method = "complete")

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
which(med == max(med))+1 # 14
which(m.si == max(m.si))+1 # 9
k <- which(med == max(med))+1
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)

# Cut tree to number of clusters
# ---------------------------------------
MPG <- cutree(grp, k=k)
length(MPG) # 3120, for each sample
MPG.silhouette <- si[,3]
length(MPG.silhouette) # 3120, for each sample

# Define MPG names and colours
MPG.colour.scheme.df <- data.frame(
	hex = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
					'#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
					'#808080', '#000000'),
	descript = c("light blue", "blue", "light green", "green", "pink", "red",
							 "light orange", "orange", "light purple", "purple", "yellow", "brown",
							 "grey", "black"))

MPG.colour <- factor(MPG, levels=c(1:14),
										 labels=MPG.colour.scheme.df$hex)

h2 <- heatmap.3(otu.mat2, 
								diss = FALSE, 
								hclust.FUN = hclust,
								dendrogram = "col", 
								cluster.by.col = TRUE, 
								Rowv =FALSE, 
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(MPG.colour), 
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								labCol = NULL, keysize = 2, mapsize = 5, 
								dist.FUN = dist, 
								hclust.col = hclust(dist.p, method = "complete"), 
								main = "",
								breaks = c(0,0.001,
													 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
													 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0)
)

# Manually map colours to OTU names
# ---------------------------------------
# hex = c('#a6cee3','#1f78b4',
# 				'#b2df8a','#33a02c',
# 				'#fb9a99','#e31a1c',
# 				'#fdbf6f','#ff7f00',
# 				'#cab2d6','#6a3d9a',
# 				'#ffff99','#b15928',
# 				'#808080', '#000000')
# descript = c("light blue", "blue",
# 						 "light green", "green", 
# 						 "pink", "red",
# 						 "light orange", "orange", 
# 						 "light purple", "purple", 
# 						 "yellow", "brown",
# 						 "grey", "black")

MPG.colour.scheme.df$MPG.names <- 
	c("Moraxella d253", "Haemophilus f579", 
		"Haemophilus influenzae bc0d", "Corynebacterium cb50",
		"Others, rare/Streptococcus 3575", "Streptococcus 4060",
		"Corynebacterium, rare", "Staphylococcus 29eb", 
		"Alloiococcus dd2e/Corynebacterium cb50", "Streptococcus, rare",
		"Moraxella, rare", "Moraxellaceae 6028",
		"Haemophilus, rare", "Moraxellaceae, rare")

MPG.colour.scheme.df$hex2 <- 
	c("#e31a1c", "#ffd700",
		"#b29600", "#a6cee3",
		"#808080", "#fdbf6f",
		"#1f78b4", "#cab2db", 
		"#b2df8a", "#b15928",
		"#fb9a99", "#ff00ff",
		"#e5c100", "#ffccff")

MPG.colour.scheme.df$descript2 <- 
	c("red", "gold",
		"darker gold", "light blue",
		"grey", "light orange",
		"blue", "light purple", 
		"light green", "brown",
		"pink", "magenta",
		"dark gold", "light magenta")
	
MPG.colour <- factor(MPG, levels=c(1:14),
										 labels=MPG.colour.scheme.df$hex2)

# Finalise heatmap
# ---------------------------------------
pdf("output/2018-02-06_common_hclust.pdf", width=16, height=10)
h2 <- heatmap.3(otu.mat2, 
								diss = FALSE, 
								hclust.FUN = hclust,
								dendrogram = "col", 
								cluster.by.col = TRUE, 
								Rowv =FALSE, 
								trace="none", 
								density.info="none",
								ColIndividualColors = as.character(MPG.colour), 
								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								# color.FUN = seqPalette,
								labCol = NULL, keysize = 2, mapsize = 5, 
								dist.FUN = dist, 
								hclust.col = hclust(dist.p, method = "complete"), 
								main = "",
								breaks = c(0,0.001,
													 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
													 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0))
dev.off()

# ------------------------------------------------------------------------------
# Get top taxa 
# ------------------------------------------------------------------------------
top.taxa <- top_taxa(cas.biom.method1.css, n=10)
top.taxa.tax.df <- tax_table(cas.biom.method1.css)[
	which(row.names(tax_table(cas.biom.method1.css)) %in% top.taxa),]
top.taxa.tax.df[order(match(row.names(top.taxa.tax.df), top.taxa))]

# ------------------------------------------------------------------------------
# Integrate MPG information with metadata
# ------------------------------------------------------------------------------
MPG.colour.scheme.df$MPG <- c(1:nrow(MPG.colour.scheme.df))
MPG.colour.scheme.df <- data.frame(MPG.colour.scheme.df,
										 MPG.freq = data.frame(table(MPG))$Freq)
MPG.colour.scheme.df <- 
	MPG.colour.scheme.df[order(MPG.colour.scheme.df$MPG.freq, decreasing=T),]
MPG.colour.scheme.df$MPG.name2 <-factor(c("Moraxella MPG",
																		"Mixed MPG",
																		"Streptococcus MPG1",
																		"Staphylococcus MPG",
																		"Alloiococcus MPG",
																		"Haemophilus MPG1",
																		"Corynebacterium MPG",
																		"Moraxellaceae MPG1",
																		"Haemophilus MPG2",
																		"Moraxella rare MPG",
																		"Corynebacterium rare MPG",
																		"Haemophilus rare MPG",
																		"Moraxellaceae rare MPG",
																		"Streptococcus rare MPG"))
	
MPG.df <- data.frame(LogNo = names(MPG), 
										 MPG = MPG,
										 stringsAsFactors = F)
MPG.df$MPG.name.old <- factor(MPG.df$MPG, 
													levels = MPG.colour.scheme.df$MPG,
													labels = MPG.colour.scheme.df$MPG.names)
MPG.df$MPG.name <- factor(MPG.df$MPG, 
															levels = MPG.colour.scheme.df$MPG,
															labels = MPG.colour.scheme.df$MPG.name2)
MPG.df$MPG.colour.name <- factor(MPG.df$MPG, 
															levels = MPG.colour.scheme.df$MPG,
															labels = MPG.colour.scheme.df$descript2)
MPG.df$MPG.colour.hex <- factor(MPG.df$MPG, 
																levels = MPG.colour.scheme.df$MPG,
																labels = MPG.colour.scheme.df$hex2)
	
# Generate metadata for each log/sample
log.table.df <- data.frame(sample_data(cas.biom.method1.css))
log.table.df$LogNo_old <- log.table.df$LogNo
log.table.df$LogNo <- row.names(log.table.df)
log.table.df <- merge(log.table.df, MPG.df, by="LogNo", all=T)

# ------------------------------------------------------------------------------
# Plot dominance of MPGs by ESVs
# ------------------------------------------------------------------------------

# major.otus <- c("dd2e15140d20c84712456e4c12ece08b",
# 								"cb50e79d177c169ec5b9d5e7d308874e",
# 								"f579ddd8d7c6a3d812ca33d72110f0d6",
# 								"bc0d904a0469d05c65b270998488ee91",
# 								"d253ca966efcd811d057dc2aa6430774",
# 								"4060107ffdf52e44f9d72343b7332609",
# 								"357551f644064bb2b49614b1caab2fa5",
# 								"29eb8e14ce61e9c4549de0f9cf200c19")
# cas.biom.method1.css.major <- prune_taxa(major.otus, 
# 																					 cas.biom.method1.css.common)

otu.tax.col <- t(otu.mat2)
otu.tax.col <- otu.tax.col[match(log.table.df$LogNo, row.names(otu.tax.col)),]
log.esv.df <- data.frame(LogNo = log.table.df$LogNo, 
												 MPG.name = log.table.df$MPG.name, 
												 otu.tax.col)

log.esv.df.m <- log.esv.df %>%
	melt(id.vars = c("LogNo", "MPG.name"),
			 variable.name = "ESV",
			 value.name = "Rel.Abundance") %>%
	mutate(ESV = factor(ESV),
				 ESV.simp = factor(as.numeric(ESV)))

# Colors to match ESVs
levels(log.esv.df.m$ESV)
MPG.colour.scheme.df
ESV.colour.hex <- c("#b2df8a", "#a6cee3", "grey", 
										"grey", "#b29600", "#ffd700",
										"#e31a1c", "#ff00ff", "grey",
										"#cab2db", "#ff7f00", "#fdbf6f", 
										"grey", "grey", "#1f78b4", 
										"#e5c100", "#fb9a99", "#ffccff", 
										"grey", "#b15928", "grey")

# Custom labels
levels(log.esv.df.m$ESV)
ESV.names <- c("Alloiococcus ESV", "Corynebacterium ESV", "E. coli ESV", 
							 "Gemellaceae ESV", "Haemophilus ESV2", "Haemophilus ESV1", 
							 "Moraxella ESV", "Moraxellaceae ESV1", "Pseudomonas ESV",
							 "Staphylococcus ESV", "Streptococcus ESV2", "Streptococcus ESV1", 
							 "Streptococcus ESV3", "Alloiococcus rare ESVs", "Corynebacterium rare ESVs", 
							 "Haemophilus rare ESVs", "Moraxella rare ESVs", "Moraxellaceae rare ESVs", 
							 "Staphylococcus rare ESVs", "Streptococcus rare ESVs", "Other ESVs")

pdf("output/2018-02-06_MPG_by_ESV.pdf", width=16, height=8)
log.esv.df.m %>% 
	filter(MPG.name %in% c("Haemophilus MPG1", "Haemophilus MPG2", "Moraxella MPG", 
												 "Streptococcus MPG1", "Mixed MPG", "Staphylococcus MPG", 
												 "Alloiococcus MPG", "Corynebacterium MPG")) %>%
	ggplot(aes(x=ESV.simp, y=Rel.Abundance, fill=ESV)) + 
	geom_boxplot() + 
	scale_fill_manual(values = ESV.colour.hex,
										labels = paste0(c(1:21), " - ", ESV.names)) + 
	facet_wrap(~MPG.name, nrow=2) + 
	labs(x = "ESVs",
			 y = "Relative abundance",
			 fill = "ESVs") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	guides(fill = guide_legend(ncol = 1)) 
dev.off()

pdf("output/2018-02-06_MPG_by_ESV_all.pdf", width=16, height=8)
log.esv.df.m %>% 
	ggplot(aes(x=ESV.simp, y=Rel.Abundance, fill=ESV)) + 
	geom_boxplot() + 
	scale_fill_manual(values = ESV.colour.hex,
										labels = paste0(c(1:21), " - ", ESV.names)) + 
	facet_wrap(~MPG.name, nrow=2) + 
	labs(x = "ESVs",
			 y = "Relative abundance",
			 fill = "ESVs") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	guides(fill = guide_legend(ncol = 1)) 
dev.off()

save.image("output/2018-02-27_core_common_microbiome.Rdata")
