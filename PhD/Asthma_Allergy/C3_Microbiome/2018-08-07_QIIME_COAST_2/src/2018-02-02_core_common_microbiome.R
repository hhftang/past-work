# ==============================================================================
# 2018-01-31 - Core microbiome analysis
# ==============================================================================

rm(list=ls())
dev.off()

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

load("output/2018-01-31_process_phyloseq_1.Rdata")

# ------------------------------------------------------------------------------
# Get common ESVs
# ------------------------------------------------------------------------------

# # Method 1: recommended by microbiome package, with or without css
# # ---------------------------------------
# coast.biom.method1.css.core # 13 taxa, 2926 samples
# # coast.biom.method1.css.core.rel <- microbiome::transform(coast.biom.method1.css.core, "compositional")
# # coast.biom.method1.css.core.rel <- subset_samples(coast.biom.method1.css.core.rel, 
# # 																								sample_sums(coast.biom.method1.css.core.rel) != 0)
# # coast.biom.method1.css.core.rel 
# otu1.mat <- otu_table(coast.biom.method1.css.core)@.Data
# tax1.df <- data.frame(tax_table(coast.biom.method1.css.core)@.Data)
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
# pdf("output/2018-01-31_core_hclust.pdf", width=16, height=10)
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
coast.biom.method1.css # 12464 taxa, 2926 samples
mean_per_OTU <- apply(otu_table(coast.biom.method1.css)@.Data,1,mean)
presence_per_OTU <- apply(otu_table(coast.biom.method1.css)@.Data,1,function(x){sum(x!=0)})
max_per_OTU <- apply(otu_table(coast.biom.method1.css)@.Data,1,max)
indiv.name <- row.names(otu_table(coast.biom.method1.css))

common.otus <- indiv.name[which(mean_per_OTU>0.001 & 
																	presence_per_OTU > 0.20*ncol(otu_table(coast.biom.method1.css)) & 
																	max_per_OTU > 0.5)]
length(common.otus) # 18 otus

coast.biom.method1.css.common <- prune_taxa(common.otus, coast.biom.method1.css)
coast.biom.method1.css.common # 18 taxa, 2926 samples

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
coast.biom.method1.css.common.across.cohorts <- prune_taxa(common.otus.across.cohorts, 
																													 coast.biom.method1.css)
rowMeans(otu_table(coast.biom.method1.css.common.across.cohorts)@.Data)
df <- data.frame(rowMeans(otu_table(coast.biom.method1.css.common.across.cohorts)@.Data))
df$esv <- row.names(df)
colnames(df) <- c("mean.rel.abund.coast", "esv")
write.table(df, "output/common.otus.rel.abund.txt", row.names=F, quote=F, col.names=T)


# Extract rel. abundance for each major ESV
# ---------------------------------------
# Get moraxella spp.
coast.biom.moraxella <- subset_taxa(coast.biom.method1.css, Genus == "g__Moraxella")
coast.biom.moraxella # 193 otus
moraxella.rare <- row.names(otu_table(coast.biom.moraxella))[-which(
	row.names(otu_table(coast.biom.moraxella)) %in% 
		c("d253ca966efcd811d057dc2aa6430774"))] 
coast.biom.moraxella.rare <- prune_taxa(moraxella.rare, coast.biom.moraxella)
coast.biom.moraxella.rare # 192 otus
read.sum.moraxella.rare <- sample_sums(coast.biom.moraxella.rare)

# Moraxellaceae (excl. moraxella spp.)
coast.biom.moraxellaceae <- subset_taxa(coast.biom.method1.css, Family == "f__Moraxellaceae")
coast.biom.moraxellaceae # 590 otus
moraxellaceae.rare <- row.names(otu_table(coast.biom.moraxellaceae))[-which(
	row.names(otu_table(coast.biom.moraxellaceae)) %in% 
		c(row.names(otu_table(coast.biom.moraxella)),
			"6028872a168fbd6bc3381494da3b116f", "a5a08d96bf11fdeed28fc530d7f5d8c7"))] 
coast.biom.moraxellaceae.rare <- prune_taxa(moraxellaceae.rare, coast.biom.moraxellaceae)
coast.biom.moraxellaceae.rare # 395 otus (590-192-1-2)
read.sum.moraxellaceae.rare <- sample_sums(coast.biom.moraxellaceae.rare)

# Streptococcus
coast.biom.streptococcus <- subset_taxa(coast.biom.method1.css, Genus == "g__Streptococcus")
coast.biom.streptococcus # 323 otus
streptococcus.rare <- row.names(otu_table(coast.biom.streptococcus))[-which(
	row.names(otu_table(coast.biom.streptococcus)) %in% 
		c("357551f644064bb2b49614b1caab2fa5", "4060107ffdf52e44f9d72343b7332609",
			"b069c8413c9eb478d4faa117c4372bd3", "a3a3b519e04bf11edbb1d82497c2b34c",
			"be1bac8c38e726214223a5ac13184ab7"))] 
coast.biom.streptococcus.rare <- prune_taxa(streptococcus.rare, coast.biom.streptococcus)
coast.biom.streptococcus.rare # 318 otus
read.sum.streptococcus.rare <- sample_sums(coast.biom.streptococcus.rare)

# Alloiococcus
coast.biom.alloiococcus <- subset_taxa(coast.biom.method1.css, Genus == "g__Alloiococcus")
coast.biom.alloiococcus # 132 otus
alloiococcus.rare <- row.names(otu_table(coast.biom.alloiococcus))[-which(
	row.names(otu_table(coast.biom.alloiococcus)) %in% 
		c("dd2e15140d20c84712456e4c12ece08b"))] 
coast.biom.alloiococcus.rare <- prune_taxa(alloiococcus.rare, coast.biom.alloiococcus)
coast.biom.alloiococcus.rare # 131 otus
read.sum.alloiococcus.rare <- sample_sums(coast.biom.alloiococcus.rare)

# Corynebacterium
coast.biom.corynebacterium <- subset_taxa(coast.biom.method1.css, Genus == "g__Corynebacterium")
coast.biom.corynebacterium # 305 otus
corynebacterium.rare <- row.names(otu_table(coast.biom.corynebacterium))[-which(
	row.names(otu_table(coast.biom.corynebacterium)) %in% 
		c("cb50e79d177c169ec5b9d5e7d308874e"))] 
coast.biom.corynebacterium.rare <- prune_taxa(corynebacterium.rare, coast.biom.corynebacterium)
coast.biom.corynebacterium.rare # 304 otus
read.sum.corynebacterium.rare <- sample_sums(coast.biom.corynebacterium.rare)

# Staphylococcus
coast.biom.staphylococcus <- subset_taxa(coast.biom.method1.css, Genus == "g__Staphylococcus")
coast.biom.staphylococcus # 114 otus
staphylococcus.rare <- row.names(otu_table(coast.biom.staphylococcus))[-which(
	row.names(otu_table(coast.biom.staphylococcus)) %in% 
		c("29eb8e14ce61e9c4549de0f9cf200c19"))] 
coast.biom.staphylococcus.rare <- prune_taxa(staphylococcus.rare, coast.biom.staphylococcus)
coast.biom.staphylococcus.rare # 113 otus
read.sum.staphylococcus.rare <- sample_sums(coast.biom.staphylococcus.rare)

# Veillonella
coast.biom.veillonella <- subset_taxa(coast.biom.method1.css, Genus == "g__Veillonella")
coast.biom.veillonella # 195 otus
veillonella.rare <- row.names(otu_table(coast.biom.veillonella))[-which(
	row.names(otu_table(coast.biom.veillonella)) %in% 
		c("fb81574d887fa68041a84122487a0e47"))] 
coast.biom.veillonella.rare <- prune_taxa(veillonella.rare, coast.biom.veillonella)
coast.biom.veillonella.rare # 194 otus
read.sum.veillonella.rare <- sample_sums(coast.biom.veillonella.rare)

# Haemophilus
coast.biom.haemophilus <- subset_taxa(coast.biom.method1.css, Genus == "g__Haemophilus")
coast.biom.haemophilus # 189 otus
haemophilus.rare <- row.names(otu_table(coast.biom.haemophilus))[-which(
	row.names(otu_table(coast.biom.haemophilus)) %in% 
		c("f579ddd8d7c6a3d812ca33d72110f0d6", "bc0d904a0469d05c65b270998488ee91"))] 
coast.biom.haemophilus.rare <- prune_taxa(haemophilus.rare, coast.biom.haemophilus)
coast.biom.haemophilus.rare # 187 otus
read.sum.haemophilus.rare <- sample_sums(coast.biom.haemophilus.rare)

# Others
read.sum.others.rare <- 1 - sample_sums(coast.biom.method1.css.common) -
	read.sum.moraxella.rare - read.sum.moraxellaceae.rare - read.sum.streptococcus.rare -
	read.sum.alloiococcus.rare - read.sum.corynebacterium.rare - read.sum.staphylococcus.rare -
	read.sum.haemophilus.rare
read.sum.others.rare[which(read.sum.others.rare < 0)] <- 0

other.rare <- row.names(otu_table(coast.biom.method1.css))[-which(
	row.names(otu_table(coast.biom.method1.css)) %in% 
		c(row.names(otu_table(coast.biom.method1.css.common)),
				 moraxella.rare, moraxellaceae.rare,
				 streptococcus.rare, alloiococcus.rare,
				 corynebacterium.rare, staphylococcus.rare,
				 haemophilus.rare))]
coast.biom.other <- prune_taxa(other.rare, coast.biom.method1.css)
coast.biom.other
top.other.rare <- top_taxa(coast.biom.other, n=10)
tax_table(coast.biom.other)[
	which(row.names(tax_table(coast.biom.other)) %in% top.other.rare),]

# Create OTU/ESV table
# ---------------------------------------
otu.mat.common <- otu_table(coast.biom.method1.css.common)@.Data
otu.mat <- otu_table(coast.biom.method1.css.common)@.Data

tax.df <- data.frame(tax_table(coast.biom.method1.css.common)@.Data)
tax.df$taxa <- paste0(tax.df$Genus, " ", tax.df$Species)
tax.df$taxa <- gsub(pattern="[a-z]__",replacement="", tax.df$taxa)
tax.df$taxa <- gsub(pattern="NA",replacement="", tax.df$taxa)
tax.df$taxa <- trimws(tax.df$taxa, which = c("right"))
tax.df$taxa[13] <- "Moraxellaceae"
tax.df$taxa[14] <- "Gemellaceae"
tax.df$taxa[15] <- "Moraxellaceae"
tax.df$taxa[16] <- "Neisseriaceae"
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
which(med == max(med))+1 # 12
which(m.si == max(m.si))+1 # 12
k <- which(m.si == max(m.si))+1
si <-  cluster::silhouette(cutree(grp, k=k), dist.p)

# Cut tree to number of clusters
# ---------------------------------------
MPG <- cutree(grp, k=k)
length(MPG) # 2926, for each sample
MPG.silhouette <- si[,3]
length(MPG.silhouette) # 2926, for each sample

# Define MPG names and colours
MPG.colour.scheme.df <- data.frame(
	hex = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
					'#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),
	descript = c("light blue", "blue", "light green", "green", "pink", "red",
							 "light orange", "orange", "light purple", "purple", "yellow", "brown"))

MPG.colour <- factor(MPG, levels=c(1:12),
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
MPG.colour.scheme.df$MPG.names <- 
	c("Corynebacterium cb50", "Alloiococcus dd2e/Corynebacterium cb50", 
		"Streptococcus 3575/Others", "Moraxella d253",
		"Staphylococcus 29eb", "Moraxella, rare",
		"Haemophilus f579", "Streptococcus 4060",
		"Moraxellaceae 6028", "Haemophilus influenzae bc0d",
		"Haemophilus, rare", "Moraxellaceae a5a0")

MPG.colour.scheme.df$hex2 <- 
	c("#a6cee3", "#b2df8a",
		"#ff7f00", "#e31a1c",
		"#cab2db", "#fb9a99",
		"#ffd700", "#fdbf6f",
		"#ff00ff", "#b29600",
		"#e5c100", "#ffccff")

MPG.colour.scheme.df$descript2 <- 
	c("light blue", "light green",
		"orange", "red",
		"light purple", "pink",
		"gold", "light orange",
		"magenta", "darker gold",
		"dark gold", "light magenta")

MPG.colour <- factor(MPG, levels=c(1:12),
										 labels=MPG.colour.scheme.df$hex2)

# Finalise heatmap
# ---------------------------------------
pdf("output/2018-01-25_common_hclust.pdf", width=16, height=10)
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
top.taxa <- top_taxa(coast.biom.method1.css, n=10)
top.taxa.tax.df <- tax_table(coast.biom.method1.css)[
	which(row.names(tax_table(coast.biom.method1.css)) %in% top.taxa),]
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
																		"Alloiococcus MPG",
																		"Streptococcus MPG1",
																		"Streptococcus MPG2",
																		"Haemophilus MPG1",
																		"Corynebacterium MPG",
																		"Moraxellaceae MPG1",
																		"Haemophilus MPG2",
																		"Staphylococcus MPG",
																		"Moraxellaceae MPG2",
																		"Moraxella rare MPG",
																		"Haemophilus rare MPG"))
	
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
log.table.df <- data.frame(sample_data(coast.biom.method1.css))
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
# coast.biom.method1.css.major <- prune_taxa(major.otus, 
# 																					 coast.biom.method1.css.common)

otu.tax.col <- t(otu.mat2)
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
										"#e31a1c", "#ff00ff", "#ffccff",
										"grey", "grey", "#cab2db",
										"#ff7f00", "#fdbf6f", "grey",
										"grey", "grey", "grey",
										"grey", "grey", "#e5c100", 
										"#fb9a99", "grey", "grey", 
										"grey", "grey")

# Custom labels
levels(log.esv.df.m$ESV)
ESV.names <- c("Alloiococcus ESV", "Corynebacterium ESV", "E. coli ESV", 
							 "Gemellaceae ESV", "Haemophilus ESV2", "Haemophilus ESV1", 
							 "Moraxella ESV", "Moraxellaceae ESV1", "Moraxellaceae ESV2",
							 "Neisseria ESV", "Neisseriaceae ESV", "Staphylococcus ESV",
							 "Streptococcus ESV2", "Streptococcus ESV1", "Streptococcus ESV3", 
							 "Streptococcus ESV4", "Streptococcus ESV5", "Veillonella ESV",
							 "Alloiococcus rare ESVs", "Corynebacterium rare ESVs", "Haemophilus rare ESVs",
							 "Moraxella rare ESVs", "Moraxellaceae rare ESVs", "Staphylococcus rare ESVs",
							 "Streptococcus rare ESVs", "Other ESVs")

pdf("output/2018-02-02_MPG_by_ESV.pdf", width=16, height=8)
log.esv.df.m %>% 
	filter(MPG.name %in% c("Haemophilus MPG1", "Haemophilus MPG2", "Moraxella MPG", 
												 "Streptococcus MPG1", "Streptococcus MPG2", "Staphylococcus MPG", 
												 "Alloiococcus MPG", "Corynebacterium MPG")) %>%
	ggplot(aes(x=ESV.simp, y=Rel.Abundance, fill=ESV)) + 
	geom_boxplot() + 
	scale_fill_manual(values = ESV.colour.hex,
										labels = paste0(c(1:26), " - ", ESV.names)) + 
	facet_wrap(~MPG.name, nrow=2) + 
	labs(x = "ESVs",
			 y = "Relative abundance",
			 fill = "ESVs") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	guides(fill = guide_legend(ncol = 1)) 
dev.off()


pdf("output/2018-02-02_MPG_by_ESV_all.pdf", width=16, height=8)
log.esv.df.m %>% 
	ggplot(aes(x=ESV.simp, y=Rel.Abundance, fill=ESV)) + 
	geom_boxplot() + 
	scale_fill_manual(values = ESV.colour.hex,
										labels = paste0(c(1:26), " - ", ESV.names)) + 
	facet_wrap(~MPG.name, nrow=3) + 
	labs(x = "ESVs",
			 y = "Relative abundance",
			 fill = "ESVs") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) +
	guides(fill = guide_legend(ncol = 1)) 
dev.off()


save.image("output/2018-02-02_core_common_microbiome.Rdata")



# Cluster only by these top 6
# ---------------------------------------
# coast.biom.top.taxa <- prune_taxa(top.taxa, coast.biom.method1.css)
# coast.biom.top.taxa # 6 taxa, 2926 samples
# 
# otu.mat3 <- otu_table(coast.biom.top.taxa)@.Data
# 
# tax.df <- data.frame(tax_table(coast.biom.top.taxa)@.Data)
# tax.df$taxa <- paste0(tax.df$Genus, " ", tax.df$Species)
# tax.df$taxa <- gsub(pattern="[a-z]__",replacement="", tax.df$taxa)
# tax.df$taxa <- gsub(pattern="NA",replacement="", tax.df$taxa)
# tax.df$taxa <- trimws(tax.df$taxa, which = c("right"))
# tax.df$taxa <- paste0(tax.df$taxa, ", ", row.names(tax.df))
# rownames(otu.mat3) <- tax.df$taxa
# 
# dist.p <- bcdist(t(otu.mat3))
# 
# # Define MPGs
# # ---------------------------------------
# # MPG is defined using hclust with bray curtis distance and complete linkage
# dist.mat <- as.matrix(dist.p)
# 
# # Run hclust on distance matrix from OTUs
# grp <- hclust(dist.p, method = "complete")
# 
# # Choose number of clusters based on silhouette
# b <- list()
# # Define empty vectors for filling by number of clusters
# med <- trim.med <- m.si <- trim.m.si <- rep(NA,20)
# # Cut tree by number of clusters
# # Begin with number of clusters = 2
# b <- boxplot(cluster::silhouette(cutree(grp, k=2), dist.p)[,3], 
# 						 at = 1, xlim = c(0,21), ylim = c(-1,1), 
# 						 xlab = "Number of clusters", 
# 						 ylab = "Silhouette")
# # Iterate for number of clusters = 2 to 21
# for(i in 2:21){
# 	si <- cluster::silhouette(cutree(grp, k=i), dist.p)
# 	y <- si[,3]
# 	m.si[i-1] = mean(y)
# 	trim.m.si[i-1] = mean(y[y>0])
# 	med[i-1] = median(y)
# 	trim.med[i-1] = median(y[y>0])
# 	b[[i-1]] = boxplot(si[,3], add = TRUE, at = i-1)
# }
# # Mean silhouette
# lines(1:20, m.si, lwd = 2)
# # Mean silhouette (positive only)
# lines(1:20, trim.m.si, col = 2, lwd = 2)
# axis(1, at=1:20, labels = 2:21, cex.axis=1.2)
# 
# # Choose k (number of clusters) that maximises median/mean silhouette
# which(med == max(med))+1 # 5
# which(m.si == max(m.si))+1 # 5
# k <- which(m.si == max(m.si))+1
# si <-  cluster::silhouette(cutree(grp, k=k), dist.p)
# 
# MPG <- cutree(grp, k=k)
# length(MPG) # 2926, for each sample
# MPG.silhouette <- si[,3]
# length(MPG.silhouette) # 2926, for each sample
# 
# # Define MPG names and colours
# MPG.colour.scheme.df <- data.frame(
# 	hex = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
# 					'#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),
# 	descript = c("light blue", "blue", "light green", "green", "pink", "red",
# 							 "light orange", "orange", "light purple", "purple", "yellow", "brown"))
# 
# MPG.colour <- factor(MPG, levels=c(1:k),
# 										 labels=MPG.colour.scheme.df$hex[1:k])
# 
# h3 <- heatmap.3(otu.mat3, 
# 								diss = FALSE, 
# 								hclust.FUN = hclust,
# 								dendrogram = "col", 
# 								cluster.by.col = TRUE, 
# 								Rowv =FALSE, 
# 								trace="none", 
# 								density.info="none",
# 								ColIndividualColors = as.character(MPG.colour), 
# 								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
# 								# color.FUN = seqPalette,
# 								labCol = NULL, keysize = 2, mapsize = 5, 
# 								dist.FUN = dist, 
# 								hclust.col = hclust(dist.p, method = "complete"), 
# 								main = "",
# 								breaks = c(0,0.001,
# 													 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
# 													 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0)
# )
# 
# # Manually map colours to OTU names
# # MPG.colour.scheme.df$MPG.names <- 
# # 	c("Corynebacterium cb50", "Alloiococcus dd2e/Corynebacterium cb50", 
# # 		"Streptococcus 3575/Others", "Moraxella d253",
# # 		"Staphylococcus 29eb", "Moraxella, rare",
# # 		"Haemophilus f579", "Streptococcus 4060",
# # 		"Moraxellaceae 6028", "Haemophilus influenzae bc0d",
# # 		"Haemophilus, rare", "Moraxellaceae a5a0")
# 
# MPG.colour.scheme.df <- MPG.colour.scheme.df[1:k,]
# MPG.colour.scheme.df$MPG.names <- 
# 	c("Alloiococcus dd2e/Corynebacterium cb50", 
# 		"Streptococcus 3575", 
# 		"Streptococcus 4060",
# 		"Moraxella d253",
# 		"Haemophilus f579")
# 
# # MPG.colour.scheme.df$hex2 <- 
# # 	c("#a6cee3", "#7fc97f",
# # 		"orange3", "#e41a1c",
# # 		"#beaed4", "pink",
# # 		"gold", "orange",
# # 		"magenta", "gold3",
# # 		"khaki", "magenta3")
# 
# MPG.colour.scheme.df$hex2 <-
# 	c("#7fc97f",
# 		"orange3", 
# 		"orange",
# 		"#e41a1c",
# 		"gold")
# 
# MPG.colour <- factor(MPG, levels=c(1:k),
# 										 labels=MPG.colour.scheme.df$hex2[1:k])
# 
# pdf("output/2018-01-25_common_hclust_top6.pdf", width=16, height=10)
# h3 <- heatmap.3(otu.mat3, 
# 								diss = FALSE, 
# 								hclust.FUN = hclust,
# 								dendrogram = "col", 
# 								cluster.by.col = TRUE, 
# 								Rowv =FALSE, 
# 								trace="none", 
# 								density.info="none",
# 								ColIndividualColors = as.character(MPG.colour), 
# 								color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
# 								# color.FUN = seqPalette,
# 								labCol = NULL, keysize = 2, mapsize = 5, 
# 								dist.FUN = dist, 
# 								hclust.col = hclust(dist.p, method = "complete"), 
# 								main = "",
# 								breaks = c(0,0.001,
# 													 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
# 													 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0))
# dev.off()