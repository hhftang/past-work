# ==============================================================================
# 2018-01-29 - Microbiome data processing and analysis
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
											"gtools")

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

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("output/2018-01-29_import_biom_format.Rdata")

# ------------------------------------------------------------------------------
# Quality control
# ------------------------------------------------------------------------------

coast.biom.orig <- coast.biom

# Pick read count threshold per sample
# ---------------------------------------
# Get total number of reads per sample
read.sum.df <- data.frame(SampleID = names(sample_sums(coast.biom.orig)),
													read.sum = sample_sums(coast.biom.orig), 
													row.names = NULL)
# sample_sums(coast.biom.orig)[1]
# sum(coast.otu[,1])

# merge with metadata
read.sum.merge.df <- merge(read.sum.df, spec.subj.table.df, by="SampleID")

pdf("output/2018-01-05_coast_negcontrol_read.sums.pdf", width=6, height=3.75)
read.sum.merge.df %>%
	mutate(type = factor(type, 
											 levels=c("COAST", "failed", "negcontrol", "CAS", "misc"),
											 labels=c("COAST", "Failed", "Negative control", "CAS", "Miscellaneous"))) %>%
	ggplot(aes(x= type,y=log10(read.sum))) + geom_boxplot() + 
	geom_hline(yintercept=log10(10000), linetype="dashed") + 
	geom_hline(yintercept=log10(4000), linetype="dashed") +
	geom_hline(yintercept=log10(3000), linetype="dashed") + 
	xlab("Type of sample") + 
	ylab("Total number of reads (log 10)")
dev.off()

df <- with(read.sum.merge.df, data.frame(type, read.sum))
# coast vs. negcontrol 
a <- with(df, prop.table(table(type, read.sum < 3000),1)) # 6% vs. 76%
b <- with(df, prop.table(table(type, read.sum < 10000),1)) # 12% vs. 88%
c <- with(df, prop.table(table(type, read.sum < 4000),1)) # 7% vs. 79%
write.table(a, "output/2018-01-23_proportions_read_threshold_3000.txt", quote=F)
write.table(b, "output/2018-01-23_proportions_read_threshold_10000.txt", quote=F)
write.table(c, "output/2018-01-23_proportions_read_threshold_4000.txt", quote=F)

# Compare number of reads across different runs
pdf("output/2018-01-05_runs_read.sums.pdf", width=8, height=5)
read.sum.merge.df %>%
	# filter(type == "coast") %>%
	ggplot(aes(x=factor(RunName), y=log10(read.sum))) + geom_boxplot() + 
	geom_hline(yintercept=log10(10000), linetype="dashed") + 
	geom_hline(yintercept=log10(3000), linetype="dashed") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	xlab("Run") + 
	ylab("Total number of reads (log 10)") + 
	facet_grid(type~.)
dev.off()

# Remove samples with total reads below threshold
# ---------------------------------------
length(which(sample_sums(coast.biom) < 4000)) # 442 NP to remove
# of which, 269 are COAST, 166 are negcontrol, 5 (all) failed, 2 CAS.
coast.biom <- subset_samples(coast.biom, sample_sums(coast.biom) >= 4000)
coast.biom # 12464 taxa, 3473 samples

# Remove positive and negative controls
# ---------------------------------------
# dim(with(spec.subj.table.df[grepl("coast", spec.subj.table.df$SampleID),],
# 			 data.frame(SampleID, specimen)))
# dim(with(spec.subj.table.df[!grepl("VI", spec.subj.table.df$specimen),],
# 			 data.frame(SampleID, specimen)))

coast.biom.coast <- subset_samples(coast.biom, 
																	 sample_data(coast.biom)$type == "COAST")
coast.biom.coast # 12464 taxa, 3416 samples
which(sample_sums(coast.biom.coast) == 0)

# Sum repeated samples
# ---------------------------------------
# Create new biom file without any metadata besides LogNo
coast.biom.simp <- coast.biom
sample_data(coast.biom.simp) <- sample_data(
	data.frame(LogNo = sample_data(coast.biom)$LogNo,
						 dummy = 0,
						 row.names = row.names(sample_data(coast.biom))))
dim(otu_table(coast.biom.simp)) # 12464 x 3473

# Modify metadata of samples from same specimen
sample_data(coast.biom)[which(sample_data(coast.biom)$dna %in% c(230, 1369)),]
sample_data(coast.biom)[which(sample_data(coast.biom)$dna %in% c(19, 3116)),]
coast.biom.simp.meta.df <- sample_data(coast.biom.simp)
coast.biom.simp.meta.df$LogNo <- as.character(coast.biom.simp.meta.df$LogNo)
coast.biom.simp.meta.df$LogNo[which(
	coast.biom.simp.meta.df$LogNo == 1369)] <- 230
coast.biom.simp.meta.df$LogNo[which(
	coast.biom.simp.meta.df$LogNo == 3116)] <- 19
coast.biom.simp.meta.df$LogNo <- factor(coast.biom.simp.meta.df$LogNo)
sample_data(coast.biom.simp) <- sample_data(coast.biom.simp.meta.df)

# merge samples (sum reads)
coast.biom.simp.sum <- merge_samples(coast.biom.simp, "LogNo")
dim(otu_table(coast.biom.simp.sum)) # 3062 x 12464

# Transpose matrix back to taxa_as_rows
otu_table(coast.biom.simp.sum) <- t(as.matrix(
	otu_table(coast.biom.simp.sum)))
dim(otu_table(coast.biom.simp.sum)) # 12464 x 3062
 
# create new dataframe for new metadata
df <- spec.subj.table.df[,6:ncol(spec.subj.table.df)]
df$specimen <- NULL
dim(df)
dim(unique(df))
df2 <- unique(df)
df2 <- df2[order(df2$LogNo),]
length(row.names(sample_data(coast.biom.simp.sum))) #3316
length(df2$LogNo) #3332

# truncate metadata dataframe to contain only those samples in the new summed phyloseq object
miss.sum <- mixedsort(setdiff(df2$LogNo, row.names(sample_data(coast.biom.simp.sum))))
length(miss.sum)
with(df2[df2$LogNo %in% miss.sum, ], data.frame(LogNo, subj, type, visittype))
# Neg controls and 9 COAST LogNo are missing (7 overlap with 9 missing COAST samples)
# Also, 2 COAST LogNo = two samples 1369 (=230) and 3116 (=19)

df3 <- df2[which(!(df2$LogNo %in% miss.sum)), ]
row.names(df3) <- df3$LogNo

which(row.names(sample_data(coast.biom.simp.sum)) != row.names(sample_data(df3)))

# create new phyloseq object with proper metadata
sample_data(coast.biom.simp.sum) <- sample_data(df3)

# check metadata example
sample_data(coast.biom.simp.sum)[sample_data(coast.biom.simp.sum)$LogNo == 456,]$prob_f3
sample_data(coast.biom)[sample_data(coast.biom)$LogNo == 456,]$prob_f3

coast.biom.simp.sum # 12464 taxa, 3062 samples, 100 features 

# Repeat process for COAST-only samples
# ---------------------------------------
# Create new biom file without any metadata besides LogNo
coast.biom.coast.simp <- coast.biom.coast
sample_data(coast.biom.coast.simp) <- sample_data(
	data.frame(LogNo = sample_data(coast.biom.coast)$LogNo,
						 dummy = 0,
						 row.names = row.names(sample_data(coast.biom.coast))))
dim(otu_table(coast.biom.coast.simp)) # 12464 x 3416

# Modify metadata of samples from same specimen
sample_data(coast.biom.coast)[which(sample_data(coast.biom.coast)$dna %in% c(230, 1369)),]
sample_data(coast.biom.coast)[which(sample_data(coast.biom.coast)$dna %in% c(19, 3116)),]
coast.biom.coast.simp.meta.df <- sample_data(coast.biom.coast.simp)
coast.biom.coast.simp.meta.df$LogNo <- as.character(coast.biom.coast.simp.meta.df$LogNo)
coast.biom.coast.simp.meta.df$LogNo[which(
	coast.biom.coast.simp.meta.df$LogNo == 1369)] <- 230
coast.biom.coast.simp.meta.df$LogNo[which(
	coast.biom.coast.simp.meta.df$LogNo == 3116)] <- 19
sample_data(coast.biom.coast.simp) <- sample_data(coast.biom.coast.simp.meta.df)

# merge samples (sum reads)
coast.biom.coast.simp.sum <- merge_samples(coast.biom.coast.simp, "LogNo")
dim(otu_table(coast.biom.coast.simp.sum)) # 3007 x 12464

# Transpose matrix back to taxa_as_rows
otu_table(coast.biom.coast.simp.sum) <- t(as.matrix(
	otu_table(coast.biom.coast.simp.sum)))
dim(otu_table(coast.biom.coast.simp.sum)) # 12464 x 3007

# create new dataframe for new metadata
df <- spec.subj.table.df[,6:ncol(spec.subj.table.df)]
df$specimen <- NULL
dim(df)
dim(unique(df))
df2 <- unique(df)
df2 <- df2[order(df2$LogNo),]
length(row.names(sample_data(coast.biom.coast.simp.sum))) #3007
length(df2$LogNo) #3332

# truncate metadata dataframe to contain only those samples in the new summed phyloseq object
miss.sum <- mixedsort(setdiff(df2$LogNo, row.names(sample_data(coast.biom.coast.simp.sum))))
miss.sum
with(df2[df2$LogNo %in% miss.sum, ], data.frame(LogNo, subj, type, visittype))
# Neg controls and 9 COAST LogNo are missing (7 overlap with 9 missing COAST samples)
# Also, 2 COAST LogNo = two samples 1369 (=230) and 3116 (=19)

df3 <- df2[which(!(df2$LogNo %in% miss.sum)), ]
row.names(df3) <- df3$LogNo

# create new phyloseq object with proper metadata
sample_data(coast.biom.coast.simp.sum) <- sample_data(df3)

# check metadata example
sample_data(coast.biom.coast.simp.sum)[sample_data(coast.biom.coast.simp.sum)$LogNo == 456,]$prob_f3
sample_data(coast.biom.coast)[sample_data(coast.biom.coast)$LogNo == 456,]$prob_f3

coast.biom.coast.simp.sum # 12464 taxa, 3007 samples 
# (vs. 3147 from original dataset)

# It is likely that 8 of the 10 missing samples correspond to the 9 zero-read
# samples from before (one - individual 706 - had two zero-read samples)
# The other two samples are 1369 (= 230) and 3116 (= 19)

coast.biom.sum <- coast.biom.simp.sum
coast.biom.coast.sum <- coast.biom.coast.simp.sum
# 
# # Repeat threshold graphs
# # ---------------------------------------
# read.sum.df <- data.frame(LogNo = names(sample_sums(coast.biom.sum)),
# 													read.sum = sample_sums(coast.biom.sum), 
# 													row.names = NULL)
# 
# read.sum.merge.df <- merge(read.sum.df, spec.subj.table.df, by="LogNo")
# 
# pdf("output/2018-01-05_coast_negcontrol_read.sums.2.pdf", width=6, height=3.75)
# read.sum.merge.df %>%
# 	mutate(type = factor(type, 
# 											 levels=c("COAST", "failed", "negcontrol", "CAS", "misc"),
# 											 labels=c("COAST", "Failed", "Negative control", "CAS", "Miscellaneous"))) %>%
# 	ggplot(aes(x= type,y=log10(read.sum))) + geom_boxplot() + 
# 	geom_hline(yintercept=log10(10000), linetype="dashed") + 
# 	geom_hline(yintercept=log10(4000), linetype="dashed") +
# 	geom_hline(yintercept=log10(3000), linetype="dashed") + 
# 	xlab("Type of sample") + 
# 	ylab("Total number of reads (log 10)")
# dev.off()
# 
# df <- with(read.sum.merge.df, data.frame(type, read.sum))
# # coast vs. negcontrol 
# a <- with(df, prop.table(table(type, read.sum < 3000),1)) # 3% vs. 63%
# b <- with(df, prop.table(table(type, read.sum < 10000),1)) # 7% vs. 78%
# c <- with(df, prop.table(table(type, read.sum < 4000),1)) # 4% vs. 65%
# write.table(a, "output/2018-01-23_proportions_read_threshold_3000_2.txt", quote=F)
# write.table(b, "output/2018-01-23_proportions_read_threshold_10000_2.txt", quote=F)
# write.table(c, "output/2018-01-23_proportions_read_threshold_4000_2.txt", quote=F)

# ------------------------------------------------------------------------------
# Further quality control - read count and healthy samples too close to infectious
# ------------------------------------------------------------------------------
# From Teo et al 2017
# We therefore removed 142 NP samples with <3000 taxonomy-assigned reads, and 
# also 226 healthy samples that did not fulfil the criterion of >4 weeks after 
# an illness episode, leaving 3014 samples (1018 healthy samples and 1996 ARI 
# samples) for further analysis. 

# Reclassify samples based on score (see Shumei powerpoint)
# ---------------------------------------
with(sample_data(coast.biom.coast.sum), table(visittype, Category, useNA="ifany"))
# Sick samples: 32 with score < 5 --> treat as mild; 
# 	2 with score NA --> remove
# Well/Sick: 134 with score >= 5 --> reclassify as Sick
# 	11 with score NA --> remove
# 	1 with core 0 --> treat as well
# 	Rest (302): treat as mild

meta.df <- sample_data(coast.biom.coast.sum)
meta.df$category2 <- meta.df$visittype
meta.df$category2[which(meta.df$visittype == "Sick" & 
													meta.df$Category == "less5")] <- "Mild"
meta.df$category2[which(meta.df$visittype == "Sick" & 
													is.na(meta.df$Category))] <- "remove"
meta.df$category2[which(meta.df$visittype == "Well/Sick")] <- "Mild"
meta.df$category2[which(meta.df$visittype == "Well/Sick" & 
													is.na(meta.df$Category))] <- "remove"
meta.df$category2[which(meta.df$visittype == "Well/Sick" & 
													meta.df$Category == 0)] <- "Well"
meta.df$category2[which(meta.df$visittype == "Well/Sick" & 
													meta.df$Category == "more5")] <- "Sick"
table(meta.df$category2)
sample_data(coast.biom.coast.sum) <- meta.df

coast.biom.coast.qc <- subset_samples(coast.biom.coast.sum, 
																			 sample_data(coast.biom.coast.sum)$category2 != "remove")
coast.biom.coast.qc # 12464 taxa, 3123 samples

# Remove healthy samples within 4 weeks after illness
# ---------------------------------------
# Sort by individual, age of sample
sample.age.df <- with(sample_data(coast.biom.coast.qc), # do NOT use coast.biom.coast.qc2
											data.frame(subj, LogNo, category2, age.day))
sample.age.df <- sample.age.df[order(sample.age.df$subj, sample.age.df$age),]
row.names(sample.age.df) <- NULL

# Get list of subjects
subj <- with(sample.age.df, unique(subj))
subj <- subj[!is.na(subj)]
length(subj) # 283

sample.age.df$bad <- 0

# For each subject
for (i in 1:length(subj)) {
	# Subset dataframe
	sample.age.df.ind <- sample.age.df[which(sample.age.df$subj == subj[i]),]
	# Get indices of all controls and infectious
	c <- which(sample.age.df.ind$category2 == "Well")
	k <- which(sample.age.df.ind$category2 == "Sick" | sample.age.df.ind$category2 == "Mild")
	for (j in 1:length(c)) {
		# Get index of control
		c.focus <- c[j]
		# get all infectious category before it
		k.2 <- k[which(k <  c.focus)]
		# Get largest k if it is a finite value
		k.focus <- ifelse(is.finite(max(k.2)), max(k.2), NA)
		# Test if c.focus age is within 28 days of k.focus age, if k.focus is finite (not NA)
		if (!is.na(k.focus)){
			if (sample.age.df.ind[c.focus, "age.day"] - 
					sample.age.df.ind[k.focus, "age.day"] <= 28) {
				# get log.no of c.focus
				l.focus <- sample.age.df.ind[c.focus, "LogNo"]
				# write to bad column
				sample.age.df[which(sample.age.df$LogNo == l.focus), "bad"] <- 1
			}}
	}	
}	

# Get logs of all bad control samples
bad.logs <- sample.age.df$LogNo[which(sample.age.df$bad == 1)] 
length(bad.logs) # 72

# Get logs of all good samples
not.bad.logs <- sample.age.df$LogNo[which(sample.age.df$bad == 0)] 
length(not.bad.logs) # 3051

coast.biom.coast.qc2 <- subset_samples(coast.biom.coast.qc, 
																			sample_data(coast.biom.coast.qc)$LogNo %in% not.bad.logs)
coast.biom.coast.qc2 #12464 taxa, 3051 samples, 166 sample variables)



# Copy-number correction
# ---------------------------------------
# Read counts were corrected for OTU-specific copy number using Picrust v1.0 
# (Langille et al., 2013) using the pre-computed copy number estimates for 
# Greengenes OTUs version 13_05

# To apply to qiime2 pipeline? ignore for now

# Normalisation
# ---------------------------------------

# ** QC BEFORE NORMALISATION **

# From Teo 2017:
# "We normalized the copy-number-corrected OTU read counts using cumulative sum 
# scaling (CSS) (Paulson et al., 2013). Briefly, for each sample, the OTU counts 
# were divided by the cumulative sum of counts up to the smallest percentile for 
# which sample-specific count distributions were largely invariant (98.9th 
# percentile for our data). 

# Conventional method of "correction" - using relative abundance or total
# sum normalisation (TSS); TSS normalizes count data by dividing feature read 
# counts by the total number of reads in each sample, i.e., converts feature 
# counts to appropriately scaled ratios. TSS has been shown to incorrectly bias 
# differential abundance estimates in RNAseq data derived through 
# high-throughput technologies since a few measurements (e.g., taxa or genes) 
# are sampled preferentially as sequencing yield increases, and have an undue 
# influence on normalized counts. 

# Use metagenomeSeq package
# cumNormStatFast: Cumulative sum scaling percentile selection
# Calculates the percentile for which to sum counts up to and scale by. Faster 
# version than available in cumNormStat. Deviates from methods described in 
# Nature Methods by making use of row means for reference.

ql <- cumNormStatFast(otu_table(coast.biom.coast.qc3)@.Data)
ql #  0.8892617

# cumNorm: Cumulative sum scaling normalization
# Calculates each column's quantile and calculates the sum up to and 
# including that quantile.

# phyloseq_transform_css is a more convenient solution that wraps multiple
# metagenomeSeq packages into one function that retains phyloseq format
# perform css transformation WITHOUT log2 transform

coast.biom.coast.qc3.css <- phyloseq_transform_css(coast.biom.coast.qc3, norm=T, log=F)
otu_table(coast.biom.coast.qc3)@.Data[1:10,1:10]
otu_table(coast.biom.coast.qc3.css)@.Data[1:10,1:10]
sample_sums(coast.biom.coast.qc3)[1:10]
sample_sums(coast.biom.coast.qc3.css)[1:10]

# Get relative abundances for qc'ed + core
coast.biom.coast.qc3.rel <- microbiome::transform(coast.biom.coast.qc3, "compositional")
otu_table(coast.biom.coast.qc3.rel)@.Data[1:10,1:10]
coast.biom.coast.qc3.css.rel.core <- core(coast.biom.coast.qc3.rel, detection = 0.01/100, prevalence = 50/100)
pca <- prcomp(t(otu_table(coast.biom.coast.qc3.css.rel.core)))
meta.df <- data.frame(sample_data(coast.biom.coast.qc3.css.rel.core))
pca.df <- data.frame(pca$x)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$category2))

# Get relative abundances for qc'ed + css + core
# This gives the same as above 
coast.biom.coast.qc3.css.rel <- microbiome::transform(coast.biom.coast.qc3.css, "compositional")
otu_table(coast.biom.coast.qc3.css.rel)@.Data[1:10,1:10]
coast.biom.coast.qc3.css.rel.core <- core(coast.biom.coast.qc3.css.rel, detection = 0.01/100, prevalence = 50/100)
pca <- prcomp(t(otu_table(coast.biom.coast.qc3.css.rel.core)))
meta.df <- data.frame(sample_data(coast.biom.coast.qc3.css.rel.core))
pca.df <- data.frame(pca$x)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$category2))

# ------------------------------------------------------------------------------
# Identify core microbiome
# ------------------------------------------------------------------------------

# Method 1: recommended by microbiome package, with or without css
# ---------------------------------------
coast.biom.coast.qc3.css.rel.core # 13 taxa, 2930 samples, 735 features
# coast.biom.coast.qc3.css.rel.core.rel <- microbiome::transform(coast.biom.coast.qc3.css.rel.core, "compositional")
# coast.biom.coast.qc3.css.rel.core.rel <- subset_samples(coast.biom.coast.qc3.css.rel.core.rel, 
# 																								sample_sums(coast.biom.coast.qc3.css.rel.core.rel) != 0)
# coast.biom.coast.qc3.css.rel.core.rel #  taxa, 3120 samples, 735 features
otu1.mat <- otu_table(coast.biom.coast.qc3.css.rel.core)@.Data
tax1.df <- data.frame(tax_table(coast.biom.coast.qc3.css.rel.core)@.Data)
tax1.df$taxa <- paste0(tax1.df$Genus, " ", tax1.df$Species)
tax1.df$taxa <- gsub(pattern="[a-z]__",replacement="", tax1.df$taxa)
tax1.df$taxa <- gsub(pattern="NA",replacement="", tax1.df$taxa)
tax1.df$taxa[11] <- "Gemellaceae"
tax1.df$taxa <- trimws(tax1.df$taxa, which = c("right"))
tax1.df$taxa <- paste0(tax1.df$taxa, ", ", row.names(tax1.df))
tax1.df

rownames(otu1.mat) <- tax1.df$taxa

dist1.p <- bcdist(t(otu1.mat))

pdf("output/2018-01-25_core_hclust.pdf", width=16, height=10)
h1 <- heatmap.3(otu1.mat, 
								diss = FALSE, 
								hclust.FUN = hclust,
								dendrogram = "col", 
								cluster.by.col = TRUE, 
								Rowv =FALSE, 
								trace="none", 
								density.info="none",
								# ColIndividualColors = as.character(otu.mat$MPG.colour), 
								# color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								color.FUN = seqPalette,
								labCol = NULL, keysize = 2, mapsize = 5, 
								dist.FUN = dist, 
								hclust.col = hclust(dist1.p, method = "complete"), 
								main = ""
								#breaks = c(0,0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
)
dev.off()

# Method 2: Shumei's method, with css
# ---------------------------------------
# average abundance > 0.1% & present in >20% of samples and dominate (>50%) at least 1 sample
coast.biom.coast.qc3.css.rel <- microbiome::transform(coast.biom.coast.qc3.css, "compositional")
mean_per_OTU <- apply(otu_table(coast.biom.coast.qc3.css.rel)@.Data,1,mean)
presence_per_OTU <- apply(otu_table(coast.biom.coast.qc3.css)@.Data,1,function(x){sum(x!=0)})
max_per_OTU <- apply(otu_table(coast.biom.coast.qc3.css.rel)@.Data,1,max)

indiv.name <- row.names(otu_table(coast.biom.coast.qc3.css.rel))

# 13 otus
common.otus <- indiv.name[which(mean_per_OTU>0.001 & 
																	presence_per_OTU > 602 & 
																	max_per_OTU > 0.5)]

coast.biom.coast.qc3.css.rel.common <- prune_taxa(common.otus, coast.biom.coast.qc3.css.rel)
coast.biom.coast.qc3.css.rel.common # 13 taxa, 3120 samples, 735 features

otu2.mat <- otu_table(coast.biom.coast.qc3.css.rel.common)@.Data
tax2.df <- data.frame(tax_table(coast.biom.coast.qc3.css.rel.common)@.Data)
tax2.df$taxa <- paste0(tax2.df$Genus, tax2.df$Species)
tax2.df$taxa <- gsub(pattern="[a-z]__",replacement="", tax2.df$taxa)
tax2.df$taxa <- c("Moraxella genus, d253ca966efcd811d057dc2aa6430774",
									"Streptococcus genus, 4060107ffdf52e44f9d72343b7332609",
									"Haemophilus genus, f579ddd8d7c6a3d812ca33d72110f0d6",
									"Corynebacterium genus, cb50e79d177c169ec5b9d5e7d308874e",
									"Alloiococcus genus, dd2e15140d20c84712456e4c12ece08b",
									"Haemophilus influenzae, bc0d904a0469d05c65b270998488ee91",
									"Streptococcus genus, 357551f644064bb2b49614b1caab2fa5",
									"Staphylococcus genus, 29eb8e14ce61e9c4549de0f9cf200c19",
									"Moraxellaceae family, 6028872a168fbd6bc3381494da3b116f",
									"Streptococcus genus, a3a3b519e04bf11edbb1d82497c2b34c",
									"Gemellaceae family, d8009096d67cd7b4ce2853718bb17013",
									"Escherichia coli, d2a4add6029e32135562e2da9d4ed2e1",
									"Pseudomonas genus, 092523a4aed84953633a5f48bc2f87b9")

rownames(otu2.mat) <- tax2.df$taxa

dist2.p <- bcdist(t(otu2.mat))

pdf("output/2018-01-25_common_css_hclust.pdf", width=16, height=10)
h2 <- heatmap.3(otu2.mat, 
								diss = FALSE, 
								hclust.FUN = hclust,
								dendrogram = "col", 
								cluster.by.col = TRUE, 
								Rowv =FALSE, 
								trace="none", 
								density.info="none",
								# ColIndividualColors = as.character(otu.mat$MPG.colour), 
								# color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								color.FUN = seqPalette,
								labCol = NULL, keysize = 2, mapsize = 5, 
								dist.FUN = dist, 
								hclust.col = hclust(dist2.p, method = "complete"), 
								main = ""
								#breaks = c(0,0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
)
dev.off()

# Method 3: Shumei's method, without css
# ---------------------------------------
coast.biom.coast.qc3.rel # 23441 taxa, 3120 samples, 735 features
mean_per_OTU <- apply(otu_table(coast.biom.coast.qc3.rel)@.Data,1,mean)
presence_per_OTU <- apply(otu_table(coast.biom.coast.qc3)@.Data,1,function(x){sum(x!=0)})
max_per_OTU <- apply(otu_table(coast.biom.coast.qc3.rel)@.Data,1,max)

indiv.name <- row.names(otu_table(coast.biom.coast.qc3.rel))

# 13 otus
common.otus <- indiv.name[which(mean_per_OTU>0.001 & 
																	presence_per_OTU > 602 & 
																	max_per_OTU > 0.5)]

coast.biom.coast.qc3.rel.common <- prune_taxa(common.otus, coast.biom.coast.qc3.rel)
coast.biom.coast.qc3.rel.common # 13 taxa, 3120 samples, 735 features

otu3.mat <- otu_table(coast.biom.coast.qc3.rel.common)@.Data
tax3.df <- data.frame(tax_table(coast.biom.coast.qc3.rel.common)@.Data)
tax3.df$taxa <- paste0(tax3.df$Genus, tax3.df$Species)
tax3.df$taxa <- gsub(pattern="[a-z]__",replacement="", tax3.df$taxa)
tax3.df$taxa <- c("Moraxella genus, d253ca966efcd811d057dc2aa6430774",
									"Streptococcus genus, 4060107ffdf52e44f9d72343b7332609",
									"Haemophilus genus, f579ddd8d7c6a3d812ca33d72110f0d6",
									"Corynebacterium genus, cb50e79d177c169ec5b9d5e7d308874e",
									"Alloiococcus genus, dd2e15140d20c84712456e4c12ece08b",
									"Haemophilus influenzae, bc0d904a0469d05c65b270998488ee91",
									"Streptococcus genus, 357551f644064bb2b49614b1caab2fa5",
									"Staphylococcus genus, 29eb8e14ce61e9c4549de0f9cf200c19",
									"Moraxellaceae family, 6028872a168fbd6bc3381494da3b116f",
									"Streptococcus genus, a3a3b519e04bf11edbb1d82497c2b34c",
									"Gemellaceae family, d8009096d67cd7b4ce2853718bb17013",
									"Escherichia coli, d2a4add6029e32135562e2da9d4ed2e1",
									"Pseudomonas genus, 092523a4aed84953633a5f48bc2f87b9")

rownames(otu3.mat) <- tax3.df$taxa

dist3.p <- bcdist(t(otu3.mat))

pdf("output/2018-01-25_common_no_css_hclust.pdf", width=16, height=10)
h3 <- heatmap.3(otu3.mat, 
								diss = FALSE, 
								hclust.FUN = hclust,
								dendrogram = "col", 
								cluster.by.col = TRUE, 
								Rowv =FALSE, 
								trace="none", 
								density.info="none",
								# ColIndividualColors = as.character(otu.mat$MPG.colour), 
								# color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
								color.FUN = seqPalette,
								labCol = NULL, keysize = 2, mapsize = 5, 
								dist.FUN = dist, 
								hclust.col = hclust(dist3.p, method = "complete"), 
								main = ""
								#breaks = c(0,0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
)
dev.off()
# No difference from with css

# -----------------------------------------------------------------------------
# Taxonomy confidence
# ------------------------------------------------------------------------------

taxon.conf.df <- read.table("data/taxonomy_confidence.tsv", header=F)
colnames(taxon.conf.df) <- c("ESV", "confidence")

save.image("output/2018-01-25_import_biom_format.Rdata")