# ==============================================================================
# 2018-02-06 - Microbiome data processing and analysis - Method 1
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

load("output/2018-02-06_import_biom_format.Rdata")

# ------------------------------------------------------------------------------
# Quality control
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# QC Method 1: 
# Read-sum threshold qc of samples 
# --> Sum repeated samples/extractions for each log 
# --> Well/mild/sick relabel and qc of logs

# QC Method 2: 
# Well/mild/sick relabel and qc of *samples* with problem LogNo 
# --> Read-sum threshold qc of samples 
# --> Sum repeated samples/extractions for each log

# QC Method 3: 
# Well/mild/sick relabel and qc of *samples* with problem LogNo 
# --> Sum repeated samples (should give original-8 missing) 
# --> Read-sum threshold qc of *logs*
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# QC Method 1: 
# Read-sum threshold qc of samples 
# --> Sum repeated samples/extractions for each log 
# --> Well/mild/sick relabel and qc of logs
# ------------------------------------------------------------------------------

# Make back-up
cas.biom.orig <- cas.biom

# Read sum threshold per sample
# ---------------------------------------
cas.biom # 23441 taxa, 3439 samples

# Get total number of reads per sample
read.sum.df <- data.frame(SampleID = names(sample_sums(cas.biom)),
												read.sum = sample_sums(cas.biom), 
												row.names = NULL)
# Merge with metadata
read.sum.merge.df <- merge(read.sum.df, 
													 data.frame(sample_data(cas.biom), 
													 					 SampleID = names(sample_sums(cas.biom))), 
													 by="SampleID")

pdf("output/2018-02-06_cas_negcontrol_read.sums_method1.pdf", width=6, height=3.75)
read.sum.merge.df %>%
	mutate(type = factor(type, 
											 levels=c("CAS", "negcontrol", "poscontrol"),
											 labels=c("CAS", "Negative control", "Positive control"))) %>%
	ggplot(aes(x= type,y=log10(read.sum))) + geom_boxplot() + 
		geom_hline(yintercept=log10(10000), linetype="dashed") + 
		geom_hline(yintercept=log10(4000), linetype="dashed") +
		geom_hline(yintercept=log10(3000), linetype="dashed", col="red") + 
		xlab("Type of sample") + 
		ylab("Total number of reads (log 10)")
dev.off()

# cas vs. negcontrol 
a <- with(read.sum.merge.df, prop.table(table(type, read.sum < 3000),1)) # 4% vs. 78%
b <- with(read.sum.merge.df, prop.table(table(type, read.sum < 4000),1)) # 5% vs. 84%
c <- with(read.sum.merge.df, prop.table(table(type, read.sum < 10000),1)) # 11% vs. 96%
write.table(a, "output/2018-02-06_proportions_read_threshold_3000_method1.txt", quote=F)
write.table(b, "output/2018-02-06_proportions_read_threshold_4000_method1.txt", quote=F)
write.table(c, "output/2018-02-06_proportions_read_threshold_10000_method1.txt", quote=F)

# Compare number of reads across different runs
pdf("output/2018-02-06_runs_read.sums_method1.pdf", width=8, height=5)
	read.sum.merge.df %>%
		# filter(type == "cas") %>%
		ggplot(aes(x=factor(RunName), y=log10(read.sum))) + geom_boxplot() + 
		geom_hline(yintercept=log10(3000), linetype="dashed") +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	xlab("Run") + 
	ylab("Total number of reads (log 10)") + 
	facet_grid(type~.)
dev.off()

# Remove samples with total reads below threshold
# ---------------------------------------
length(which(sample_sums(cas.biom) < 3000)) # 173 samples to remove
with(sample_data(cas.biom)[which(sample_sums(cas.biom) < 3000),],
		 table(type))
# Of the 173 samples, 133 are CAS, 40 are negcontrol.
cas.biom.tr <- subset_samples(cas.biom, sample_sums(cas.biom) >= 3000)
cas.biom.tr # 23441 taxa, 3266 samples (3439-173)

# Sum repeated samples for each log
# ---------------------------------------
# Create new biom file without any metadata besides LogNo
length(sample_data(cas.biom.tr)$LogNo)
sample_data(cas.biom.tr)$LogNo[duplicated(sample_data(cas.biom.tr)$LogNo) &
															 	sample_data(cas.biom.tr)$type == "CAS"]
# No duplicated LogNo - each sample corresponds to one log

# Limit to cas-only Logs
# ---------------------------------------
cas.biom.tr.cas <- subset_samples(cas.biom.tr, 
																	sample_data(cas.biom.tr)$type == "CAS")
cas.biom.tr.cas # 23441 taxa, 3249 logs

# Remove healthy logs within 4 weeks after illness
# ---------------------------------------
sum(sample_data(cas.biom)$healthy.not.within.4wk.inf == 0) # 137
sum(sample_data(cas.biom.tr.cas)$healthy.not.within.4wk.inf == 0) # 129

cas.biom.tr.cas.qc <- subset_samples(
	cas.biom.tr.cas, 
	sample_data(cas.biom.tr.cas)$healthy.not.within.4wk.inf == 1
	)
cas.biom.tr.cas.qc # 23441 taxa, 3120 logs (3249-129)

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

# Conventional method of "normalisation" - using relative abundance or total
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

ql <- cumNormStatFast(otu_table(cas.biom.tr.cas.qc)@.Data)
ql #  0.9223881

# cumNorm: Cumulative sum scaling normalization
# Calculates each column's quantile and calculates the sum up to and 
# including that quantile.

# phyloseq_transform_css is a more convenient solution that wraps multiple
# metagenomeSeq packages into one function that retains phyloseq format
# perform css transformation WITHOUT log2 transform

cas.biom.tr.cas.qc.css <- phyloseq_transform_css(cas.biom.tr.cas.qc, norm=T, log=F)

# otu_table(cas.biom.tr.cas.qc)@.Data[1:10,1:10]
# otu_table(cas.biom.tr.cas.qc.css)@.Data[1:10,1:10]
# sample_sums(cas.biom.tr.cas.qc)[1:10]
# sample_sums(cas.biom.tr.cas.qc.css)[1:10]

# Get relative abundances and pca for qc'ed + core
cas.biom.tr.cas.qc.rel <- 
	microbiome::transform(cas.biom.tr.cas.qc, "compositional")
cas.biom.tr.cas.qc.rel.core <- 
	core(cas.biom.tr.cas.qc.rel, detection = 0.01/100, prevalence = 50/100)

pca <- prcomp(t(otu_table(cas.biom.tr.cas.qc.rel.core)))
meta.df <- data.frame(sample_data(cas.biom.tr.cas.qc.rel.core))
pca.df <- data.frame(pca$x)

pdf("output/2018-02-06_pca_method1_core.pdf", width=8, height=5)
pca.df %>% 
	ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$Category))
dev.off()

# Get relative abundances and pca for qc'ed + css + core
# This gives the same as above 
cas.biom.tr.cas.qc.css.rel <- 
	microbiome::transform(cas.biom.tr.cas.qc.css, "compositional")
cas.biom.tr.cas.qc.css.rel.core <- 
	core(cas.biom.tr.cas.qc.css.rel, detection = 0.01/100, prevalence = 50/100)

pca <- prcomp(t(otu_table(cas.biom.tr.cas.qc.css.rel.core)))
meta.df <- data.frame(sample_data(cas.biom.tr.cas.qc.css.rel.core))
pca.df <- data.frame(pca$x)

pdf("output/2018-02-06_pca_method1_css_core.pdf", width=8, height=5)
pca.df %>% 
	ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$Category))
dev.off()

cas.biom.method1 <- cas.biom.tr.cas.qc.rel
cas.biom.method1 # 23441 taxa, 3120 logs
cas.biom.method1.css <- cas.biom.tr.cas.qc.css.rel
cas.biom.method1.css # 23441 taxa, 3120 logs
cas.biom.method1.css.core <- cas.biom.tr.cas.qc.rel.core
cas.biom.method1.css.core # 9 taxa, 3120 logs

save.image("output/2018-02-06_process_phyloseq_1.Rdata")
