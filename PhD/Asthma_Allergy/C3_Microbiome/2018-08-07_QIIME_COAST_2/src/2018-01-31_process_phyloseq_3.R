# ==============================================================================
# 2018-01-29 - Microbiome data processing and analysis - Method 3
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
# QC Method 3: 
# Well/mild/sick relabel and qc of *samples* with problem LogNo 
# --> Sum repeated samples (should give original-8 missing) 
# --> Read-sum threshold qc of *logs*
# ------------------------------------------------------------------------------

# Make back-up
coast.biom.orig <- coast.biom

# Relabel samples, and remove ambiguous samples based on score (see Shumei Aug31 powerpoint)
# ---------------------------------------
with(sample_data(coast.biom), 
		 table(visittype, Category, useNA="ifany"))
# Sick samples: 
# 	33 with score < 5 --> treat as mild; 
# 	4 with score NA --> remove
# Well/Sick: 
# 	159 with score >= 5 --> reclassify as Sick
# 	11 with score NA --> remove
# 	1 with core 0 --> treat as well
# 	Rest (360): treat as mild

# Relabelling
with(sample_data(coast.biom), 
		 table(visittype, visittype2, useNA="ifany"))
# Find ambiguous samples
length(which(sample_data(coast.biom)$visittype2 == "remove")) 
# 15 as per 2018-01-29_import_biom_format
# Remove ambiguous (include NAs which are esentially all non COAST samples
# dim(sample_data(coast.biom)[which(is.na(sample_data(coast.biom)$visittype2)),])
#  dim(sample_data(coast.biom)[which(sample_data(coast.biom)$type != "COAST"),])
coast.biom.qc <- 
	subset_samples(coast.biom, 
								 (is.na(sample_data(coast.biom)$visittype2) |
								  	sample_data(coast.biom)$visittype2 != "remove"))
coast.biom.qc # 12464 taxa, 3900 samples (3915-15)

# Remove healthy samples within 4 weeks after illness
# ---------------------------------------
coast.biom.qc2 <- subset_samples(
	coast.biom.qc, 
	sample_data(coast.biom.qc)$healthy.not.within.4wk.inf == 1
)
coast.biom.qc2 # 12464 taxa, 3815 samples (3900-85)

# Sum repeated samples for each log
# ---------------------------------------
# Create new biom file without any metadata besides LogNo
coast.biom.qc2.simp <- coast.biom.qc2
sample_data(coast.biom.qc2.simp) <- 
	sample_data(data.frame(LogNo = sample_data(coast.biom.qc2)$LogNo,
												 dummy = 0, # required - otherwise merge_samples won't work
												 row.names = row.names(sample_data(coast.biom.qc2))))
dim(otu_table(coast.biom.qc2.simp)) # 12464 x 3815

# Modify metadata of extractions from same sample (see Shumei's Aug31 powerpoint)
# Samples 230 and 1369 are extractions from same sample
# Samples 19 and 3116 are extractions from same sample
sample_data(coast.biom.qc2)[which(sample_data(coast.biom.qc2)$LogNo %in% c(230, 1369)),]
sample_data(coast.biom.qc2)[which(sample_data(coast.biom.qc2)$LogNo %in% c(19, 3116)),]
# Identical in almost all metadata except LogNo 

meta.df <- sample_data(coast.biom.qc2.simp)
meta.df$LogNo <- as.character(meta.df$LogNo)
meta.df$LogNo[which(meta.df$LogNo == 1369)] <- 230
meta.df$LogNo[which(meta.df$LogNo == 3116)] <- 19
meta.df$LogNo <- factor(meta.df$LogNo)
sample_data(coast.biom.qc2.simp) <- sample_data(meta.df)

# Merge samples - sum reads for each LogNo
coast.biom.qc2.simp.sum <- 
	merge_samples(coast.biom.qc2.simp, "LogNo")
dim(otu_table(coast.biom.qc2.simp.sum)) # 3258 x 12464

# Transpose matrix back to taxa_as_rows
otu_table(coast.biom.qc2.simp.sum) <- t(as.matrix(
	otu_table(coast.biom.qc2.simp.sum)))
dim(otu_table(coast.biom.qc2.simp.sum)) # 12464 x 3258

# Create new LogNo-centric dataframe for new metadata
df <- spec.subj.table.df[,6:ncol(spec.subj.table.df)]
# Remove all specimen/sample-specific metadata
df$specimen <- NULL
dim(df) # 3936
# Number of unique rows in dataframe (some are duplicated due to multiple samples per LogNo)
dim(unique(df)) # 3362 
# Get LogNo-centric dataframe
df2 <- unique(df)
df2 <- df2[order(df2$LogNo),]
length(row.names(sample_data(coast.biom.qc2.simp.sum))) #3258
length(df2$LogNo) # 3362
# Off by 104

# Truncate metadata dataframe to contain only those samples in the new summed phyloseq object
# Get those LogNo's found in df2 but not in target sample_data
miss.sum <- mixedsort(
	setdiff(df2$LogNo,row.names(sample_data(coast.biom.qc2.simp.sum))))
length(miss.sum) # 104
# Correct offset number

# Get total number of reads per sample
read.sum.df <- data.frame(SampleID = names(sample_sums(coast.biom)),
													read.sum = sample_sums(coast.biom), 
													row.names = NULL)
# Merge with metadata
read.sum.merge.df <- merge(read.sum.df, 
													 data.frame(sample_data(coast.biom), 
													 					 SampleID = names(sample_sums(coast.biom))), 
													 by="SampleID")

df3 <- with(read.sum.merge.df[read.sum.merge.df$LogNo %in% miss.sum,], 
						data.frame(LogNo, subj, type, visittype, read.sum, visittype2, healthy.not.within.4wk.inf))
# Most of these are visittype2 == "remove" or healthy.not.within.4wk.inf == "0"
with(df3, table(read.sum<4000, visittype2, healthy.not.within.4wk.inf)) 

df3[which(df3$read.sum >= 4000 &
						df3$visittype2 != "remove" & 
						df3$healthy.not.within.4wk.inf != 0),]
# The only samples that do not abide by the obvious filter contidions are the 
# ones that are replicate extractions from the same sample.
# LogNo 1369 = 230, 3116 = 19
# These should already be summed up

df4 <- df2[which(!(df2$LogNo %in% miss.sum)), ]
dim(df4) # 3258
row.names(df4) <- df4$LogNo

coast.biom.qc2.sum <- coast.biom.qc2.simp.sum
which(row.names(sample_data(coast.biom.qc2.sum)) != row.names(sample_data(df4)))
# Rownames match

# Create new phyloseq object with proper metadata
sample_data(coast.biom.qc2.sum) <- sample_data(df4)

# check metadata for selected LogNo
sample_data(coast.biom.qc2.sum)[sample_data(coast.biom.qc2.sum)$LogNo == 456,]$prob_f3
sample_data(coast.biom.qc2)[sample_data(coast.biom.qc2)$LogNo == 456,]$prob_f3

sample_sums(coast.biom.qc2)[names(sample_sums(coast.biom.qc2)) == 
															 	row.names(sample_data(coast.biom.qc2))[
															 		which(sample_data(coast.biom.qc2)$LogNo == 19)]] # 52319
sample_sums(coast.biom.qc2)[names(sample_sums(coast.biom.qc2)) == 
															 	row.names(sample_data(coast.biom.qc2))[
															 		which(sample_data(coast.biom.qc2)$LogNo == 3116)]] # 70564
sample_sums(coast.biom.qc2.sum)[names(sample_sums(coast.biom.qc2.sum)) == 19] # 122883

sample_sums(coast.biom.qc2)[names(sample_sums(coast.biom.qc2)) == 
															 	row.names(sample_data(coast.biom.qc2))[
															 		which(sample_data(coast.biom.qc2)$LogNo == 230)]] # 138430
sample_sums(coast.biom.qc2)[names(sample_sums(coast.biom.qc2)) == 
															 	row.names(sample_data(coast.biom.qc2))[
															 		which(sample_data(coast.biom.qc2)$LogNo == 1369)]] # 129598
sample_sums(coast.biom.qc2.sum)[names(sample_sums(coast.biom.qc2.sum)) == 230] # 268028
# Correct summations

coast.biom.qc2.sum # 12464 taxa, 3258 logs

# Read sum threshold per **LOG**
# ---------------------------------------
# Get total number of reads per log
read.sum.df <- data.frame(SampleID = names(sample_sums(coast.biom.qc2.sum)),
													read.sum = sample_sums(coast.biom.qc2.sum), 
													row.names = NULL)
# Merge with metadata
read.sum.merge.df <- merge(read.sum.df, 
													 data.frame(sample_data(coast.biom.qc2.sum), 
													 					 SampleID = names(sample_sums(coast.biom.qc2.sum))), 
													 by="SampleID")

pdf("output/2018-01-05_coast_negcontrol_read.sums_method3.pdf", width=6, height=3.75)
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

# coast vs. negcontrol 
a <- with(read.sum.merge.df, prop.table(table(type, read.sum < 3000),1)) # 3% vs. 73%
b <- with(read.sum.merge.df, prop.table(table(type, read.sum < 4000),1)) # 4% vs. 76%
c <- with(read.sum.merge.df, prop.table(table(type, read.sum < 10000),1)) # 7% vs. 87%
write.table(a, "output/2018-01-23_proportions_read_threshold_3000_method3.txt", quote=F)
write.table(b, "output/2018-01-23_proportions_read_threshold_4000_method3.txt", quote=F)
write.table(c, "output/2018-01-23_proportions_read_threshold_10000_method3.txt", quote=F)

# Compare number of reads across different runs
# No longer works / makes sense as we are looking per log, not per sample
# (Logs can be split across multiple runs)
# pdf("output/2018-01-05_runs_read.sums_method3.pdf", width=8, height=5)
# read.sum.merge.df %>%
# 	# filter(type == "coast") %>%
# 	ggplot(aes(x=factor(RunName), y=log10(read.sum))) + geom_boxplot() + 
# 	geom_hline(yintercept=log10(4000), linetype="dashed") + 
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
# 	xlab("Run") + 
# 	ylab("Total number of reads (log 10)") + 
# 	facet_grid(type~.)
# dev.off()

# Remove samples with total reads below threshold
# ---------------------------------------
length(which(sample_sums(coast.biom.qc2.sum) < 4000)) # 270 to remove
with(sample_data(coast.biom.qc2.sum)[which(sample_sums(coast.biom.qc2.sum) < 4000),],
		 table(type))
# Of the 270 samples, 121 are COAST, 142 are negcontrol, 5 are failed, 2 are CAS.
coast.biom.qc2.sum.tr <- subset_samples(coast.biom.qc2.sum, sample_sums(coast.biom.qc2.sum) >= 4000)
coast.biom.qc2.sum.tr # 12464 taxa, 2988 samples (3258-270)

# Limit to COAST-only Logs
# ---------------------------------------
coast.biom.qc2.sum.tr.coast <- subset_samples(coast.biom.qc2.sum.tr, 
																							sample_data(coast.biom.qc2.sum.tr)$type == "COAST")
coast.biom.qc2.sum.tr.coast # 12464 taxa, 2930 samples

# Copy-number correction
# ---------------------------------------
# To apply to qiime2 pipeline? ignore for now

# Normalisation
# ---------------------------------------
ql <- cumNormStatFast(otu_table(coast.biom.qc2.sum.tr.coast)@.Data)
ql #  0.8892617

coast.biom.qc2.sum.tr.coast.css <- 
	phyloseq_transform_css(coast.biom.qc2.sum.tr.coast, norm=T, log=F)

# otu_table(coast.biom.qc2.tr.sum.coast)@.Data[1:10,1:10]
# otu_table(coast.biom.qc2.tr.sum.coast.css)@.Data[1:10,1:10]
# sample_sums(coast.biom.qc2.tr.sum.coast)[1:10]
# sample_sums(coast.biom.qc2.tr.sum.coast.css)[1:10]

# Get relative abundances and pca for qc'ed + core
coast.biom.qc2.sum.tr.coast.rel <- 
	microbiome::transform(coast.biom.qc2.sum.tr.coast, "compositional")
coast.biom.qc2.sum.tr.coast.rel.core <- 
	core(coast.biom.qc2.sum.tr.coast.rel, detection = 0.01/100, prevalence = 50/100)

pca <- prcomp(t(otu_table(coast.biom.qc2.sum.tr.coast.rel.core)))
meta.df <- data.frame(sample_data(coast.biom.qc2.sum.tr.coast.rel.core))
pca.df <- data.frame(pca$x)

pdf("output/2018-01-31_pca_method3_core.pdf", width=8, height=5)
pca.df %>% 
	ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$visittype2))
dev.off()

# Get relative abundances and pca for qc'ed + css + core
# This gives the same as above 
coast.biom.qc2.sum.tr.coast.css.rel <- 
	microbiome::transform(coast.biom.qc2.sum.tr.coast.css, "compositional")
coast.biom.qc2.sum.tr.coast.css.rel.core <- 
	core(coast.biom.qc2.sum.tr.coast.css.rel, detection = 0.01/100, prevalence = 50/100)

pca <- prcomp(t(otu_table(coast.biom.qc2.sum.tr.coast.css.rel.core)))
meta.df <- data.frame(sample_data(coast.biom.qc2.sum.tr.coast.css.rel.core))
pca.df <- data.frame(pca$x)

pdf("output/2018-01-31_pca_method3_css_core.pdf", width=8, height=5)
pca.df %>% 
	ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$visittype2))
dev.off()

coast.biom.method3 <- microbiome::transform(coast.biom.qc2.sum.tr.coast, "compositional")
coast.biom.method3 # 12464 taxa, 2930 logs
coast.biom.method3.css <- coast.biom.qc2.sum.tr.coast.css.rel
coast.biom.method3.css # 12464 taxa, 2930 logs
coast.biom.method3.css.core <- coast.biom.qc2.sum.tr.coast.css.rel.core
coast.biom.method3.css.core # 13 taxa, 2930 logs

# Gains an extra 4 logs

save.image("output/2018-01-31_process_phyloseq_3.Rdata")
