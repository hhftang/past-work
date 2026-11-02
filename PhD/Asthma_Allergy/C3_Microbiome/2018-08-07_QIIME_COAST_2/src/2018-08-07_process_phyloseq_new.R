# ==============================================================================
# 2018-08-07 - Microbiome data processing and analysis - New Method
# ==============================================================================

# METHOD UPDATED TO INCORPORATE SUGGESTIONS FROM LOUISE AND KAT:
# 
# lmj [4:33 PM]
# Hi Howard, Most of the samples are duplicated are as the library was run 
# across two flow cells as the first flow cell had low yield ( 11 and 25) or on 
# request from Shu Mei as she found  found systematic Ns in certain bases across 
# all sample (16, 29 and 37). For the first issue you could pool the data from 
# the duplicates and for the second issue I believe the re-run fixed the N issue 
# and so only use the data from the second run of the library. The third reason 
# for a duplicate sample is that I chose 94 samples that seemed to amplify well 
# at the PCR stage (>1 ng/ul) but did not sequence well (<10, 000 reads). The 
# DNA from these samples was re-amplified and an new library was prepared (37). 
# I have annotated you file wiht these notes and also sent across a file that 
# includes information on all af the samples from  which a  new library was 
# made. Hope this helps. Louise
# 
# kholt [4:46 PM]
# Thanks very much @lmj. So @hhftang hopefully this clarifies things... 
# to reiterate, I agree with Louise that 
# 
# - for the ones that were just reruns of same lib due to low yield (libs 11 & 
# 25), pool both runs before subjecting to read count threshold

# - for the libraries that were resequenced due to Ns in the first run (libs 16, 
# 29, 37), use only the second (good, N-free) run 

# - for the samples that were re-amplified due to low yield 
# (repeated_samples_lmj.csv), use only the second run


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
# Quality control - new method
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
# New Method: 

# For reruns of libs 11 and 25, pool both runs then read-sum threshold - "was run twice"
# For resequenced libs 16, 29, 37, use only second (good, N-free) run - "systemic Ns"
# For re-amplified samples (repeated_samples_lmj.csv), use only second run - "re-amplified, sequenced in lib 37"

# Other samples:
# Read-sum threshold qc of samples 
# --> Sum repeated samples/extractions for each log 
# --> Well/mild/sick relabel and qc of logs
# ------------------------------------------------------------------------------

# Make back-up
coast.biom.orig <- coast.biom
coast.biom # 12464 taxa, 3915 samples
sample_data(coast.biom)$sampleid <- rownames(sample_data(coast.biom))


# Get louise annotations of repeated samples
# ---------------------------------------
repeated.samples.annot.df <- read.csv("data/repeated_samples_lmj.csv", 
																			stringsAsFactors = F)
head(repeated.samples.annot.df)

# Get SampleIDs of samples that should be removed
# ---------------------------------------
# Get reamplified and resequenced in lib 37)
with(repeated.samples.annot.df[
	which(repeated.samples.annot.df$comments == "sample had amplicon > 1ng/ul BUT < 10,000 reads so sample was re-amplified and sequenced in library 37"),],
	table(LogNo, RunName != "COAST16S_037"))

# Get second N-free run
with(repeated.samples.annot.df[
	which(repeated.samples.annot.df$comments == "same library was resequenced as ShuMei found systematic Ns in certain bases across samples"),],
	table(LogNo, RunName %in% c("COAST16S_016_repeat", "COAST16S_029_repeat", "COAST16S_037_repeat")))
with(repeated.samples.annot.df,
	table(LogNo, RunName %in% c("COAST16S_016_repeat", "COAST16S_029_repeat", "COAST16S_037_repeat")))

with(repeated.samples.annot.df, table(RunName, comments))

# Get problematic SampleIDs (not LogNo because there are duplicates)
samples.to.remove <- repeated.samples.annot.df$sampleid[
	which( 
		repeated.samples.annot.df$comments == "sample had amplicon > 1ng/ul BUT < 10,000 reads so sample was re-amplified and sequenced in library 37" 
		#		&	repeated.samples.annot.df$RunName != "COAST16S_037") |
		#	(repeated.samples.annot.df$comments == "same library was resequenced as ShuMei found systematic Ns in certain bases across samples" &
			| repeated.samples.annot.df$RunName %in% c("COAST16S_016", "COAST16S_029", "COAST16S_037"))]

length(samples.to.remove) # 330
unique(length(samples.to.remove)) # also 272
dim(repeated.samples.annot.df) # 965
dim(sample_data(coast.biom)) # 3915

# Remove these samples from coast.biom object
# ---------------------------------------
s <- sample_data(coast.biom)$sampleid[
	which(!('%in%'(sample_data(coast.biom)$sampleid, samples.to.remove)))]
coast.biom.prune <- subset_samples(coast.biom, sampleid %in% s)
coast.biom.prune # 12464 taxa, 3643 samples

# Test e.g. LogNo 1121
df <- data.frame(sample_data(coast.biom.prune))
df[which(df$LogNo == 1121),]

# Modify metadata of extractions from same sample (see Shumei's Aug31 powerpoint)
# ---------------------------------------
# Samples 230 and 1369 are extractions from same sample
# Samples 19 and 3116 are extractions from same sample
sample_data(coast.biom.prune)[which(sample_data(coast.biom.prune)$LogNo %in% c(230, 1369)),]
sample_data(coast.biom.prune)[which(sample_data(coast.biom.prune)$LogNo %in% c(19, 3116)),]
# Identical in almost all metadata except LogNo 

meta.df <- sample_data(coast.biom.prune)
meta.df$LogNo <- as.character(meta.df$LogNo)
meta.df$LogNo[which(meta.df$LogNo == 1369)] <- 230
meta.df$LogNo[which(meta.df$LogNo == 3116)] <- 19
meta.df$LogNo <- factor(meta.df$LogNo)
sample_data(coast.biom.prune) <- sample_data(meta.df)

# Get SampleIDs of samples that should be summed prior to QC
# ---------------------------------------
samples.to.sum <- repeated.samples.annot.df$sampleid[
	which(repeated.samples.annot.df$comments %in% 
					c("library was run twice across two flow cell due to low clustering leading to low read number on first run"))]
length(samples.to.sum) # 358
unique(length(samples.to.sum)) # 358

# Create subset of samples that need to be summed prior to QC
# ---------------------------------------
s <- sample_data(coast.biom.prune)$sampleid[
	which('%in%'(sample_data(coast.biom.prune)$sampleid, samples.to.sum))]
coast.biom.prune.tosum <- subset_samples(coast.biom.prune, sampleid %in% s)

s <- sample_data(coast.biom.prune)$sampleid[
	which(!('%in%'(sample_data(coast.biom.prune)$sampleid, samples.to.sum)))]
coast.biom.prune.others <- subset_samples(coast.biom.prune, sampleid %in% s)

coast.biom.prune.tosum # 358 samples
coast.biom.prune.others # 3227 samples 

# Sum repeated samples for each log before QC
# ---------------------------------------
coast.biom.prune.tosum.simp <- coast.biom.prune.tosum
sample_data(coast.biom.prune.tosum.simp) <- 
	sample_data(data.frame(LogNo = sample_data(coast.biom.prune.tosum.simp)$LogNo,
												 dummy = 0, # required - otherwise merge_samples won't work
												 row.names = row.names(sample_data(coast.biom.prune.tosum.simp))))
dim(otu_table(coast.biom.prune.tosum.simp)) # 12464 x 358
coast.biom.prune.tosum.simp.summed <- merge_samples(coast.biom.prune.tosum.simp, "LogNo")
coast.biom.prune.tosum.simp.summed # 180 samples
sample_data(coast.biom.prune.tosum.simp.summed)$LogNo <-
	rownames(sample_data(coast.biom.prune.tosum.simp.summed))

coast.biom.prune.tosum.summed <- coast.biom.prune.tosum.simp.summed
meta.df <- merge(data.frame(sample_data(coast.biom.prune.tosum.simp.summed)), 
				data.frame(sample_data(coast.biom.prune.tosum) %>% select(-sampleid) %>% distinct()), 
				by = "LogNo", all.x=T, all.y=F)
rownames(meta.df) <- meta.df$LogNo
sample_data(coast.biom.prune.tosum.summed) <- sample_data(meta.df)

# Test e.g. LogNo == 1121
sample_data(coast.biom.prune.tosum.summed)[
	which(sample_data(coast.biom.prune.tosum.summed)$LogNo == 1121),]
rownames(sample_data(coast.biom.prune.tosum.summed))

sample_data(coast.biom.prune.tosum.summed)$sampleid

sample_data(coast.biom.prune.others)[
	which(sample_data(coast.biom.prune.others)$sampleid == "REDACTED_SAMPLEID"),]
rownames(sample_data(coast.biom.prune.others))

sum(duplicated(sample_data(coast.biom.prune.others)$LogNo))

# Merge phyloseq objects
# ---------------------------------------

# merge_phyloseq function is bugged? susequent LogNo do not match  
sample_data(coast.biom.prune.tosum.summed)$dummy <- NULL
sample_data(coast.biom.prune.tosum.summed)$sampleid <- NA
sample_names(coast.biom.prune.tosum.summed) <- 
	as.character(sample_data(coast.biom.prune.tosum.summed)$Description)

df <- sample_data(coast.biom.prune.others)
col <- which(colnames(df) == "LogNo")
sample_data(coast.biom.prune.others) <- df[, c(col, (1:ncol(df))[-col])]

# Fix: convert numerical LogNo in coast.biom.prune.tosum.summed to factor type variable first
str(sample_data(coast.biom.prune.others)$LogNo)
str(sample_data(coast.biom.prune.tosum.summed)$LogNo)
sample_data(coast.biom.prune.tosum.summed)$LogNo <- 
	factor(sample_data(coast.biom.prune.tosum.summed)$LogNo)

coast.biom.prune.remerge <- merge_phyloseq(coast.biom.prune.tosum.summed, 
													coast.biom.prune.others)

coast.biom.prune.remerge # 3407 samples
# Note that rownames are somewhat "schizophrenic" - 
# LogNo for those from coast.biom.prune.tosum.summed,
# SampleIDs for those from coast.biom.prune.others

# Test e.g. LogNo 1121
sample_data(coast.biom.prune.remerge)[
	which(sample_data(coast.biom.prune.remerge)$LogNo == 1121),]
# No problem!

# Check read counts and get filter criteria
# ---------------------------------------
# Get total number of reads per sample (needs simple input)
read.sum.remerge.df <- data.frame(UniqueID = names(sample_sums(coast.biom.prune.remerge)),
													read.sum = sample_sums(coast.biom.prune.remerge), 
													row.names = NULL,
													stringsAsFactors = F)
# Merge with metadata
read.sum.remerge.merge.df <- merge(read.sum.remerge.df, 
													 data.frame(sample_data(coast.biom.prune.remerge), 
													 					 UniqueID = names(sample_sums(coast.biom.prune.remerge)),
													 					 stringsAsFactors = F), 
													 by="UniqueID")

pdf("output/2018-08-07_coast_negcontrol_read.sums_newmethod.pdf", width=6, height=3.75)
read.sum.remerge.merge.df %>%
	mutate(type = factor(type, 
											 levels=c("COAST", "failed", "negcontrol", "CAS", "misc"),
											 labels=c("COAST", "Failed", "Negative control", "CAS", "Miscellaneous"))) %>%
	ggplot(aes(x= type,y=log10(read.sum))) + geom_boxplot() + 
	geom_hline(yintercept=log10(10000), linetype="dashed") + 
	geom_hline(yintercept=log10(4000), linetype="dashed", col="red") +
	geom_hline(yintercept=log10(3000), linetype="dashed") + 
	xlab("Type of sample") + 
	ylab("Total number of reads (log 10)")
dev.off()

# coast vs. negcontrol 
a <- with(read.sum.remerge.merge.df, prop.table(table(type, read.sum < 3000),1)) # 4% vs. 76%
b <- with(read.sum.remerge.merge.df, prop.table(table(type, read.sum < 4000),1)) # 5% vs. 79%
c <- with(read.sum.remerge.merge.df, prop.table(table(type, read.sum < 10000),1)) # 9% vs. 88%

write.table(a, "output/2018-08-07_proportions_read_threshold_3000_newmethod.txt", quote=F)
write.table(b, "output/2018-08-07_proportions_read_threshold_4000_newmethod.txt", quote=F)
write.table(c, "output/2018-08-07_proportions_read_threshold_10000_newmethod.txt", quote=F)

# Compare number of reads across different runs
pdf("output/2018-08-07_runs_read.sums_newmethod.pdf", width=8, height=5)
read.sum.remerge.merge.df %>%
	# filter(type == "coast") %>%
	ggplot(aes(x=factor(RunName), y=log10(read.sum))) + geom_boxplot() + 
	geom_hline(yintercept=log10(4000), linetype="dashed") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	xlab("Run") + 
	ylab("Total number of reads (log 10)") + 
	facet_grid(type~.)
dev.off()

# Remove samples with total reads below threshold
# ---------------------------------------
length(which(sample_sums(coast.biom.prune.remerge) < 4000)) # 338 samples to remove
with(sample_data(coast.biom)[which(sample_sums(coast.biom.prune.remerge) < 4000),],
		 table(type)) #
# Of the 338 samples, 326 are COAST, 9 are negcontrol, 1 misc, 2 are CAS.
coast.biom.prune.remerge.tr <- subset_samples(coast.biom.prune.remerge, sample_sums(coast.biom.prune.remerge) >= 4000)
coast.biom.prune.remerge.tr # 12464 taxa, 3069 samples (3407 - 338)

# Sum repeated samples for each log after QC
# ---------------------------------------
# Simplify
coast.biom.prune.remerge.tr.simp <- coast.biom.prune.remerge.tr
sample_data(coast.biom.prune.remerge.tr.simp) <- 
	sample_data(data.frame(LogNo = sample_data(coast.biom.prune.remerge.tr)$LogNo,
												 dummy = 0, # required - otherwise merge_samples won't work
												 row.names = row.names(sample_data(coast.biom.prune.remerge.tr))))
dim(otu_table(coast.biom.prune.remerge.tr.simp)) # 12464 x 3069

# Merge samples - sum reads for each LogNo
coast.biom.prune.remerge.tr.sum.simp <- merge_samples(coast.biom.prune.remerge.tr.simp, "LogNo")
coast.biom.prune.remerge.tr.sum.simp # 3058 samples
sample_data(coast.biom.prune.remerge.tr.sum.simp)$LogNo <-
	rownames(sample_data(coast.biom.prune.remerge.tr.sum.simp))

# Merge with rest of metadata
coast.biom.prune.remerge.tr.sum <- coast.biom.prune.remerge.tr.sum.simp

sample_data(coast.biom.prune.remerge.tr)[which(sample_data(coast.biom.prune.remerge.tr)$LogNo == 19),]

df <- data.frame(sample_data(coast.biom.prune.remerge.tr) %>% 
								 	select(-sampleid, -BarcodeSequence,-RunName, -dna) %>% distinct())
dim(df) # 3058 x 170

df$LogNo[which(duplicated(df$LogNo))] # none
df[which(df$LogNo == 1121),]

x <- sample_data(df)
sample_names(x) <- df$LogNo
sample_data(coast.biom.prune.remerge.tr.sum) <- x

coast.biom.prune.remerge.tr.sum.simp
coast.biom.prune.remerge.tr.sum # 12464 taxa, 3058 samples


coast.biom.prune.remerge.tr.sum.coast <- subset_samples(coast.biom.prune.remerge.tr.sum, 
																					sample_data(coast.biom.prune.remerge.tr.sum)$type == "COAST")
coast.biom.prune.remerge.tr.sum.coast # 12464 taxa, 3003 logs

# Relabel logs, and remove ambiguous logs based on score (see Shumei Aug31 powerpoint)
# ---------------------------------------
with(sample_data(coast.biom.prune.remerge.tr.sum.coast), 
		 table(visittype, Category, useNA="ifany"))
# Sick samples: 
# 	32 with score < 5 --> treat as mild; 
# 	2 with score NA --> remove
# Well/Sick: 
# 	134 with score >= 5 --> reclassify as Sick
# 	10 with score NA --> remove
# 	1 with core 0 --> treat as well
# 	Rest (295): treat as mild

# Relabelling
with(sample_data(coast.biom.prune.remerge.tr.sum.coast), 
		 table(visittype, visittype2, useNA="ifany"))
# Find ambiguous logs
length(which(sample_data(coast.biom.prune.remerge.tr.sum.coast)$visittype2 == "remove")) 
# 12 = 2 + 10
# Remove ambiguous
coast.biom.prune.remerge.tr.sum.coast.qc <- 
	subset_samples(coast.biom.prune.remerge.tr.sum.coast, 
								 sample_data(coast.biom.prune.remerge.tr.sum.coast)$visittype2 != "remove")
coast.biom.prune.remerge.tr.sum.coast.qc # 12464 taxa, 2991 logs (3003 - 12)

# Remove healthy logs within 4 weeks after illness
# ---------------------------------------
coast.biom.prune.remerge.tr.sum.coast.qc2 <- subset_samples(
	coast.biom.prune.remerge.tr.sum.coast.qc, 
	sample_data(coast.biom.prune.remerge.tr.sum.coast.qc)$healthy.not.within.4wk.inf == 1
)
coast.biom.prune.remerge.tr.sum.coast.qc2 # 12464 taxa, 2922 logs (2991 - 69)

# # ------------------------------------------------------------------------------
# # Quality control - original method
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# # QC Method 1: 
# # Read-sum threshold qc of samples 
# # --> Sum repeated samples/extractions for each log 
# # --> Well/mild/sick relabel and qc of logs
# 
# # QC Method 2: 
# # Well/mild/sick relabel and qc of *samples* with problem LogNo 
# # --> Read-sum threshold qc of samples 
# # --> Sum repeated samples/extractions for each log
# 
# # QC Method 3: 
# # Well/mild/sick relabel and qc of *samples* with problem LogNo 
# # --> Sum repeated samples (should give original-8 missing) 
# # --> Read-sum threshold qc of *logs*
# # ------------------------------------------------------------------------------
# 
# # ------------------------------------------------------------------------------
# # QC Method 1: 
# # Read-sum threshold qc of samples 
# # --> Sum repeated samples/extractions for each log 
# # --> Well/mild/sick relabel and qc of logs
# # ------------------------------------------------------------------------------
# 
# # Make back-up
# coast.biom.orig <- coast.biom
# 
# # Read sum threshold per sample
# # ---------------------------------------
# coast.biom # 12464 taxa, 3915 samples
# 
# # Get total number of reads per sample
# read.sum.df <- data.frame(SampleID = names(sample_sums(coast.biom)),
# 													read.sum = sample_sums(coast.biom), 
# 													row.names = NULL)
# # Merge with metadata
# read.sum.merge.df <- merge(read.sum.df, 
# 													 data.frame(sample_data(coast.biom), 
# 													 					 SampleID = names(sample_sums(coast.biom))), 
# 													 by="SampleID")
# 
# pdf("output/2018-01-05_coast_negcontrol_read.sums_method1.pdf", width=6, height=3.75)
# read.sum.merge.df %>%
# 	mutate(type = factor(type, 
# 											 levels=c("COAST", "failed", "negcontrol", "CAS", "misc"),
# 											 labels=c("COAST", "Failed", "Negative control", "CAS", "Miscellaneous"))) %>%
# 	ggplot(aes(x= type,y=log10(read.sum))) + geom_boxplot() + 
# 	geom_hline(yintercept=log10(10000), linetype="dashed") + 
# 	geom_hline(yintercept=log10(4000), linetype="dashed", col="red") +
# 	geom_hline(yintercept=log10(3000), linetype="dashed") + 
# 	xlab("Type of sample") + 
# 	ylab("Total number of reads (log 10)")
# dev.off()
# 
# # coast vs. negcontrol 
# a <- with(read.sum.merge.df, prop.table(table(type, read.sum < 3000),1)) # 6% vs. 76%
# b <- with(read.sum.merge.df, prop.table(table(type, read.sum < 4000),1)) # 7% vs. 79%
# c <- with(read.sum.merge.df, prop.table(table(type, read.sum < 10000),1)) # 12% vs. 88%
# write.table(a, "output/2018-01-23_proportions_read_threshold_3000_method1.txt", quote=F)
# write.table(b, "output/2018-01-23_proportions_read_threshold_4000_method1.txt", quote=F)
# write.table(c, "output/2018-01-23_proportions_read_threshold_10000_method1.txt", quote=F)
# 
# # Compare number of reads across different runs
# pdf("output/2018-01-05_runs_read.sums_method1.pdf", width=8, height=5)
# read.sum.merge.df %>%
# 	# filter(type == "coast") %>%
# 	ggplot(aes(x=factor(RunName), y=log10(read.sum))) + geom_boxplot() + 
# 	geom_hline(yintercept=log10(4000), linetype="dashed") +
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
# 	xlab("Run") + 
# 	ylab("Total number of reads (log 10)") + 
# 	facet_grid(type~.)
# dev.off()
# 
# # Remove samples with total reads below threshold
# # ---------------------------------------
# length(which(sample_sums(coast.biom) < 4000)) # 442 samples to remove
# with(sample_data(coast.biom)[which(sample_sums(coast.biom) < 4000),],
# 		 table(type))
# # Of the 442 samples, 269 are COAST, 166 are negcontrol, 5 are failed, 2 are CAS.
# coast.biom.tr <- subset_samples(coast.biom, sample_sums(coast.biom) >= 4000)
# coast.biom.tr # 12464 taxa, 3473 samples (3915-442)
# 
# # Sum repeated samples for each log
# # ---------------------------------------
# # Create new biom file without any metadata besides LogNo
# coast.biom.tr.simp <- coast.biom.tr
# sample_data(coast.biom.tr.simp) <- 
# 	sample_data(data.frame(LogNo = sample_data(coast.biom.tr)$LogNo,
# 												 dummy = 0, # required - otherwise merge_samples won't work
# 												 row.names = row.names(sample_data(coast.biom.tr))))
# dim(otu_table(coast.biom.tr.simp)) # 12464 x 3473
# 
# # Modify metadata of extractions from same sample (see Shumei's Aug31 powerpoint)
# # Samples 230 and 1369 are extractions from same sample
# # Samples 19 and 3116 are extractions from same sample
# sample_data(coast.biom.tr)[which(sample_data(coast.biom.tr)$LogNo %in% c(230, 1369)),]
# sample_data(coast.biom.tr)[which(sample_data(coast.biom.tr)$LogNo %in% c(19, 3116)),]
# # Identical in almost all metadata except LogNo 
# 
# meta.df <- sample_data(coast.biom.tr.simp)
# meta.df$LogNo <- as.character(meta.df$LogNo)
# meta.df$LogNo[which(meta.df$LogNo == 1369)] <- 230
# meta.df$LogNo[which(meta.df$LogNo == 3116)] <- 19
# meta.df$LogNo <- factor(meta.df$LogNo)
# sample_data(coast.biom.tr.simp) <- sample_data(meta.df)
# 
# # Merge samples - sum reads for each LogNo
# coast.biom.tr.simp.sum <- merge_samples(coast.biom.tr.simp, "LogNo")
# dim(otu_table(coast.biom.tr.simp.sum)) # 3062 x 12464
# 
# # Transpose matrix back to taxa_as_rows
# otu_table(coast.biom.tr.simp.sum) <- t(as.matrix(
# 	otu_table(coast.biom.tr.simp.sum)))
# dim(otu_table(coast.biom.tr.simp.sum)) # 12464 x 3062
# 
# # Create new LogNo-centric dataframe for new metadata
# df <- spec.subj.table.df[,6:ncol(spec.subj.table.df)]
# # Remove all specimen/sample-specific metadata
# df$specimen <- NULL
# dim(df) # 3936
# # Number of unique rows in dataframe (some are duplicated due to multiple samples per LogNo)
# dim(unique(df)) # 3362 
# # Get LogNo-centric dataframe
# df2 <- unique(df)
# df2 <- df2[order(df2$LogNo),]
# length(row.names(sample_data(coast.biom.tr.simp.sum))) #3062
# length(df2$LogNo) # 3332
# # Off by 300
# 
# # Truncate metadata dataframe to contain only those samples in the new summed phyloseq object
# # Get those LogNo's found in df2 but not in target sample_data
# miss.sum <- mixedsort(
# 	setdiff(df2$LogNo,row.names(sample_data(coast.biom.tr.simp.sum))))
# length(miss.sum) # 300
# # Correct offset number
# 
# df3 <- with(read.sum.merge.df[read.sum.merge.df$LogNo %in% miss.sum,], 
# 						data.frame(LogNo, subj, type, visittype, read.sum))
# # Most of these are neg controls, and COAST samples with read.sum < 4000, 
# 
# df3[which(df3$read.sum >= 4000),]
# # The only samples exceeding 4000 are the ones that are replicate extractions from 
# # the same sample.
# # LogNo 1369 = 230, 3116 = 19
# # These should already be summed up
# 
# df4 <- df2[which(!(df2$LogNo %in% miss.sum)), ]
# dim(df4) # 3062
# row.names(df4) <- df4$LogNo
# 
# coast.biom.tr.sum <- coast.biom.tr.simp.sum
# which(row.names(sample_data(coast.biom.tr.sum)) != row.names(sample_data(df4)))
# # Rownames match
# 
# # Create new phyloseq object with proper metadata
# sample_data(coast.biom.tr.sum) <- sample_data(df4)
# 
# # check metadata for selected LogNo
# sample_data(coast.biom.tr.sum)[sample_data(coast.biom.tr.sum)$LogNo == 456,]$prob_f3
# sample_data(coast.biom.tr)[sample_data(coast.biom.tr)$LogNo == 456,]$prob_f3
# 
# sample_sums(coast.biom.tr)[names(sample_sums(coast.biom.tr)) == 
# 													 	row.names(sample_data(coast.biom.tr))[
# 													 		which(sample_data(coast.biom.tr)$LogNo == 19)]] # 52319
# sample_sums(coast.biom.tr)[names(sample_sums(coast.biom.tr)) == 
# 													 	row.names(sample_data(coast.biom.tr))[
# 													 		which(sample_data(coast.biom.tr)$LogNo == 3116)]] # 70564
# sample_sums(coast.biom.tr.sum)[names(sample_sums(coast.biom.tr.sum)) == 19] # 122883
# 
# sample_sums(coast.biom.tr)[names(sample_sums(coast.biom.tr)) == 
# 													 	row.names(sample_data(coast.biom.tr))[
# 													 		which(sample_data(coast.biom.tr)$LogNo == 230)]] # 138430
# sample_sums(coast.biom.tr)[names(sample_sums(coast.biom.tr)) == 
# 													 	row.names(sample_data(coast.biom.tr))[
# 													 		which(sample_data(coast.biom.tr)$LogNo == 1369)]] # 129598
# sample_sums(coast.biom.tr.sum)[names(sample_sums(coast.biom.tr.sum)) == 230] # 268028
# # Correct summations
# 
# coast.biom.tr.sum # 12464 taxa, 3062 logs
# 
# # Limit to COAST-only Logs
# # ---------------------------------------
# coast.biom.tr.sum.coast <- subset_samples(coast.biom.tr.sum, 
# 																					sample_data(coast.biom.tr.sum)$type == "COAST")
# coast.biom.tr.sum.coast # 12464 taxa, 3007 logs
# 
# # Relabel logs, and remove ambiguous logs based on score (see Shumei Aug31 powerpoint)
# # ---------------------------------------
# with(sample_data(coast.biom.tr.sum.coast), 
# 		 table(visittype, Category, useNA="ifany"))
# # Sick samples: 
# # 	32 with score < 5 --> treat as mild; 
# # 	2 with score NA --> remove
# # Well/Sick: 
# # 	134 with score >= 5 --> reclassify as Sick
# # 	10 with score NA --> remove
# # 	1 with core 0 --> treat as well
# # 	Rest (295): treat as mild
# 
# # Relabelling
# with(sample_data(coast.biom.tr.sum.coast), 
# 		 table(visittype, visittype2, useNA="ifany"))
# # Find ambiguous logs
# length(which(sample_data(coast.biom.tr.sum.coast)$visittype2 == "remove")) 
# # 12 = 2 + 10
# # Remove ambiguous
# coast.biom.tr.sum.coast.qc <- 
# 	subset_samples(coast.biom.tr.sum.coast, 
# 								 sample_data(coast.biom.tr.sum.coast)$visittype2 != "remove")
# coast.biom.tr.sum.coast.qc # 12464 taxa, 2995 logs (3007 - 12)
# 
# # Remove healthy logs within 4 weeks after illness
# # ---------------------------------------
# coast.biom.tr.sum.coast.qc2 <- subset_samples(
# 	coast.biom.tr.sum.coast.qc, 
# 	sample_data(coast.biom.tr.sum.coast.qc)$healthy.not.within.4wk.inf == 1
# )
# coast.biom.tr.sum.coast.qc2 # 12464 taxa, 2926 logs (2995 - 69)

coast.biom.tr.sum.coast.qc2 <- coast.biom.prune.remerge.tr.sum.coast.qc2

# ------------------------------------------------------------------------------
# Normalise and finalise 
# ------------------------------------------------------------------------------

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

ql <- cumNormStatFast(otu_table(coast.biom.tr.sum.coast.qc2)@.Data)
ql #  0.8892617

# cumNorm: Cumulative sum scaling normalization
# Calculates each column's quantile and calculates the sum up to and 
# including that quantile.

# phyloseq_transform_css is a more convenient solution that wraps multiple
# metagenomeSeq packages into one function that retains phyloseq format
# perform css transformation WITHOUT log2 transform

coast.biom.tr.sum.coast.qc2.css <- phyloseq_transform_css(coast.biom.tr.sum.coast.qc2, norm=T, log=F)

# otu_table(coast.biom.tr.sum.coast.qc2)@.Data[1:10,1:10]
# otu_table(coast.biom.tr.sum.coast.qc2.css)@.Data[1:10,1:10]
# sample_sums(coast.biom.tr.sum.coast.qc2)[1:10]
# sample_sums(coast.biom.tr.sum.coast.qc2.css)[1:10]

# Get relative abundances and pca for qc'ed + core
coast.biom.tr.sum.coast.qc2.rel <- 
	microbiome::transform(coast.biom.tr.sum.coast.qc2, "compositional")
coast.biom.tr.sum.coast.qc2.rel.core <- 
	core(coast.biom.tr.sum.coast.qc2.rel, detection = 0.01/100, prevalence = 50/100)

pca <- prcomp(t(otu_table(coast.biom.tr.sum.coast.qc2.rel.core)))
meta.df <- data.frame(sample_data(coast.biom.tr.sum.coast.qc2.rel.core))
pca.df <- data.frame(pca$x)

pdf("output/2018-01-31_pca_method1_core.pdf", width=8, height=5)
pca.df %>% 
	ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$visittype2))
dev.off()

# Get relative abundances and pca for qc'ed + css + core
# This gives the same as above 
coast.biom.tr.sum.coast.qc2.css.rel <- 
	microbiome::transform(coast.biom.tr.sum.coast.qc2.css, "compositional")
coast.biom.tr.sum.coast.qc2.css.rel.core <- 
	core(coast.biom.tr.sum.coast.qc2.css.rel, detection = 0.01/100, prevalence = 50/100)

pca <- prcomp(t(otu_table(coast.biom.tr.sum.coast.qc2.css.rel.core)))
meta.df <- data.frame(sample_data(coast.biom.tr.sum.coast.qc2.css.rel.core))
pca.df <- data.frame(pca$x)

pdf("output/2018-01-31_pca_method1_css_core.pdf", width=8, height=5)
pca.df %>% 
	ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$visittype2))
dev.off()

coast.biom.method1 <- microbiome::transform(coast.biom.tr.sum.coast.qc2, "compositional")
coast.biom.method1 # 12464 taxa, 2926 logs
coast.biom.method1.css <- coast.biom.tr.sum.coast.qc2.css.rel
coast.biom.method1.css # 12464 taxa, 2926 logs
coast.biom.method1.css.core <- coast.biom.tr.sum.coast.qc2.css.rel.core
coast.biom.method1.css.core # 13 taxa, 2926 logs

save.image("output/2018-01-31_process_phyloseq_1.Rdata")
