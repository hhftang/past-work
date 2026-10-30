# ==============================================================================
# 2018-01-29 - Import .biom file
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

load("data/2016-10-26_vir_respchart_var_gen.Rdata")
load("data/2017-06-19_master_table_gen.Rdata")

# Information on pipeline for derivation of CAS biom data
# ---------------------------------------

# 1) import.sh: imported fastq files (CASAVA 1.8) to QIIME2-compatible format
# "SampleData[PairedEndSequencesWithQuality]"

# 2) dada.sh: dereplicates, denoises, joins paired sequences, and removes chimeras
# trimlength (5') 10, trunclength (3') 150* (i.e. no truncation), for both forward 
# and reverse reads;
# then infers exact sequence variants (SVs), taking into account noise, quality etc.
# Produces "100% OTU" output (each unique SV is given a unique label); 
# no clustering of (e.g. 99%) similar variants by uclust or the like
# (i.e. no closed reference or de novo OTU picking)
# See https://www.nature.com/articles/ismej2017119 and 
# https://www.nature.com/articles/nmeth.3869 for explanations

# 3) filthuman.sh: filter out human reads based on human NCBI reference (GRCh38)
# Exclude sequences with >=97% identity and >=97% alignment to human reference

# 4) merge_output.sh: merge files for each runs into a single QIIME2 feature table, 
# and a single list of representative sequences

# 5) classify.sh: use a SciKitlearn-based taxonomy classifier provided by the 
# creators of QIIME2. Classifier trained on Greengenes 13.8 99% OTUs, 
# where sequences were trimmed to only include 250 bases from 16S V4 region 
# bound by 515F/806R primer pair (same as our data)
# Start a highest-resolved level (SVs) --> classify / cluster into taxons (OTUs)

# Collate metadata
# ---------------------------------------

# Viral / respchart variables by sample
# Merge viral and respchart datasets
# ----

virus_resp_all.df <- merge(virus_sample_all.df, respchart.df, by=c("SUBJ_ID", "LogNo"), 
													 all=T)
# virus_resp_all.df <- virus_resp_all.df[which(virus_resp_all.df$year < 4),]

# Merge ages (no NAs)
virus_resp_all.df$age.x[is.na(virus_resp_all.df$age.x)] <- 
	virus_resp_all.df$age.y[is.na(virus_resp_all.df$age.x)]
virus_resp_all.df$age.y[is.na(virus_resp_all.df$age.y)] <- 
	virus_resp_all.df$age.x[is.na(virus_resp_all.df$age.y)]
with(virus_resp_all.df[is.na(virus_resp_all.df$age.x),], 
		 data.frame(age.x, age.y))

colnames(virus_resp_all.df)[colnames(virus_resp_all.df)=="age.y"] <- "age"
virus_resp_all.df$age.x <- NULL

# Use respchart values instead of viral
colnames(virus_resp_all.df)[colnames(virus_resp_all.df)=="uri.y"] <- "uri"
colnames(virus_resp_all.df)[colnames(virus_resp_all.df)=="lri.y"] <- "lri"
colnames(virus_resp_all.df)[colnames(virus_resp_all.df)=="wlri.y"] <- "wlri"
colnames(virus_resp_all.df)[colnames(virus_resp_all.df)=="flri.y"] <- "flri"
colnames(virus_resp_all.df)[colnames(virus_resp_all.df)=="wflri.y"] <- "wflri"
virus_resp_all.df$uri.x <- NULL
virus_resp_all.df$lri.x <- NULL
virus_resp_all.df$wlri.x <- NULL
virus_resp_all.df$flri.x <- NULL
virus_resp_all.df$wflri.x <- NULL

# Variables by subject
# ----
dim(master.df)

# merge by subj
subj.sample.table.df <- merge(virus_resp_all.df, master.df, 
															by="SUBJ_ID", all=T)

# get select metadata from miseq
miseq.meta.df <- read.table("data/metadata.tsv",
											stringsAsFactors = F, header=T)
colnames(miseq.meta.df)[which(colnames(miseq.meta.df) == "Subjid")] <- "SUBJ_ID"

# merge by specimen ID (only use those that have miseq data)
spec.subj.table.df <- merge(subj.sample.table.df, 
														miseq.meta.df, by=c("SUBJ_ID", "LogNo"), 
														all.x=F, all.y=T)
dim(spec.subj.table.df) # 3458 samples
spec.subj.table.df2 <- merge(subj.sample.table.df, 
														miseq.meta.df, by=c("SUBJ_ID", "LogNo"), 
														all=T)
dim(spec.subj.table.df2) # 5369 samples (i.e. 1911 had no NP sampling data)

# Fix / Harmonize some variables
# Category
spec.subj.table.df$Category.x <- as.character(spec.subj.table.df$Category.x)
spec.subj.table.df$Category.y <- as.character(spec.subj.table.df$Category.y)

spec.subj.table.df$Category.x[is.na(spec.subj.table.df$Category.x)] <- 
	spec.subj.table.df$Category.y[is.na(spec.subj.table.df$Category.x)]
spec.subj.table.df$Category.y[is.na(spec.subj.table.df$Category.y)] <- 
	spec.subj.table.df$Category.x[is.na(spec.subj.table.df$Category.y)]
spec.subj.table.df[is.na(spec.subj.table.df$Category.y),"Description"]

spec.subj.table.df$Category.y[
	which(spec.subj.table.df$Category.y == "Infectious")] <- "infectious"
spec.subj.table.df$Category <- factor(spec.subj.table.df$Category.y)

spec.subj.table.df$Category.x <- NULL
spec.subj.table.df$Category.y <- NULL

# Gender / sex
# (Do not need to fill all values as long as there is a gender for each given individual)
spec.subj.table.df$gender[is.na(spec.subj.table.df$gender)] <- 
	spec.subj.table.df$sex[is.na(spec.subj.table.df$gender)]
spec.subj.table.df$sex[is.na(spec.subj.table.df$sex)] <- 
	spec.subj.table.df$gender[is.na(spec.subj.table.df$sex)]
spec.subj.table.df[is.na(spec.subj.table.df$sex),"Description"]
spec.subj.table.df[is.na(spec.subj.table.df$gender),"Description"]
spec.subj.table.df[which(spec.subj.table.df$sex != spec.subj.table.df$gender),"Description"]
spec.subj.table.df$gender <- NULL

# Age
spec.subj.table.df$age[is.na(spec.subj.table.df$age)] <- 
	spec.subj.table.df$Age[is.na(spec.subj.table.df$age)]
spec.subj.table.df$Age[is.na(spec.subj.table.df$Age)] <- 
	spec.subj.table.df$age[is.na(spec.subj.table.df$Age)]
spec.subj.table.df[is.na(spec.subj.table.df$age),"Description"]
spec.subj.table.df[is.na(spec.subj.table.df$Age),"Description"]
spec.subj.table.df[which(spec.subj.table.df$age != spec.subj.table.df$Age),"Description"]
spec.subj.table.df$Age <- NULL

# Year
with(spec.subj.table.df, data.frame(year, year_new))
spec.subj.table.df[which(spec.subj.table.df$year != spec.subj.table.df$year_new),"Description"]
spec.subj.table.df$year_new2 <- as.numeric(ceiling(spec.subj.table.df$age/365))
spec.subj.table.df[which(spec.subj.table.df$year_new != spec.subj.table.df$year_new2),"Description"]
spec.subj.table.df$year <- NULL
spec.subj.table.df$year_new <- NULL
spec.subj.table.df$year <- factor(spec.subj.table.df$year_new2)
spec.subj.table.df$year_new2 <- NULL

# Date as number of days from first day of 1996
# First, retrieve correct DOB from premaster.df
premaster.snip.df <- with(premaster.df, data.frame(SUBJ_ID, 
																									 dob))

spec.subj.table.df <- merge(spec.subj.table.df, premaster.snip.df, 
														by="SUBJ_ID", all.x=T, all.y=F)

spec.subj.table.df$dob.x <- NULL
spec.subj.table.df$dob <- spec.subj.table.df$dob.y
spec.subj.table.df$dob.y <- NULL

# Next fill missing date entries
spec.subj.table.df$Date2 <- spec.subj.table.df$dob + spec.subj.table.df$age
with(spec.subj.table.df, which(Date != Date2))

spec.subj.table.df$Date <- spec.subj.table.df$Date2
spec.subj.table.df$Date2 <- NULL

# Finally, get date.days
spec.subj.table.df$date.day <- as.integer(difftime(as.POSIXct(spec.subj.table.df$Date), 
																									 as.POSIXct(as.Date("1996-01-01")), 
																									 units="days"))
# Date as year
spec.subj.table.df$date.year <- factor(
	as.numeric(ceiling(spec.subj.table.df$date.day/365)))

# Duplicated sample IDs
length(spec.subj.table.df$SampleID) # 3458
length(unique(spec.subj.table.df$SampleID)) # 3458

# Classify CAS vs. negative control vs. positive controls
spec.subj.table.df$type <- "CAS"
spec.subj.table.df[!grepl("MG", spec.subj.table.df$Description),]$type <- "negcontrol"
spec.subj.table.df[grepl("EC", spec.subj.table.df$Description),]$type <- "negcontrol"
spec.subj.table.df[which(spec.subj.table.df$Description %in% 
												 	c(1345, 1438, 1554, 1444, 1154, 1139)),]$type <- "poscontrol"
dim(spec.subj.table.df) # 3458 samples
dim(spec.subj.table.df[which(spec.subj.table.df$type == "CAS"),]) # 3382 samples


# Rearrange columns
refcols <- c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description", "RunName", "LogNo")
spec.subj.table.df <- spec.subj.table.df[, 
																				 c(refcols, 
																				 	setdiff(colnames(spec.subj.table.df), refcols))]

# Mark healthy samples within 4 weeks after illness
# ---------------------------------------
sample.age.df <- with(spec.subj.table.df, 
											data.frame(SUBJ_ID, LogNo, Category, age))
sample.age.df <- sample.age.df[order(sample.age.df$SUBJ_ID, sample.age.df$age),]
row.names(sample.age.df) <- NULL

# Get list of subjects
subj <- with(sample.age.df, unique(SUBJ_ID))
subj <- subj[!is.na(subj)]
length(subj) # 244

sample.age.df$bad <- 0

# For each subject
for (i in 1:length(subj)) {
	# Subset dataframe
	sample.age.df.ind <- sample.age.df[which(sample.age.df$SUBJ_ID == subj[i]),]
	# Get indices of all controls and infectious
	c <- which(sample.age.df.ind$Category == "control")
	k <- which(sample.age.df.ind$Category == "infectious")
	for (j in 1:length(c)) {
		# Get index of control
		c.focus <- c[j]
		# get all infectious category before it
		k.2 <- k[which(k <  c.focus)]
		# Get largest k if it is a finite value
		k.focus <- ifelse(is.finite(max(k.2)), max(k.2), NA)
		# Test if c.focus age is within 28 days of k.focus age, if k.focus is finite (not NA)
		if (!is.na(k.focus)){
			if (sample.age.df.ind[c.focus, "age"] - 
					sample.age.df.ind[k.focus, "age"] <= 28) {
				# get log.no of c.focus
				l.focus <- sample.age.df.ind[c.focus, "LogNo"]
				# write to bad column
				sample.age.df[which(sample.age.df$LogNo == l.focus), "bad"] <- 1
			}}
	}	
}		

# Get logs of all bad control samples
bad.logs <- sample.age.df$LogNo[which(sample.age.df$bad == 1)] 
length(bad.logs) # 137

# Get logs of all good samples
not.bad.logs <- sample.age.df$LogNo[which(sample.age.df$bad == 0)] 
length(not.bad.logs) # 3321

spec.subj.table.df$healthy.not.within.4wk.inf <- 0
spec.subj.table.df$healthy.not.within.4wk.inf[
	which(spec.subj.table.df$LogNo %in% not.bad.logs)] <- 1

# 137 bad vs. 3321 good

# Load CAS biom (Amplicon data)
# ---------------------------------------
# Note that the read_phyloseq function from the microbiome package
# requires biom files to contain taxonomic data,
# which QIIME2 does not normally do when exporting feature-table.qza to feature-table.biom.

# Therefore we need to add taxonomic data manually.
# See export_biom.sh in CAS project folder for more details on conversion process

# import_biom is old way of importing (via phyloseq) e.g.:
# otu.biom <- import_biom("data/feature-table-taxa.biom")

spec.subj.table.df$Description[which(spec.subj.table.df$Description == "")] <- NA
spec.subj.table.df$RunName2 <- NULL

# Rewrite metadata file to .csv which can then be read by read_phyloseq
write.csv(spec.subj.table.df, "data/metadata_full.csv",
					quote = F, row.names = F)
write.table(spec.subj.table.df[,1:7], "data/metadata_full.tsv", sep="\t",
						quote = F, row.names = F)
# Import biom to  phyloseq object using read_phyloseq with type="biom"
# NOTE: 
# read_phyloseq(otu.file = NULL, taxonomy.file = NULL, metadata.file = NULL)
# otu.file = A biom file with '.biom' extension
# taxonomy.file	= NULL the latest version has taxonomic information within the biom
# metadata.file	= A simple metadata/mapping file with .csv extension
cas.biom <- read_phyloseq(otu.file = "data/feature-table-taxa.biom", 
														metadata.file = "data/metadata_full.csv",
														taxonomy.file = NULL, 
														type = "biom")

cas.biom # 23441 taxa, 3439 samples

# Samples missing microbiome
miss.micro <- mixedsort(setdiff(unique(spec.subj.table.df$SampleID), unique(row.names(sample_data(cas.biom)))))
miss.micro 
with(spec.subj.table.df[spec.subj.table.df$SampleID %in% miss.micro,], 
		 data.frame(SampleID, RunName, LogNo, Category, type))

# All are control samples that were given unique IDs previously
# NOTE: some control samples must have been summed together during QIIME2 analysis because
# they were not given unique IDs

save.image("output/2018-02-06_import_biom_format.Rdata")

