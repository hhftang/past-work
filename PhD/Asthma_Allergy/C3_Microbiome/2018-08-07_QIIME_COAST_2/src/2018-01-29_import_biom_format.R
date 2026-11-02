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

sample.table.df <- read.csv("data/melbourne_by_sample.csv",
														stringsAsFactors = F, header=T)
subject.table.df <- read.csv("data/melbourne_by_subject.csv",
														 stringsAsFactors = F, header=T)

sample.vir.table.df <- read.csv("data/melb_vir_red.csv",
																 stringsAsFactors = F, header=T)
sample.vir.table.df$RV.type[which(sample.vir.table.df$RV.type == "")] <- NA
sample.vir.table.df$RV.type <- 
	gsub(", ", "/", sample.vir.table.df$RV.type)

subject.postBD.table.df <- read.csv("data/melbourne_spirometry_postBD.csv",
																 stringsAsFactors = F, header=T)
colnames(subject.postBD.table.df)[
	2:length(colnames(subject.postBD.table.df))] <-
	sub("^", "postBD_", colnames(subject.postBD.table.df)[
		2:length(colnames(subject.postBD.table.df))])
subject.postBD.table.df.w <- reshape(subject.postBD.table.df,
				idvar = "subj",
				timevar = "postBD_whatyear",
				direction = "wide")

subject.preBD.table.df <- read.csv("data/melbourne_spirometry_preBD.csv",
														 stringsAsFactors = F, header=T)
colnames(subject.preBD.table.df)[
	2:length(colnames(subject.preBD.table.df))] <-
	sub("^", "preBD_", colnames(subject.preBD.table.df)[
		2:length(colnames(subject.preBD.table.df))])
subject.preBD.table.df.w <- reshape(subject.preBD.table.df,
																		 idvar = "subj",
																		 timevar = "preBD_whatyear",
																		 direction = "wide")

# merge by sample
sample.table.df <- merge(sample.table.df, sample.vir.table.df, 
												 by=c("dna","subj","specimen"), all=T)
# 3147 x 27

# merge by subj
subj.sample.table.df <- merge(sample.table.df, subject.table.df, 
															by="subj", all=T) # 3147 x 97
subj.sample.table.df <- merge(subj.sample.table.df, subject.postBD.table.df.w,
															by="subj", all=T) # 3147 x 129
subj.sample.table.df <- merge(subj.sample.table.df, subject.preBD.table.df.w,
															by="subj", all=T) # 3147 x 161

subj.sample.table.df$LogNo <- subj.sample.table.df$dna
# 3147 samples, 161 features

# get select metadata from miseq
miseq.meta.df <- read.table("data/metadata_simp.tsv",
														sep = "\t", 
														stringsAsFactors = F, header=T)
# 3936 samples, 7 features
miseq.meta.df$specimen <- miseq.meta.df$Description
miseq.meta.df$specimen[which(miseq.meta.df$specimen == "")] <- NA
miseq.meta.df$LogNo <- gsub("redo", "", miseq.meta.df$LogNo)

# merge by specimen ID
spec.subj.table.df <- merge(subj.sample.table.df, 
														miseq.meta.df, by=c("specimen", "LogNo"),
														all.x=F, all.y=T)

# Rearrange columns
refcols <- c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description", "RunName", "LogNo")
spec.subj.table.df <- spec.subj.table.df[, 
																				 c(refcols, 
																				 	setdiff(colnames(spec.subj.table.df), refcols))]
dim(spec.subj.table.df) # 3936 samples, 172 features

# Fix / Harmonize some variables
# DOB
which(spec.subj.table.df$DOB.x != spec.subj.table.df$DOB.y)
spec.subj.table.df$DOB <- spec.subj.table.df$DOB.x
spec.subj.table.df$DOB.x <- NULL
spec.subj.table.df$DOB.y <- NULL

# Sick vs. well/sick vs. well samples
spec.subj.table.df$score <- 0
spec.subj.table.df$score <- (spec.subj.table.df$fever == 1) + 
	ifelse(spec.subj.table.df$severec == 3, 3,
				 ifelse(spec.subj.table.df$cough == 2, 2,
				 			 ifelse(spec.subj.table.df$mildc == 1, 1, 0))) +
	ifelse(spec.subj.table.df$rhin5 == 2, 2,
				 ifelse(spec.subj.table.df$rhin0 == 1, 1, 0)) +	
	(spec.subj.table.df$hoarse == 1) +
	3*(spec.subj.table.df$apnea == 3) +
	5*(spec.subj.table.df$wheezing == 5) +
	5*(spec.subj.table.df$cyanosis == 5) +
	5*(spec.subj.table.df$retract == 5) +
	5*(spec.subj.table.df$dyspnea == 5) +
	(spec.subj.table.df$lasting == 1)

spec.subj.table.df$Category <- 
	ifelse(spec.subj.table.df$score >= 5, "more5",
				 ifelse(spec.subj.table.df$score < 5, 
				 			 ifelse(spec.subj.table.df$score == 0, "0","less5"), NA))
spec.subj.table.df[which(spec.subj.table.df$visittype == "Sick" &
												 	spec.subj.table.df$Category == "0" | "less5"),]
with(spec.subj.table.df, table(visittype, Category, useNA="ifany"))

spec.subj.table.df[which(spec.subj.table.df$visittype == "Sick" &
												 	spec.subj.table.df$Category == "less5" & 
												 	spec.subj.table.df$subj == 222),]
# Duplicated sample IDs
length(spec.subj.table.df$SampleID) # 3936
length(unique(spec.subj.table.df$SampleID)) # 3936

# Classify CAS vs. negative control vs. positive controls
spec.subj.table.df$type <- "COAST"
spec.subj.table.df[grepl("water", spec.subj.table.df$SampleID),]$type <- "negcontrol"
spec.subj.table.df[grepl("neg", spec.subj.table.df$SampleID),]$type <- "negcontrol"
spec.subj.table.df[grepl("NEG", spec.subj.table.df$SampleID),]$type <- "negcontrol"
spec.subj.table.df[grepl("blank", spec.subj.table.df$SampleID),]$type <- "negcontrol"
spec.subj.table.df[grepl("Blank", spec.subj.table.df$SampleID),]$type <- "negcontrol"
spec.subj.table.df[grepl("balnk", spec.subj.table.df$SampleID),]$type <- "negcontrol"
spec.subj.table.df[grepl("BLANK", spec.subj.table.df$SampleID),]$type <- "negcontrol"
spec.subj.table.df[grepl("fail", spec.subj.table.df$SampleID),]$type <- "failed"
spec.subj.table.df[grepl("noproduct", spec.subj.table.df$LogNo),]$type <- "failed"
spec.subj.table.df[grepl("BNE018v01", spec.subj.table.df$SampleID),]$type <- "misc"
# spec.subj.table.df[!grepl("VI", spec.subj.table.df$specimen),]$type <- "negcontrol"
spec.subj.table.df[grepl("CAS", spec.subj.table.df$SampleID),]$type <- "CAS"
dim(spec.subj.table.df[which(spec.subj.table.df$type == "COAST"),]) # 3694 samples, 103 features

with(spec.subj.table.df, table(visittype, type, useNA="ifany"))

with(spec.subj.table.df[which(is.na(spec.subj.table.df$visittype)),], 
		 data.frame(SampleID, LogNo, Description, type))

# Rearrange columns
refcols <- c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description", "RunName", "LogNo")
spec.subj.table.df <- spec.subj.table.df[, 
																				 c(refcols, 
																				 	setdiff(colnames(spec.subj.table.df), refcols))]

# Age in days
colnames(spec.subj.table.df)[colnames(spec.subj.table.df) == "age.days"] <-
	"age.day"

# Age in years
spec.subj.table.df$age.year <- ceiling((spec.subj.table.df$age.day)/365)

# Date of sampling, by days from first day of 1999
spec.subj.table.df$date.day <- as.integer(difftime(
	as.POSIXct(as.Date(spec.subj.table.df$colldt, format="%m/%d/%y")),
	as.POSIXct(as.Date("1999-01-01")),
	units="days"))

# Date of sampling, by years from 1999
spec.subj.table.df$date.year <- factor(
	as.numeric(ceiling(spec.subj.table.df$date.day/365)))

# Fix neg control LogNo so that entries are unique for each
# ---------------------------------------
spec.subj.table.df$RunName2 <- substr(spec.subj.table.df$RunName, 10, 12)
with(spec.subj.table.df[which(spec.subj.table.df$type == "negcontrol"),], 
		 data.frame(SampleID, LogNo, RunName))
sampleids <- spec.subj.table.df[which(spec.subj.table.df$type == "negcontrol"),]$SampleID
runname2 <- spec.subj.table.df[which(spec.subj.table.df$type == "negcontrol"),]$RunName2
spec.subj.table.df[which(spec.subj.table.df$type == "negcontrol"),]$LogNo <- 
	paste0(substr(sampleids,1, nchar(sampleids) - 6),"_",runname2)
with(spec.subj.table.df[which(spec.subj.table.df$type == "negcontrol"),], 
		 length(unique(LogNo)))
with(spec.subj.table.df[which(spec.subj.table.df$type == "negcontrol"),], 
		 length(LogNo))

# Well/mild/sick relabelling
# ---------------------------------------
with(spec.subj.table.df, table(visittype, Category, useNA="ifany"))
# Sick samples: 
# 	33 with score < 5 --> treat as mild; 
# 	4 with score NA --> remove
# Well/Sick: 
# 	159 with score >= 5 --> reclassify as Sick
# 	11 with score NA --> remove
# 	1 with core 0 --> treat as well
# 	Rest (360): treat as mild

spec.subj.table.df$visittype2 <- spec.subj.table.df$visittype
spec.subj.table.df$visittype2[which(spec.subj.table.df$visittype == "Sick" & 
													spec.subj.table.df$Category == "less5")] <- "Mild"
spec.subj.table.df$visittype2[which(spec.subj.table.df$visittype == "Sick" & 
													is.na(spec.subj.table.df$Category))] <- "remove"
spec.subj.table.df$visittype2[which(spec.subj.table.df$visittype == "Well/Sick")] <- "Mild"
spec.subj.table.df$visittype2[which(spec.subj.table.df$visittype == "Well/Sick" & 
													is.na(spec.subj.table.df$Category))] <- "remove"
spec.subj.table.df$visittype2[which(spec.subj.table.df$visittype == "Well/Sick" & 
													spec.subj.table.df$Category == 0)] <- "Well"
spec.subj.table.df$visittype2[which(spec.subj.table.df$visittype == "Well/Sick" & 
													spec.subj.table.df$Category == "more5")] <- "Sick"

# Mark healthy samples within 4 weeks after illness
# ---------------------------------------
sample.age.df <- with(spec.subj.table.df, 
											data.frame(subj, LogNo, visittype2, age.day))
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
	c <- which(sample.age.df.ind$visittype2 == "Well")
	k <- which(sample.age.df.ind$visittype2 == "Sick" | sample.age.df.ind$visittype2 == "Mild")
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
length(bad.logs) # 85

# Get logs of all good samples
not.bad.logs <- sample.age.df$LogNo[which(sample.age.df$bad == 0)] 
length(not.bad.logs) # 3851

spec.subj.table.df$healthy.not.within.4wk.inf <- 0
spec.subj.table.df$healthy.not.within.4wk.inf[
	which(spec.subj.table.df$LogNo %in% not.bad.logs)] <- 1

# 85 bad vs. 3851 good

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
coast.biom <- read_phyloseq(otu.file = "data/feature-table-taxa.biom", 
														metadata.file = "data/metadata_full.csv",
														taxonomy.file = NULL, 
														type = "biom")

coast.biom # 12464 taxa, 3915 samples, 171 features

# Samples missing microbiome
miss.micro <- mixedsort(setdiff(unique(spec.subj.table.df$SampleID), unique(row.names(sample_data(coast.biom)))))
miss.micro 
with(spec.subj.table.df[spec.subj.table.df$SampleID %in% miss.micro,], 
		 data.frame(SampleID, RunName, LogNo, visittype, Category, type))

# 9 COAST "samples" from 8 sample logs (two samples from log 706) produced very 
# few reads (in single digits; on inspection of imported qiime files for each 
# sample) --> after dada/filtering, 0 reads

save.image("output/2018-01-29_import_biom_format.Rdata")

