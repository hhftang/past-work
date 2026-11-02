# Create metadata files for each run - COAST data
# ==============================================================================

rm(list=ls())

# load libraries
desired_packages <- c("reshape2", "dplyr", "tidyr", "Rcpp")
new_packages <- desired_packages[!(desired_packages %in% 
																	 	installed.packages()[,"Package"])]
if(length(new_packages)>0) install.packages(new_packages)
for (i in 1:length(desired_packages)){
	library(desired_packages[i], character.only=T)
}

# load functions
# Script to get reverse complementary strand for given DNA string
revComp.rcpp <- cppFunction(
	"std::string comp(std::string x) {
	const int n = x.length();
	for (int i=0; i < n; ++i) {
	if (x[i] == 'A' || x[i] == 'a')  x[i] = 'T';
	else if (x[i] == 'T' || x[i] == 't')  x[i] = 'A';
	else if (x[i] == 'G' || x[i] == 'g')  x[i] = 'C';
	else  x[i] = 'G';
	}
	std::reverse(x.begin(), x.end());
	return x;
	}")
revComp.rcpp.vec <- function(x) sapply(x, revComp.rcpp)

# Load metadata files
filenames <- list.files(path="data/metadata/",pattern="*.csv")
runnames <- substr(filenames,1,nchar(filenames)-4)

for (i in runnames){
	assign(paste0(i, ".df"), read.csv(file=paste0("data/metadata/", i, ".csv"),
				 stringsAsFactors=F))
}

# Load melbourne_by_sample
coast_master.df <- read.csv(file="data/melbourne_by_sample.csv",
														stringsAsFactors=F)
coast_master.df <- coast_master.df %>% 
	rename(DNA.. = dna,
				 Specimen.. = specimen)

# Merge metadata with melbourne_by_sample
# And rename headers to be in line with 

cols.in.front <- c("SampleID", "BarcodeSequence", "LinkerPrimerSequence",
											"Description", "RunName", "LogNo", "Category", "Age", 
											"Subjid")
	
for (i in runnames){
		assign(paste0(i, "_merge.df"), 
				 merge(get(paste0(i, ".df")), coast_master.df,
				 			by=c("DNA..","Specimen.."), 
				 			all.x=T) %>%
				 	rename(SampleID = DNA..,
				 				 BarcodeSequence = reverse.primer.barcode,
				 				 Description = Specimen..,
				 				 Category = visittype,
				 				 Age = age.days,
				 				 Subjid = subj) %>%
				 	mutate(LogNo = SampleID,
				 				 LinkerPrimerSequence = NA,
				 				 RunName = i,
				 				 BarcodeSequence = revComp.rcpp.vec(BarcodeSequence))
				 )
	df <- get(paste0(i, "_merge.df"))
	df$SampleID <-	gsub(" ", "", df$SampleID, fixed = TRUE)
	df$SampleID <-	gsub("_", "", df$SampleID, fixed = TRUE)
	assign(paste0(i, "_merge.df"), 
				 df[, c(cols.in.front, setdiff(
				 	colnames(df), cols.in.front))])
}

# Manual edits for COAST16S_002_merge.df
COAST16S_002_merge.df$SampleID[
	COAST16S_002_merge.df$SampleID == "NEGCON1"] <- "neg1extraction"
COAST16S_002_merge.df$SampleID[
	COAST16S_002_merge.df$SampleID == "NEGCON2"] <- "neg2extraction"
COAST16S_002_merge.df$SampleID[
	COAST16S_002_merge.df$Description == "neg amp #1"] <- "water1"
COAST16S_002_merge.df$SampleID[
	COAST16S_002_merge.df$Description == "neg amp #2"] <- "water2"
COAST16S_002_merge.df$SampleID <- 
	paste0(COAST16S_002_merge.df$SampleID, "_AULMM")

# Manual edits for COAST16S_003_merge.df
COAST16S_003_merge.df$SampleID <- 
	paste0(COAST16S_003_merge.df$SampleID, "_AVMDW")

# Manual edits for COAST16S_004_merge.df
COAST16S_004_merge.df$SampleID[
	COAST16S_004_merge.df$SampleID == "NEGCON7"] <- "NEG_CON_7"
COAST16S_004_merge.df$SampleID <- 
	paste0(COAST16S_004_merge.df$SampleID, "_ATLDY")

# Manual edits for COAST16S_005_merge.df
COAST16S_005_merge.df$SampleID[
	COAST16S_005_merge.df$SampleID == "negex5"] <- "neg_ex_5"
COAST16S_005_merge.df$SampleID <- 
	paste0(COAST16S_005_merge.df$SampleID, "_B2GRM")

# Manual edits for COAST16S_006_merge.df
COAST16S_006_merge.df$SampleID <- 
	paste0(COAST16S_006_merge.df$SampleID, "_AY558")

# Manual edits for COAST16S_007_merge.df
COAST16S_007_merge.df <- 
	COAST16S_007_merge.df[which(COAST16S_007_merge.df$COAST.plate == "COAST16S_007"),]
COAST16S_007_merge.df$SampleID[
	COAST16S_007_merge.df$SampleID == "2803expfailurenosampleinlibrary"] <- 
	"fail2803"
COAST16S_007_merge.df$SampleID[
	COAST16S_007_merge.df$SampleID == "2837expfailurenosampleinlibrary"] <- 
	"fail2837"
COAST16S_007_merge.df$SampleID[
	COAST16S_007_merge.df$SampleID == "2871expfailurenosampleinlibrary"] <- 
	"fail2871"
COAST16S_007_merge.df$SampleID[
	COAST16S_007_merge.df$SampleID == "3143expfailurenosampleinlibrary"] <- 
	"fail3143"
COAST16S_007_merge.df$SampleID[
	COAST16S_007_merge.df$SampleID == "485samplemissingusedwaterinstead"] <- 
	"water485"
COAST16S_007_merge.df$SampleID[
	COAST16S_007_merge.df$BarcodeSequence == "ACACCTGGTGAT"] <- "water1"
COAST16S_007_merge.df$SampleID[
	COAST16S_007_merge.df$BarcodeSequence == "GCACGACAACAC"] <- "water2"
COAST16S_007_merge.df$SampleID <- 
	paste0(COAST16S_007_merge.df$SampleID, "_AY5GK")

# Manual edits for COAST16S_008_merge.df
COAST16S_008_merge.df$SampleID[
	COAST16S_008_merge.df$SampleID == "negex#10"] <- "neg_ex_10"
COAST16S_008_merge.df$SampleID[
	COAST16S_008_merge.df$SampleID == "water"] <- "water_3"
COAST16S_008_merge.df$SampleID <- 
	paste0(COAST16S_008_merge.df$SampleID, "_AY5GW")

# Manual edits for COAST16S_009_merge.df
COAST16S_009_merge.df$SampleID <- 
	paste0(COAST16S_009_merge.df$SampleID, "_AY6DA")

# Manual edits for COAST16S_010_merge.df
COAST16S_010_merge.df$SampleID[
	COAST16S_010_merge.df$SampleID == "NEGEX6"] <- "NEG_EX_6"
COAST16S_010_merge.df$SampleID <- 
	paste0(COAST16S_010_merge.df$SampleID, "_AY3TN")

# Manual edits for COAST16S_011_merge.df
COAST16S_011_merge.df$SampleID[
	COAST16S_011_merge.df$SampleID == "negex11"] <- "neg_ex_11"
COAST16S_011_merge.df$SampleID[
	COAST16S_011_merge.df$SampleID == "water1"] <- "water_1"

COAST16S_011_1_merge.df <- COAST16S_011_merge.df
COAST16S_011_1_merge.df$SampleID <- 
	paste0(COAST16S_011_1_merge.df$SampleID, "_AYEC4")
COAST16S_011_2_merge.df <- COAST16S_011_merge.df
COAST16S_011_2_merge.df$SampleID <- 
	paste0(COAST16S_011_2_merge.df$SampleID, "_AYEE6") 

# Manual edits for COAST16S_012_merge.df
COAST16S_012_merge.df$SampleID[
	COAST16S_012_merge.df$SampleID == "negex12"] <- "neg_ex_12"
COAST16S_012_merge.df$SampleID[
	COAST16S_012_merge.df$SampleID == "water2"] <- "water_2"
COAST16S_012_merge.df$SampleID <- 
	paste0(COAST16S_012_merge.df$SampleID, "_AYEBY")

# Manual edits for COAST16S_013_merge.df
COAST16S_013_merge.df <- COAST16S_013_merge.df[-c(1:2),]
COAST16S_013_merge.df$SampleID <- 
	paste0(COAST16S_013_merge.df$SampleID, "_B36V3")

# Manual edits for COAST16S_014_merge.df
COAST16S_014_merge.df$SampleID <- 
	paste0(COAST16S_014_merge.df$SampleID, "_AWLY0")

# Manual edits for COAST16S_015_merge.df
COAST16S_015_merge.df$SampleID <- 
	paste0(COAST16S_015_merge.df$SampleID, "_AWLYV")

# Manual edits for COAST16S_016_merge.df
COAST16S_016_merge.df <- COAST16S_016_merge.df[-c(1:2),]
COAST16S_016_merge.df$SampleID <- 
	paste0(COAST16S_016_merge.df$SampleID, "_B37B9")

# Manual edits for COAST16S_016_repeat_merge.df
COAST16S_016_repeat_merge.df$SampleID <- 
	paste0(COAST16S_016_repeat_merge.df$SampleID, "_B5YK2") 

# Manual edits for COAST16S_017_merge.df
COAST16S_017_merge.df$SampleID[
	COAST16S_017_merge.df$SampleID == "379noproductseelibrary020"] <- "379see020"
COAST16S_017_merge.df$SampleID <- 
	paste0(COAST16S_017_merge.df$SampleID, "_B5NGN")

# Manual edits for COAST16S_018_merge.df
COAST16S_018_merge.df$SampleID <- 
	paste0(COAST16S_018_merge.df$SampleID, "_B4HKB")

# Manual edits for COAST16S_019_merge.df
COAST16S_019_merge.df$SampleID <- 
	paste0(COAST16S_019_merge.df$SampleID, "_B4PJ4")

# Manual edits for COAST16S_020_merge.df
COAST16S_020_merge.df$SampleID <- 
	paste0(COAST16S_020_merge.df$SampleID, "_B4PHL")

# Manual edits for COAST16S_021_merge.df
COAST16S_021_merge.df$SampleID <- 
	paste0(COAST16S_021_merge.df$SampleID, "_AVE8H")

# Manual edits for COAST16S_022_merge.df
COAST16S_022_merge.df$SampleID <- 
	paste0(COAST16S_022_merge.df$SampleID, "_B32F5")

# Manual edits for COAST16S_023_merge.df
COAST16S_023_merge.df$SampleID <- 
	paste0(COAST16S_023_merge.df$SampleID, "_B2MC5")

# Manual edits for COAST16S_024_merge.df
COAST16S_024_merge.df$SampleID <- 
	paste0(COAST16S_024_merge.df$SampleID, "_B47VB")

# Manual edits for COAST16S_025_merge.df
COAST16S_025_1_merge.df <- COAST16S_025_merge.df
COAST16S_025_1_merge.df$SampleID <- 
	paste0(COAST16S_025_1_merge.df$SampleID, "_B3BDW")
COAST16S_025_2_merge.df <- COAST16S_025_merge.df
COAST16S_025_2_merge.df$SampleID <- 
	paste0(COAST16S_025_2_merge.df$SampleID, "_B37PL") 

# Manual edits for COAST16S_026_merge.df
COAST16S_026_merge.df$SampleID <- 
	paste0(COAST16S_026_merge.df$SampleID, "_B36RV")

# Manual edits for COAST16S_027_merge.df
COAST16S_027_merge.df$SampleID <- 
	paste0(COAST16S_027_merge.df$SampleID, "_B47VK")

# Manual edits for COAST16S_028_merge.df
# Note misspelt "028balnk2"
COAST16S_028_merge.df$SampleID <- 
	paste0(COAST16S_028_merge.df$SampleID, "_B37JB")

# Manual edits for COAST16S_029_merge.df
COAST16S_029_merge.df$SampleID <- 
	paste0(COAST16S_029_merge.df$SampleID, "_B37D8")

# Manual edits for COAST16S_029_repeat_merge.df
COAST16S_029_repeat_merge.df$SampleID <- 
	paste0(COAST16S_029_repeat_merge.df$SampleID, "_B5NH5") 

# Manual edits for COAST16S_030_merge.df
# Note misspelt "030balnk1"
COAST16S_030_merge.df$SampleID <- 
	paste0(COAST16S_030_merge.df$SampleID, "_B3BG2")

# Manual edits for COAST16S_031_merge.df
COAST16S_031_merge.df$SampleID <- 
	paste0(COAST16S_031_merge.df$SampleID, "_AYB4A")

# Manual edits for COAST16S_032_merge.df
COAST16S_032_merge.df$SampleID <- 
	paste0(COAST16S_032_merge.df$SampleID, "_B3796")

# Manual edits for COAST16S_033_merge.df
COAST16S_033_merge.df$SampleID <- 
	paste0(COAST16S_033_merge.df$SampleID, "_B47VV")

# Manual edits for COAST16S_034_merge.df
COAST16S_034_merge.df$SampleID <- 
	paste0(COAST16S_034_merge.df$SampleID, "_B4HG9")

# Manual edits for COAST16S_035_merge.df
COAST16S_035_merge.df$SampleID <- 
	paste0(COAST16S_035_merge.df$SampleID, "_B4JBR")

# Manual edits for COAST16S_036_merge.df
COAST16S_036_merge.df$SampleID <- 
	paste0(COAST16S_036_merge.df$SampleID, "_B37B6")

# Manual edits for COAST16S_037_merge.df
COAST16S_037_merge.df$SampleID <- 
	paste0(COAST16S_037_merge.df$SampleID, "_B9H5H")

# Manual edits for COAST16S_037_repeat_merge.df
COAST16S_037_repeat_merge.df$SampleID <- 
	paste0(COAST16S_037_repeat_merge.df$SampleID, "_B69K2") 

# Output .tsv files
for (i in runnames){
	j <- substr(i,9,nchar(i))
	write.table(get(paste0(i, "_merge.df")), 
							file=paste0("data/metadata/AGRF_CAGRF13540",j,".tsv"),
							sep="\t",
							row.names=F,
							quote=F)
}

runnames_addit <- c("COAST16S_011_1", "COAST16S_011_2", 
							"COAST16S_025_1", "COAST16S_025_2")
for (i in runnames_addit){
	j <- substr(i,9,nchar(i))
	write.table(get(paste0(i, "_merge.df")), 
							file=paste0("data/metadata/AGRF_CAGRF13540",j,".tsv"),
							sep="\t",
							row.names=F,
							quote=F)
}

# Still need to add "#" to start of first column header ("#SampleID")

