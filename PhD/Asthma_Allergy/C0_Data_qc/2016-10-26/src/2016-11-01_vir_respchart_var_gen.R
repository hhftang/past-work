#' ---
#' title: "Generation of viral and respchart variables from CAS"
#' author: "Howard Tang"
#' date: "2016-11-01"
#' output: pdf_document
#' ---
#' 
#' This protocol loads all available CAS datasets relating to viral variables,
#' and the "respchart" dataset. 
#' 
#' It then generates variables in the dataframe "viral_sample_all.df" per 
#' log number ("viral_infect.df" for infection-only samples). 
#' 
#' It merges viral_infect.df with respchart.df to create viral_resp.df. 
#' 
#' It checks for consistency between viral and respchart clinical
#' variables, before generating new ari phenotypes as per COAST collaboration,
#' and then counting all phenotype occurrences per individual 
#' ("viral_resp_sum.df"), taking into account age of exclusion (from immunol
#' dataset).
#' 
#' Finally it checks for consistency with immunol dataset
#' 
#' NOTE: key differences from older scripts:
#' -  no longer applies hrv fix - johnson lab result "HRV" and gern lab results 
#' "HRVA/B/C" are now kept separate, as per Kat's request
#' - new variable types created for year 1 = 425d (labelled "_425")
#' 
#' #' Libraries
#' ----
#'
desired_packages <- c("compare", "dplyr", "ggplot2", "reshape2", "lazyeval")
new_packages <- desired_packages[!(desired_packages %in% 
																	 	installed.packages()[,"Package"])]
if(length(new_packages)>0) install.packages(new_packages)
for (i in 1:length(desired_packages)){
	library(desired_packages[i], character.only=T)
}
#'
#' Internal functions
#' ----
#' None
#' 
#' Viral and infection data
#' ----
#' 

# Load viral and infection tables
# ----
# Viral pathogens, per NPA sample, year 1
virus_y1.df <- read.csv("data/2015-12-18_NPA_Virus_y1.csv",header=T,sep=",")
# Viral pathogens, per NPA sample, years 2 and 3
virus_y23.df <- read.csv("data/2015-12-18_NPA_Virus_y23.csv",header=T,sep=",")
virus_y2.df <- virus_y23.df[virus_y23.df$age>365 & virus_y23.df$age<=730,]
virus_y3.df <- virus_y23.df[virus_y23.df$age>730 & virus_y23.df$age<=1095,]

# HRV subtyping, per NPA sample, year 1
virus_y1_hrvtype.df <- read.csv("data/2015-12-18_NPA_Virus_y1_hrvtype.csv",
																header=T, sep=",")
# HRV subtyping, per NPA sample, years 2 and 3
virus_y23_hrvtype.df <- read.csv("data/2015-12-18_NPA_Virus_y23_hrvtype.csv",
																 header=T,sep=",")
virus_y2_hrvtype.df <- virus_y23_hrvtype.df[
	virus_y23_hrvtype.df$age>365 & virus_y23_hrvtype.df$age<=730,]
virus_y3_hrvtype.df <- virus_y23_hrvtype.df[
	virus_y23_hrvtype.df$age>730 & virus_y23_hrvtype.df$age<=1095,]

# Clinical features of infection, per sample, years 1,2,3
ari.df <- read.table("data/CAS_PNA_cyrobanked_171013.txt",header=T,sep="\t")
colnames(ari.df) <- c("subjid","DOB","sex","Specimen","Category",
											"Year","age","col.date","LogNo",
											"1st.vial","2nd.vial","infection.log.no",
											"uri","lri","flri","wlri","wflri","slri","mlri")
ari.df <- with(ari.df, data.frame(subjid,sex,Category,age,LogNo,uri,lri,flri,
																	wlri,wflri))
ari_y1.df <- ari.df[ari.df$age<=365,]
colnames(ari_y1.df)[colnames(ari_y1.df)=="LogNo"] <- "LogNoS"
ari_y2.df <- ari.df[ari.df$age>365 & ari.df$age<=730,]
colnames(ari_y2.df)[colnames(ari_y2.df)=="LogNo"] <- "LogNoG"
ari_y3.df <- ari.df[ari.df$age>730 & ari.df$age<=1095,]
colnames(ari_y3.df)[colnames(ari_y3.df)=="LogNo"] <- "LogNoG"

# Merge year 1 infection tables; fix hrv typing
# ----
# (previously hrv column had some 0 entries when clearly hrva/b/c was 
# detected in that sample)
virus_y1_typed.df <- merge(virus_y1.df,virus_y1_hrvtype.df,
													 by=c("LogNoS","subjid","Category"),all=T)
# note: virus_y1_hrvtype.df does not have age variable
virus_y1_merged.df <- merge(virus_y1_typed.df, ari_y1.df,
														by=c("LogNoS","subjid","Category","age"),all=T)

# Merge year 2 infection tables
# ----
virus_y2_typed.df <- merge(virus_y2.df, virus_y2_hrvtype.df,
													 by=c("LogNoG","subjid","Category","age"),all=T)
virus_y2_merged.df <- merge(virus_y2_typed.df, ari_y2.df,
														by=c("LogNoG","subjid","Category","age"),all=T)

# Merge year 3 infection tables
# ----
virus_y3_typed.df <- merge(virus_y3.df, virus_y3_hrvtype.df,
													 by=c("LogNoG","subjid","Category","age"),all=T)
virus_y3_merged.df <- merge(virus_y3_typed.df, ari_y3.df,
														by=c("LogNoG","subjid","Category","age"),all=T)

# Trim to only those pathogens and clinical features of interest to us
# ----
virus_y1_trim.df <- with(virus_y1_merged.df, 
												 data.frame(LogNoS, subjid, Category, age, 
												 					 RSV, Influ, HRV, HRVA, HRVB, HRVC, 
												 					 uri, lri, wlri, flri, wflri))
virus_y2_trim.df <- with(virus_y2_merged.df, 
												 data.frame(LogNoG, subjid, Category, age, 
												 					 RSV, Influ, HRV, HRVA, HRVB, HRVC, 
												 					 uri, lri, wlri, flri, wflri))
virus_y3_trim.df <- with(virus_y3_merged.df, 
												 data.frame(LogNoG, subjid, Category, age, 
												 					 RSV, Influ, HRV, HRVA, HRVB, HRVC, 
												 					 uri, lri, wlri, flri, wflri))
colnames(virus_y1_trim.df)[colnames(virus_y1_trim.df)=="LogNoS"] <- "LogNo"
colnames(virus_y2_trim.df)[colnames(virus_y2_trim.df)=="LogNoG"]  <- "LogNo"
colnames(virus_y3_trim.df)[colnames(virus_y3_trim.df)=="LogNoG"]  <- "LogNo"

# Merge trimmed tables
# ----
virus_sample_all.df <- rbind(virus_y1_trim.df, 
														 virus_y2_trim.df, 
														 virus_y3_trim.df)

# Find and fix duplicated log numbers
# ----
dup_LogNo <- virus_sample_all.df$LogNo[duplicated(virus_sample_all.df$LogNo)==1]

# These are:
# 358   619   902  1832  1885   543  1933  2487 10671
#
# Interrogating what they are can be done by:
virus_sample_pre_dup_fix.df <- 
	virus_sample_all.df[virus_sample_all.df$LogNo %in% dup_LogNo,]
#
# or:
# virus_y1_typed.df[virus_y1_typed.df$LogNoS %in% dup_LogNo,]
# virus_y2_typed.df[virus_y2_typed.df$LogNoG %in% dup_LogNo,]
# virus_y3_typed.df[virus_y3_typed.df$LogNoG %in% dup_LogNo,]
# ari_y1.df[ari_y1.df$LogNoS %in% dup_LogNo,]
# ari_y1.df[ari_y1.df$LogNoG %in% dup_LogNo,]
# ari_y1.df[ari_y1.df$LogNoG %in% dup_LogNo,]

# These are semi-automatically fixed as follows:
useful_cols <- which(colnames(virus_sample_all.df) %in%
										 	c("age","RSV","Influ","HRV","HRVA","HRVB","HRVC",
										 		"uri","lri","wlri","flri","wflri"))

for (i in dup_LogNo) {
	# Identify row index of each duplicated LogNo
	x <- which(virus_sample_all.df$LogNo==i)
	
	# If more than 1 duplicate, break loop
	if(length(x) > 2) {
		paste("Error: need to resolve more than 1 duplicate for", i)
		break
	}
	
	# Subset out rows with that LogNo
	freq_var_test <- virus_sample_all.df[x,]
	
	# For each column/variable related to pathogen/infection/age:
	for (j in useful_cols) {
		# if both rows have NAs, value is NA
		# else, if at least one NA or zero, summate and treat NAs as zeroes
		if (all(is.na(virus_sample_all.df[x,j]))) {
			freq_var_test[,j] <- NA
		} else if (any(is.na(virus_sample_all.df[x,j]) |
									 virus_sample_all.df[x,j]==0)) {
			freq_var_test[,j] <- sum(virus_sample_all.df[x,j], na.rm=T)
		}
		# else, keep original value
	}
	
	# If infection reported --- need to change Category to Infectious
	if (with(freq_var_test, sum(uri, lri, wlri, flri, na.rm=T)>0)) {
		freq_var_test$Category <- "Infectious"
	} else {
		freq_var_test$Category <- "control"
	}
	virus_sample_all.df[x,] <- freq_var_test
}
# Tidy
rm(freq_var_test, useful_cols, i, j, x)

virus_sample_post_dup_fix.df <- 
	virus_sample_all.df[virus_sample_all.df$LogNo %in% dup_LogNo,]

# Note: Log 543 has one entry in virus_y1_typed.df, and the other in ari_y2.df

# Remove duplicates
virus_sample_all.df <- unique(virus_sample_all.df)

colnames(virus_sample_all.df)[colnames(virus_sample_all.df) == "RSV"] <- "rsv"
colnames(virus_sample_all.df)[colnames(virus_sample_all.df) == "Influ"] <- "influ"
colnames(virus_sample_all.df)[colnames(virus_sample_all.df) == "HRV"] <- "hrv"
colnames(virus_sample_all.df)[colnames(virus_sample_all.df) == "HRVA"] <- "hrva"
colnames(virus_sample_all.df)[colnames(virus_sample_all.df) == "HRVB"] <- "hrvb"
colnames(virus_sample_all.df)[colnames(virus_sample_all.df) == "HRVC"] <- "hrvc"
colnames(virus_sample_all.df)[colnames(virus_sample_all.df) == "subjid"] <- "SUBJ_ID"

# Restrict to infectious samples only
# ----
# We assume control samples with positive viral counts are spurious, for now.
# (although that sparks questions about whether any of the positive viral 
# counts in the infectious samples are also spurious) 
virus_infect.df <- virus_sample_all.df[
	which(virus_sample_all.df$Category == "Infectious"),]

# Exclude entries with NAs for clinical variables
# Trust respchart more than viral dataset
virus_infect_old.df <- virus_infect.df
virus_infect.df <- virus_infect.df[which(!is.na(virus_infect.df$uri)),]

#' Respchart dataset 
#' ----

# Load Respchart dataset
# ----
# Stick to old version of respchart, because we are no longer combining
# with diary files, and there are some bugged dates in respchart_dates.csv
respchart.df <- read.csv("data/respchart.csv",stringsAsFactors=T)

respchart.df$SUBJ_ID <- factor(respchart.df$SUBJ_ID)
colnames(respchart.df)[colnames(respchart.df)=="Collection.Date"] <- 
	"Date"  
colnames(respchart.df)[colnames(respchart.df)=="URI"] <- "uri"
colnames(respchart.df)[colnames(respchart.df)=="any.LRI"] <- "lri"
colnames(respchart.df)[colnames(respchart.df)=="Log.No"] <- "LogNo"

numeric_col <- colnames(respchart.df[, !names(respchart.df) %in% 
																		 	c("SUBJ_ID", "dob", "Date","DVISIT", 
																		 		"DPHONE", "INFN_S_ID","idin")])
for (i in length(numeric_col)){
	respchart.df[, i] <- as.numeric(as.character(respchart.df[, i]))
}

# All respchart entries are documented ARIs.

# Generate new definitions of uri, cough ari, fever ari, severe ari (fever|wheeze)
# severe lri ((fever|wheeze)+lri)
# ----
# uri_new = COUGH | RUNBLNOS
# cari = COUGH 
# fari = TEMP
# sari = TEMP | WHZRATT==1
# slri = (TEMP | WHZRATT==1) + (WRCHEST==1)

respchart.df$ari <- rep(1,nrow(respchart.df))
respchart.df$uri_new <- with(respchart.df, as.numeric((COUGH==1 | RUNBLNOS==1) &
																												(WRCHEST==2)))
respchart.df$cari <- with(respchart.df, as.numeric(COUGH==1))
respchart.df$fari <- with(respchart.df, as.numeric(TEMP==1))
respchart.df$sari <- with(respchart.df, as.numeric(TEMP==1 | WHZRATT==1))
# slri already made

# Quick checks
with(respchart.df, table(WRCHEST==2, uri, useNA="ifany"))
with(respchart.df, table(WRCHEST==1, lri, useNA="ifany"))
with(respchart.df, table((WHZRATT==1 | WHZRATT==2), lri, useNA="ifany"))
with(respchart.df, table(WHZRATT==1, wlri, useNA="ifany"))
with(respchart.df, table((WRCHEST==1 & TEMP==1), flri, useNA="ifany"))
with(respchart.df, table((WHZRATT==1 & TEMP==1), wflri, useNA="ifany"))
with(respchart.df, table(lri==1 & (WHZRATT==1 | TEMP==1), slri, useNA="ifany")) 
with(respchart.df, table(COUGH==1, cari, useNA="ifany")) 
with(respchart.df, table(TEMP==1, fari, useNA="ifany")) 
# All are internally consistent

# Fix dates and get ages
# ----
# Fix missing dates of collection, using DVISIT
respchart.df$Date <- as.Date(respchart.df$Date, "%d/%m/%Y")
respchart.df$DVISIT <- as.character(respchart.df$DVISIT)
index <- which(nchar(respchart.df$DVISIT)<8) 
respchart.df[index,]$DVISIT <- paste0("0", respchart.df[index,]$DVISIT)
respchart.df$DVISIT <- as.Date(as.character(respchart.df$DVISIT), "%d%m%Y")
index <- which(is.na(respchart.df$Date)) 
respchart.df[index,]$Date <- respchart.df[index,]$DVISIT
respchart.df[which(is.na(respchart.df$Date)),]

# Get ages
respchart.df$dob <- as.Date(respchart.df$dob, "%d/%m/%Y")
respchart.df$age <- with(respchart.df, Date - dob)

# Correct age vs. year
respchart.df$year_new <- with(respchart.df, as.numeric(ceiling(age/365)))
with(respchart.df, table(year, year_new))

#' Merge viral and respchart datasets
#' ----
#' 

virus_resp.df <- merge(virus_infect.df, respchart.df, by=c("SUBJ_ID", "LogNo"), all=T)

# Consistency checks
# ----
with(virus_resp.df[which(!is.na(virus_resp.df$age.x) & 
												 	!is.na(virus_resp.df$age.y)),], sum(age.x - age.y))
with(virus_resp.df, table(is.na(age.x), is.na(age.y)))
# Consistent, except there are 1022 entries in respchart.df missing from virus_infect.df

# Use respchart ages as they are more complete (no NAs)
colnames(virus_resp.df)[colnames(virus_resp.df)=="age.y"] <- "age"
virus_resp.df$age.x <- NULL

with(virus_resp.df, table(uri.x, uri.y, useNA="ifany"))
with(virus_resp.df, table(lri.x, lri.y, useNA="ifany"))
with(virus_resp.df, table(wlri.x, wlri.y, useNA="ifany"))
with(virus_resp.df, table(flri.x, flri.y, useNA="ifany"))
# Consistent apart from the fact that many respchart entries did not have virus sampling.

with(virus_resp.df[which(virus_resp.df$age<=365),], 
		 table(uri.x, uri.y, useNA="ifany"))
# 7 entries had no viral sampling at age 1
with(virus_resp.df[which(virus_resp.df$age>365 & virus_resp.df$age<=365*2),], 
		 table(uri.x, uri.y, useNA="ifany"))
# 14 entries at age 2
with(virus_resp.df[which(virus_resp.df$age>365*2 & virus_resp.df$age<=365*3),], 
		 table(uri.x, uri.y, useNA="ifany"))
# 15 entries at age 3
with(virus_resp.df[which(virus_resp.df$age>365*3 & virus_resp.df$age<=365*4),], 
		 table(uri.x, uri.y, useNA="ifany"))
with(virus_resp.df[which(virus_resp.df$age>365*4 & virus_resp.df$age<=365*5),], 
		 table(uri.x, uri.y, useNA="ifany"))
# All entries at ages 4 and 5

with(virus_resp.df[which(virus_resp.df$age<=365*3),], table(uri.x, uri.y, useNA="ifany"))
with(virus_resp.df[which(virus_resp.df$age<=365*3),], table(lri.x, lri.y, useNA="ifany"))
with(virus_resp.df[which(virus_resp.df$age<=365*3),], table(wlri.x, wlri.y, useNA="ifany"))
with(virus_resp.df[which(virus_resp.df$age<=365*3),], table(flri.x, flri.y, useNA="ifany"))
with(virus_resp.df[which(virus_resp.df$age<=365*3),], table(wflri.x, wflri.y, useNA="ifany"))
# Total of 36 entries across all first 3 years of life are missing viral samples 

# Use respchart values instead of viral
colnames(virus_resp.df)[colnames(virus_resp.df)=="uri.y"] <- "uri"
colnames(virus_resp.df)[colnames(virus_resp.df)=="lri.y"] <- "lri"
colnames(virus_resp.df)[colnames(virus_resp.df)=="wlri.y"] <- "wlri"
colnames(virus_resp.df)[colnames(virus_resp.df)=="flri.y"] <- "flri"
colnames(virus_resp.df)[colnames(virus_resp.df)=="wflri.y"] <- "wflri"
virus_resp.df$uri.x <- NULL
virus_resp.df$lri.x <- NULL
virus_resp.df$wlri.x <- NULL
virus_resp.df$flri.x <- NULL
virus_resp.df$wflri.x <- NULL

#' Generate summary datasets per individual
#' ----
#' 

# Get exclusion ages (EXCL_AGE) from immunol dataset
# ----
load("output/2016-09-23_immuno_var_gen.Rdata")
excl_age.df <- with(immunol_5.df, data.frame(SUBJ_ID, EXCL_AGE))
excl_age.df[excl_age.df$EXCL_AGE==0,]$EXCL_AGE <- 365*6

# Merge with main dataset
virus_resp.df <- merge(virus_resp.df, excl_age.df, by="SUBJ_ID", all=T)

# Month-based summaries (0-2, 2-6, 6-12m)
# ----
monthly_var <- c("0_2m","2_6m","6_12m")
monthly_val <- c(0,2,6,12)
clin_var_to_count <- c("ari", "uri", "lri", "wlri", "flri", 
											 "uri_new", "cari", "fari", "sari", "slri", "wflri")
path_var_to_count <- c("rsv", "influ", "hrv", "hrva", "hrvb", "hrvc")

for (k in 1:3){
	subset_df <- virus_resp.df %>%
		# Filter samples by desired timeperiod, and by exclusion age greater than  
		# timeperiod (if less, then value for timeperiod cannot be calculated)
		filter(age>365/12*(monthly_val[k]) &
					 	age<=365/12*(monthly_val[k+1]) &
					 	EXCL_AGE>365/12*(monthly_val[k+1]))
	df <- subset_df %>%
		# select only desirable variables
		select(SUBJ_ID, one_of(clin_var_to_count, path_var_to_count)) %>%
		# prepare to sum by subject
		group_by(SUBJ_ID) %>%
		# summate all desirable variables
		summarise_each(funs(sum(., na.rm=T))) %>%
		# Rename with timepoint suffix
		rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_", monthly_var[k])))
	for (j in 1:length(path_var_to_count)){
		df2 <- subset_df %>%
			# Filter samples by positive finding of desired pathogen
			filter_(interp(~ var == 1, var = as.name(path_var_to_count[j]))) %>%
			select(SUBJ_ID, one_of(clin_var_to_count)) %>%
			group_by(SUBJ_ID) %>%
			summarise_each(funs(sum(., na.rm=T))) %>%
			# Rename with pathogen prefix, timepoint suffix
			rename_(.dots = setNames(names(.)[-1], 
			 											 paste0(path_var_to_count[j],"_",
			 											 			 names(.)[-1], "_", monthly_var[k])))
		# Incrementally merge with previous dataframe
		df <- merge(df, df2, by="SUBJ_ID", all=T)
		# All NAs replaced with zeroes at this stage
		df[is.na(df)] <- 0
	}
	# Give complete dataframe a useful name
	assign(paste0("virus_resp_sum_",monthly_var[k],".df"), df)
}

# Year-based summaries (1,2,3,4,5)
# ----
for (k in 1:5){
	subset_df <- virus_resp.df %>%
		filter(age>365*(k-1) &
					 	age<=365*k &
					 	EXCL_AGE>365*k)
	if (k <= 3){ 
		df <- subset_df %>%
			select(SUBJ_ID, one_of(clin_var_to_count, path_var_to_count)) %>%
			group_by(SUBJ_ID) %>%
			summarise_each(funs(sum(., na.rm=T))) %>%
			rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_", k)))
		for (j in 1:length(path_var_to_count)){
			df2 <- subset_df %>%
				filter_(interp(~ var == 1, var = as.name(path_var_to_count[j]))) %>%
				select(SUBJ_ID, one_of(clin_var_to_count)) %>%
				group_by(SUBJ_ID) %>%
				summarise_each(funs(sum(., na.rm=T))) %>%
				rename_(.dots = setNames(names(.)[-1], 
																 paste0(path_var_to_count[j],"_",
																 			 names(.)[-1], "_",k)))
			df <- merge(df, df2, by="SUBJ_ID", all=T)
			df[is.na(df)] <- 0
		}
	} else {
		df <- subset_df %>%
			select(SUBJ_ID, one_of(clin_var_to_count)) %>%
			group_by(SUBJ_ID) %>%
			summarise_each(funs(sum(., na.rm=T))) %>%
			rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_", k)))
	}
	assign(paste0("virus_resp_sum_y",k,".df"), df)
}

# Year 1 = 425d
# ----

subset_df <- virus_resp.df %>%
	filter(age<=425 & EXCL_AGE>425)
df <- subset_df %>%
	select(SUBJ_ID, one_of(clin_var_to_count, path_var_to_count)) %>%
	group_by(SUBJ_ID) %>%
	summarise_each(funs(sum(., na.rm=T))) %>%
	rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_425")))
for (j in 1:length(path_var_to_count)){
	df2 <- subset_df %>%
		filter_(interp(~ var == 1, var = as.name(path_var_to_count[j]))) %>%
		select(SUBJ_ID, one_of(clin_var_to_count)) %>%
		group_by(SUBJ_ID) %>%
		summarise_each(funs(sum(., na.rm=T))) %>%
		rename_(.dots = setNames(names(.)[-1], 
														 paste0(path_var_to_count[j],"_",
														 			 names(.)[-1], "_425")))
	df <- merge(df, df2, by="SUBJ_ID", all=T)
	df[is.na(df)] <- 0
}
virus_resp_sum_y425.df <- df

# Exclusion ages
# ----
for (k in 1:3){
	df <- eval(parse(text=paste0("virus_resp_sum_",monthly_var[k],".df")))
	df <- merge(df, excl_age.df, by="SUBJ_ID", all=T)
	for (i in 1:nrow(df)) {
		if(df[i,]$EXCL_AGE <= 365/12*monthly_val[k+1]) {
			df[i,-which(colnames(df) %in% c("SUBJ_ID","EXCL_AGE"))] <- NA
		} else {
			df[i,which(is.na(df[i,]))] <- 0
		}
	}
	assign(paste0("virus_resp_sum_",monthly_var[k],".df"), df)
}

for (k in 1:5){
	df <- eval(parse(text=paste0("virus_resp_sum_y",k,".df")))
	df <- merge(df, excl_age.df, by="SUBJ_ID", all=T)
	for (i in 1:nrow(df)) {
		if(df[i,]$EXCL_AGE <=365*k) {
			df[i,-which(colnames(df) %in% c("SUBJ_ID","EXCL_AGE"))] <- NA
		} else {
			df[i,which(is.na(df[i,]))] <- 0
		}
	}
	assign(paste0("virus_resp_sum_y",k,".df"), df)
}

df <- virus_resp_sum_y425.df
df <- merge(df, excl_age.df, by="SUBJ_ID", all=T)
for (i in 1:nrow(df)) {
	if(df[i,]$EXCL_AGE <=425) {
		df[i,-which(colnames(df) %in% c("SUBJ_ID","EXCL_AGE"))] <- NA
	} else {
		df[i,which(is.na(df[i,]))] <- 0
	}
}
virus_resp_sum_y425.df <- df

# merge all
virus_resp_sum.df <- Reduce(function(...) 
	merge(..., by=c("SUBJ_ID", "EXCL_AGE"), all=T),
														list(virus_resp_sum_0_2m.df, virus_resp_sum_2_6m.df, 
																 virus_resp_sum_6_12m.df, virus_resp_sum_y1.df,
																 virus_resp_sum_y2.df, virus_resp_sum_y3.df,
																 virus_resp_sum_y4.df, virus_resp_sum_y5.df,
																 virus_resp_sum_y425.df))

# Check that the age exclusion actually worked
with(virus_resp_sum.df, 
		 data.frame(SUBJ_ID, EXCL_AGE, ari_0_2m, ari_2_6m, ari_6_12m, 
		 					 ari_1, ari_2, ari_3, ari_4, ari_5, ari_425))

# Summary statistics across all subjects
# ----

virus_resp_sum.df$EXCL_AGE <- NULL
stat_all.df <- virus_resp_sum.df[,-1] %>%
	summarise_each(funs(sum(.,na.rm=T), mean(.,na.rm=T), 
											median = quantile(.,probs=0.5, na.rm=T))) 
stat_long.df <- stat_all.df %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	as.matrix() %>%
	matrix(ncol=11, byrow=T) %>%
	signif(digits=3) %>%
	as.data.frame(stringsAsFactors=F) 

colnames(stat_long.df) <- c("ari","uri","lri","wlri","flri",
														"uri_new","cari","fari","sari","slri", "wflri")
stat_long.df$time <- rep(c("0-2m","2-6m","6-12m","1","2","3","4","5","425"),3)
stat_long_sum.df <- stat_long.df[1:9,]
stat_long_mean.df <- stat_long.df[10:18,]
stat_long_median.df <- stat_long.df[19:27,]
write.table(stat_long_sum.df, "output/2016-11-01_stat_long_sum.txt", 
						sep="\t")
write.table(stat_long_mean.df, "output/2016-11-01_stat_long_mean.txt", 
						sep="\t")

# Compare lri, wlri/wari, fari, flri, wfari/wflri (fever+wheeze), 
# sari (fever|wheeze) and slri ((fever|wheeze) + lri)
with(stat_long_sum.df, data.frame(time, lri, wlri, fari, flri, wflri, sari, slri))
with(stat_long_mean.df, data.frame(time, lri, wlri, fari, flri, wflri, sari, slri))

# Graphs
# ----

theme_set(theme_bw())

virus_resp_sum_uri.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("uri_new")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("uri_new_0_2m","uri_new_2_6m",
																		"uri_new_6_12m","uri_new_1","uri_new_2",
																		"uri_new_3","uri_new_4","uri_new_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_uri_hist.pdf",height=5,width=8)
virus_resp_sum_uri.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of URIs in time period") +
	labs(x="Number of infections", y="Count")
dev.off()

virus_resp_sum_lri.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("lri_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("lri_0_2m","lri_2_6m",
																		"lri_6_12m","lri_1","lri_2",
																		"lri_3","lri_4","lri_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_lri_hist.pdf",height=5,width=8)
virus_resp_sum_lri.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of LRIs in time period") +
	labs(x="Number of infections", y="Count") + xlim(c(-1,8))
dev.off()

virus_resp_sum_fari.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("fari_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("fari_0_2m","fari_2_6m",
																		"fari_6_12m","fari_1","fari_2",
																		"fari_3","fari_4","fari_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_fari_hist.pdf",height=5,width=8)
virus_resp_sum_fari.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of febrile RIs in time period") +
	labs(x="Number of infections", y="Count") + xlim(c(-1,11))
dev.off()

virus_resp_sum_cari.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("cari_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("cari_0_2m","cari_2_6m",
																		"cari_6_12m","cari_1","cari_2",
																		"cari_3","cari_4","cari_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_cari_hist.pdf",height=5,width=8)
virus_resp_sum_cari.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of coughing RIs in time period") +
	labs(x="Number of infections", y="Count")
dev.off()

virus_resp_sum_wlri.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("wlri_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("wlri_0_2m","wlri_2_6m",
																		"wlri_6_12m","wlri_1","wlri_2",
																		"wlri_3","wlri_4","wlri_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_wlri_hist.pdf",height=5,width=8)
virus_resp_sum_wlri.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of wheezing RIs in time period") +
	labs(x="Number of infections", y="Count")
dev.off()

# flri
virus_resp_sum_flri.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("flri_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("flri_0_2m","flri_2_6m",
																		"flri_6_12m","flri_1","flri_2",
																		"flri_3","flri_4","flri_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_flri_hist.pdf",height=5,width=8)
virus_resp_sum_flri.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of febrile LRIs in time period") +
	labs(x="Number of infections", y="Count")
dev.off()

# wflri
virus_resp_sum_wflri.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("wflri_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("wflri_0_2m","wflri_2_6m",
																		"wflri_6_12m","wflri_1","wflri_2",
																		"wflri_3","wflri_4","wflri_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_wflri_hist.pdf",height=5,width=8)
virus_resp_sum_wflri.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of wheezy and febrile LRIs in time period") +
	labs(x="Number of infections", y="Count")
dev.off()

# sari
virus_resp_sum_sari.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("sari_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("sari_0_2m","sari_2_6m",
																		"sari_6_12m","sari_1","sari_2",
																		"sari_3","sari_4","sari_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_sari_hist.pdf",height=5,width=8)
virus_resp_sum_sari.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of severe RIs (wheeze or fever) in time period") +
	labs(x="Number of infections", y="Count")
dev.off()

# slri
virus_resp_sum_slri.df <- virus_resp_sum.df %>%
	select(SUBJ_ID, starts_with("slri_")) %>%
	select(-matches("rsv|influ|hrv|hrva|hrvb|hrvc")) %>%
	melt(variable.name="age",value.name="frequency",id.vars="SUBJ_ID") %>%
	mutate(age = factor(age, levels=c("slri_0_2m","slri_2_6m",
																		"slri_6_12m","slri_1","slri_2",
																		"slri_3","slri_4","slri_5"), 
											labels=c("0-2m","2-6m","6-12m","0-1y","1-2y","2-3y",
															 "3-4y","4-5y")))
pdf("output/2016-11-01_slri_hist.pdf",height=5,width=8)
virus_resp_sum_slri.df %>%
	filter(age!="0-1y", age!="3-4y", age!="4-5y") %>%
	ggplot(aes(x=frequency)) + geom_histogram(binwidth=1) + facet_grid(.~age) +
	ggtitle("Histogram of individuals with number of severe LRIs (wheeze or fever AND LRI) in time period") +
	labs(x="Number of infections", y="Count")
dev.off()

# Consistency check with immunol dataset
# ----
premaster.df <- merge(virus_resp_sum.df, immunol_cyto_all.df, 
											by="SUBJ_ID", all=T)

# Fix some variable names
colnames(premaster.df)[colnames(premaster.df)=="rsv_2.x"] <- "rsv_2"
colnames(premaster.df)[colnames(premaster.df)=="rsv_2.y"] <- "rsv2"
colnames(premaster.df)[colnames(premaster.df)=="rsv_3.x"] <- "rsv_3"
colnames(premaster.df)[colnames(premaster.df)=="rsv_3.y"] <- "rsv3"

# check dob of respchart.df vs. immunol
nrow(unique(merge(with(premaster.df, data.frame(SUBJ_ID, dob)), 
			with(respchart.df, data.frame(SUBJ_ID, dob)),
			by=c("SUBJ_ID", "dob"), all=T)))
nrow(premaster.df[which(!is.na(premaster.df$dob)),])

# Check infection consistencies
with(premaster.df, table(ari_1, ari1, useNA = "ifany"))
with(premaster.df, table(ari_2, ari2, useNA = "ifany"))
with(premaster.df, table(ari_3, ari3, useNA = "ifany"))
with(premaster.df, table(ari_4, ari4, useNA = "ifany"))
with(premaster.df, table(ari_5, ari5, useNA = "ifany"))
with(premaster.df, table(uri_1, uri1, useNA = "ifany"))
with(premaster.df, table(uri_2, uri2, useNA = "ifany"))
with(premaster.df, table(uri_3, uri3, useNA = "ifany"))
with(premaster.df, table(uri_4, uri4, useNA = "ifany"))
with(premaster.df, table(uri_5, uri5, useNA = "ifany"))
with(premaster.df, table(lri_1, lri1, useNA = "ifany"))
with(premaster.df, table(lri_2, lri2, useNA = "ifany"))
with(premaster.df, table(lri_3, lri3, useNA = "ifany"))
with(premaster.df, table(lri_4, lri4, useNA = "ifany"))
with(premaster.df, table(lri_5, lri5, useNA = "ifany"))

# Manually interrogating each of the discrepancies where value from immunol
# dataset exceeds respchart, suggests that most of these are due to mislabelling 
# of "year" variable, leading to count discrepancies
premaster.df[which(premaster.df$ari_1 != premaster.df$ari1),]$SUBJ_ID

# Check viral pathogen consistencies
with(premaster.df, table(hrv_1, hrv1, useNA = "ifany"))
with(premaster.df, table(hrv_2, hrv2, useNA = "ifany"))
with(premaster.df, table(hrv_3, hrv3, useNA = "ifany"))
with(premaster.df, table(hrv_425, hrv1, useNA = "ifany"))

with(virus_resp.df[which(virus_resp.df$SUBJ_ID %in%
												 	premaster.df[which(premaster.df$hrv_1 != premaster.df$hrv1),]$SUBJ_ID
												 & virus_resp.df$year_new == 1 & virus_resp.df$hrv == 1),],
		 data.frame(SUBJ_ID, dob, Date, year, year_new, hrv, lri))
with(premaster.df[which(premaster.df$hrv_1 != premaster.df$hrv1),], data.frame(SUBJ_ID, hrv1, hrv_1))


with(premaster.df, table(rsv_1, rsv1, useNA = "ifany"))
with(premaster.df, table(rsv_2, rsv2, useNA = "ifany"))
with(premaster.df, table(rsv_3, rsv3, useNA = "ifany"))


with(virus_resp.df[which(virus_resp.df$SUBJ_ID %in%
												 	premaster.df[which(
												 		((premaster.df$rsv_1 != premaster.df$rsv1) | 
												 			is.na(premaster.df$rsv_1) |
												 			is.na(premaster.df$rsv1)) &
												 			!(is.na(premaster.df$rsv_1) & is.na(premaster.df$rsv1))
												 		),]$SUBJ_ID & 
												 	virus_resp.df$year_new == 1),],
		 data.frame(SUBJ_ID, dob, Date, year, rsv, lri))
with(premaster.df[which(
	((premaster.df$rsv_1 != premaster.df$rsv1) | 
	 	is.na(premaster.df$rsv_1) |
	 	is.na(premaster.df$rsv1)) &
		!(is.na(premaster.df$rsv_1) & is.na(premaster.df$rsv1))),],
		 data.frame(SUBJ_ID, EXCL_AGE, rsv1, rsv_1))
		 	

# Many discrepancies, especially at earlier years for hrv, and
# later years for rsv - unclear why counts are so different

with(premaster.df, table(hrv_lri_1, lrihrv1, useNA="ifany"))
with(premaster.df, table(hrv_lri_2, lrihrv2, useNA="ifany"))
with(premaster.df, table(hrv_lri_3, lrihrv3, useNA="ifany"))
with(premaster.df, table(hrv_lri_425, lrihrv1, useNA="ifany"))

with(virus_resp.df[which(virus_resp.df$SUBJ_ID %in%
												 	premaster.df[which(
												 		((premaster.df$hrv_lri_1 != premaster.df$lrihrv1) | 
												 		 	is.na(premaster.df$hrv_lri_1) |
												 		 	is.na(premaster.df$lrihrv1)) &
												 			!(is.na(premaster.df$hrv_lri_1) & is.na(premaster.df$lrihrv1))
												 	),]$SUBJ_ID & 
												 	virus_resp.df$year_new == 1),],
		 data.frame(SUBJ_ID, dob, Date, year, hrv, lri))
with(premaster.df[which(
	((premaster.df$hrv_lri_1 != premaster.df$lrihrv1) | 
	 	is.na(premaster.df$hrv_lri_1) |
	 	is.na(premaster.df$lrihrv1)) &
		!(is.na(premaster.df$hrv_lri_1) & is.na(premaster.df$lrihrv1))),],
	data.frame(SUBJ_ID, EXCL_AGE, hrv1, hrv_1, lri1, lri_1, lrihrv1, hrv_lri_1))


with(premaster.df, table(rsv_lri_1, lrirsv1, useNA="ifany"))
with(premaster.df, table(rsv_lri_2, lrirsv2, useNA="ifany"))
with(premaster.df, table(rsv_lri_3, lrirsv3, useNA="ifany"))
with(premaster.df, table(rsv_lri_425, lrirsv1, useNA="ifany"))

with(virus_resp.df[which(virus_resp.df$SUBJ_ID %in%
												 	premaster.df[which(
												 		((premaster.df$rsv_lri_1 != premaster.df$lrirsv1) | 
												 		 	is.na(premaster.df$rsv_lri_1) |
												 		 	is.na(premaster.df$lrirsv1)) &
												 			!(is.na(premaster.df$rsv_lri_1) & is.na(premaster.df$lrirsv1))
												 	),]$SUBJ_ID & 
												 	virus_resp.df$year_new == 1),],
		 data.frame(SUBJ_ID, dob, Date, year, rsv, lri))
with(premaster.df[which(
	((premaster.df$rsv_lri_1 != premaster.df$lrirsv1) | 
	 	is.na(premaster.df$rsv_lri_1) |
	 	is.na(premaster.df$lrirsv1)) &
		!(is.na(premaster.df$rsv_lri_1) & is.na(premaster.df$lrirsv1))),],
	data.frame(SUBJ_ID, EXCL_AGE, rsv1, rsv_1, lri1, lri_1, lrirsv1, rsv_lri_1))


with(virus_resp.df[which(virus_resp.df$SUBJ_ID %in%
												 	premaster.df[which(
												 		((premaster.df$influ_lri_1 != premaster.df$lriinflu1) | 
												 		 	is.na(premaster.df$influ_lri_1) |
												 		 	is.na(premaster.df$lriinflu1)) &
												 			!(is.na(premaster.df$influ_lri_1) & is.na(premaster.df$lriinflu1))
												 	),]$SUBJ_ID & 
												 	virus_resp.df$year_new == 1),],
		 data.frame(SUBJ_ID, dob, Date, year, influ, lri))
with(premaster.df[which(
	((premaster.df$influ_lri_1 != premaster.df$lriinflu1) | 
	 	is.na(premaster.df$influ_lri_1) |
	 	is.na(premaster.df$lriinflu1)) &
		!(is.na(premaster.df$influ_lri_1) & is.na(premaster.df$lriinflu1))),],
	data.frame(SUBJ_ID, EXCL_AGE, influ1, influ_1, lri1, lri_1, lriinflu1, influ_lri_1))



# Many discrepancies (rsv is better) - again unclear why

# Need to check consistency for semi-derived variables (number of 
# infections caused by known pathogen), based on values of overall number of
# infections (ari, uri, lri).
# ---
with(premaster.df, table(ari_1, rsv_1, useNA = "ifany"))
with(premaster.df, table(ari_1, hrv_1, useNA = "ifany"))
with(premaster.df, table(ari_1, influ_1, useNA = "ifany"))
with(premaster.df, table(ari_2, rsv_2, useNA = "ifany"))
with(premaster.df, table(ari_2, hrv_2, useNA = "ifany"))
with(premaster.df, table(ari_2, influ_2, useNA = "ifany"))
with(premaster.df, table(ari_3, rsv_3, useNA = "ifany"))
with(premaster.df, table(ari_3, hrv_3, useNA = "ifany"))
with(premaster.df, table(ari_3, influ_3, useNA = "ifany"))
with(premaster.df, table(uri_1, rsv_uri_1, useNA = "ifany"))
with(premaster.df, table(uri_1, hrv_uri_1, useNA = "ifany"))
with(premaster.df, table(uri_1, influ_uri_1, useNA = "ifany"))
with(premaster.df, table(uri_2, rsv_uri_2, useNA = "ifany"))
with(premaster.df, table(uri_2, hrv_uri_2, useNA = "ifany"))
with(premaster.df, table(uri_2, influ_uri_2, useNA = "ifany"))
with(premaster.df, table(uri_3, rsv_uri_3, useNA = "ifany"))
with(premaster.df, table(uri_3, hrv_uri_3, useNA = "ifany"))
with(premaster.df, table(uri_3, influ_uri_3, useNA = "ifany"))
with(premaster.df, table(lri_1, rsv_lri_1, useNA = "ifany"))
with(premaster.df, table(lri_1, hrv_lri_1, useNA = "ifany"))
with(premaster.df, table(lri_1, influ_lri_1, useNA = "ifany"))
with(premaster.df, table(lri_2, rsv_lri_2, useNA = "ifany"))
with(premaster.df, table(lri_2, hrv_lri_2, useNA = "ifany"))
with(premaster.df, table(lri_2, influ_lri_2, useNA = "ifany"))
with(premaster.df, table(lri_3, rsv_lri_3, useNA = "ifany"))
with(premaster.df, table(lri_3, hrv_lri_3, useNA = "ifany"))
with(premaster.df, table(lri_3, influ_lri_3, useNA = "ifany"))

save.image("output/2016-11-01_vir_respchart_var_gen.Rdata")

