# ==============================================================================
# 2018-03-09 - Time before and after infection
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
											"reshape2", 
											"psych")

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

# Functions
# ---------------------------------------
# A better signif
good.signif <- function(x, n) {
	trimws(formatC(signif(x,digits=n), digits=n,format="fg"))
}
# Note that this converts number to character string

theme_set(theme_bw())

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("output/2018-02-20_micro_assoc_analysis_vir.Rdata")

SUBJ_ID.list <- unique(log.table.df2$SUBJ_ID) # 244

log.table.df2$healthy.no.prior.sick <- NA
log.table.df2$healthy.before.cat2 <- NA
log.table.df2$healthy.after.cat2 <- NA

# Double-check per heathy sample
# For each individual, get list of healthy samples
for (i in 1:length(SUBJ_ID.list)) {
	SUBJ_ID.target <- SUBJ_ID.list[i]
	healthy.samples <- log.table.df2$LogNo[which(log.table.df2$SUBJ_ID == SUBJ_ID.target &
																					 	log.table.df2$well.sick == "Well")]
	
	if (!identical(healthy.samples, character(0))) {
		# For each healthy sample and individual, get sick samples following/preceding it
		# Classify into appropriate timepoint category
		# Overwrite if higher-priority label exists 
		for (j in 1:length(healthy.samples)) {
			healthy.target <- healthy.samples[j]
			healthy.target.date <- log.table.df2$age[which(log.table.df2$LogNo == healthy.target)]
			
			# Find all sick samples after
			sick.samples.after <- log.table.df2$LogNo[which(log.table.df2$SUBJ_ID == SUBJ_ID.target &
																										 	log.table.df2$well.sick == "Sick" &
																										 	log.table.df2$age > healthy.target.date)]
			# If there are no sick samples after, current well sample does not precede 
			# any sick samples, so category is NA
			if (identical(sick.samples.after, character(0))) {
				log.table.df2$healthy.before.cat2[which(log.table.df2$LogNo == healthy.target)] <- NA
			} else {
				# get date of sick sample closest to current well sample
				sick.samples.after.date <- log.table.df2$age[which(log.table.df2$LogNo %in% sick.samples.after)]
				min.sick.samples.after.date <- min(as.numeric(sick.samples.after.date),na.rm=T)
				min.sick.samples.after <- log.table.df2$LogNo[which(log.table.df2$SUBJ_ID == SUBJ_ID.target &
																													 	log.table.df2$well.sick == "Sick" &
																													 	log.table.df2$age == min.sick.samples.after.date)]
				# calculate difference
				difference <- min.sick.samples.after.date - healthy.target.date
				log.table.df2$healthy.before.cat2[which(log.table.df2$LogNo == healthy.target)] <- 
					ifelse(difference >= 0 & difference < 365/52*2, "0-2w before",
								 ifelse(difference >= 365/52*2 & difference < 365/12, "2-4w before",
								 			 ifelse(difference >= 365/12 & difference < 365/12*2, "1-2m before",
								 			 			 ifelse(difference >= 365/12*2 & difference < 365/12*4, "2-4m before",
								 			 			 			 ifelse(difference >= 365/12*4 & difference < 365/12*6, "4-6m before",
								 			 			 			 			 ifelse(difference >= 365/12*6 & difference < 365, "6-12m before",
								 			 			 			 			 			 ifelse(difference >= 365/126, "More than 12m before", NA)))))))
			}
			
			sick.samples.before <- log.table.df2$LogNo[which(log.table.df2$SUBJ_ID == SUBJ_ID.target &
																												log.table.df2$well.sick == "Sick" &
																												log.table.df2$age < healthy.target.date)]
			if (identical(sick.samples.before, character(0))) {
				log.table.df2$healthy.after.cat2[which(log.table.df2$LogNo == healthy.target)] <- "No prior sickness"
			} else {
				sick.samples.before.date <- log.table.df2$age[which(log.table.df2$LogNo %in% sick.samples.before)]
				max.sick.samples.before.date <- max(as.numeric(sick.samples.before.date),na.rm=T)
				max.sick.samples.before <- log.table.df2$LogNo[which(log.table.df2$SUBJ_ID == SUBJ_ID.target &
																													 	log.table.df2$well.sick == "Sick" &
																													 	log.table.df2$age == max.sick.samples.before.date)]
				difference <- healthy.target.date - max.sick.samples.before.date
				log.table.df2$healthy.after.cat2[which(log.table.df2$LogNo == healthy.target)] <- 
					ifelse(difference > 0 & difference < 365/12*2, "1-2m after",
								 ifelse(difference >= 365/12*2 & difference < 365/12*4, "2-4m after",
								 			 ifelse(difference >= 365/12*4 & difference < 365/12*6, "4-6m after",
								 			 			 ifelse(difference >= 365/12*6 & difference < 365, "6-12m after",
								 			 			 			 ifelse(difference >= 365, "More than 12m after", NA)))))
			}
		}
	}
}

# Troubleshooting
with(log.table.df2, table(healthy.before.cat, healthy.before.cat2, useNA="ifany"))
with(log.table.df2, table(healthy.after.cat, healthy.after.cat2, useNA="ifany"))
with(log.table.df2[which(log.table.df2$healthy.before.cat == "0-1w before" & is.na(log.table.df2$healthy.before.cat2)),],
		 data.frame(SUBJ_ID, LogNo, age, well.sick, healthy.before.cat, healthy.before.cat2))

with(log.table.df2[which(log.table.df2$SUBJ_ID == 83 & log.table.df2$well.sick == "Well"),],
		 data.frame(SUBJ_ID, LogNo, age, well.sick, healthy.before.cat, healthy.before.cat2))

with(log.table.df2[which(log.table.df2$SUBJ_ID == 83 & log.table.df2$well.sick == "Sick"),],
		 data.frame(SUBJ_ID, LogNo, age, well.sick, healthy.before.cat, healthy.before.cat2))

# Get proportions
# levels(log.table.df2$MPG.name.top)
# levels(log.table.df2$MPG.colour.hex.top)
MPG.colour.hex.top.cut <- rev(c("#e31a1c", "#fdbf6f", "#ffd700", "#b29600"))

df1 <- log.table.df2 %>% 
	filter(!is.na(healthy.before.cat2)) %>%
	group_by(healthy.before.cat2, MPG.name.top) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(MPG.name.top %in% c("Moraxella MPG", "Streptococcus MPG1", 
														 "Haemophilus MPG1", "Haemophilus MPG2")) %>%
	mutate(MPG.name.top = factor(MPG.name.top, 
															 levels=rev(c("Moraxella MPG", "Streptococcus MPG1", 
															 				 "Haemophilus MPG1", "Haemophilus MPG2"))))
df1$healthy.before.cat2 <- 
	factor(df1$healthy.before.cat2,
				 levels=c("More than 12m before", "6-12m before", "4-6m before", 
				 				 "2-4m before", "1-2m before", "2-4w before", "0-2w before"))
df1$prop.lower.se <- with(df1, prop - qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df1$prop.lower.se <- ifelse(df1$prop.lower.se < 0, 0, df1$prop.lower.se)
df1$prop.upper.se <- with(df1, prop + qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df1$prop.upper.se <- ifelse(df1$prop.upper.se > 1, 1, df1$prop.upper.se)

ggplot(df1, aes(x=healthy.before.cat2, y=prop, fill=MPG.name.top)) + geom_col() +
	scale_fill_manual(values = MPG.colour.hex.top.cut)

df2 <- log.table.df2 %>% 
	filter(!is.na(healthy.after.cat2)) %>%
	group_by(healthy.after.cat2, MPG.name.top) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(MPG.name.top %in% c("Moraxella MPG", "Streptococcus MPG1", 
														 "Haemophilus MPG1", "Haemophilus MPG2")) %>%
	mutate(MPG.name.top = factor(MPG.name.top, 
															 levels=rev(c("Moraxella MPG", "Streptococcus MPG1", 
															 						 "Haemophilus MPG1", "Haemophilus MPG2"))))

df2$healthy.after.cat2 <- 
	factor(df2$healthy.after.cat2,
				 levels=c("1-2m after", "2-4m after", "4-6m after", 
				 				 "6-12m after", "More than 12m after", "No prior sickness"))
df2$prop.lower.se <- with(df2, prop - qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df2$prop.lower.se <- ifelse(df2$prop.lower.se < 0, 0, df2$prop.lower.se)
df2$prop.upper.se <- with(df2, prop + qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df2$prop.upper.se <- ifelse(df2$prop.upper.se > 1, 1, df2$prop.upper.se)

ggplot(df2, aes(x=healthy.after.cat2, y=prop, fill=MPG.name.top)) + geom_col() +
	scale_fill_manual(values = MPG.colour.hex.top.cut)

df3 <- log.table.df2 %>%
	filter(well.sick == "Sick") %>%
	group_by(MPG.name.top) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(MPG.name.top %in% c("Moraxella MPG", "Streptococcus MPG1", 
														 "Haemophilus MPG1", "Haemophilus MPG2")) %>%
	mutate(MPG.name.top = factor(MPG.name.top, 
															 levels=rev(c("Moraxella MPG", "Streptococcus MPG1", 
															 						 "Haemophilus MPG1", "Haemophilus MPG2"))))

df3$prop.lower.se <- with(df3, prop - qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df3$prop.lower.se <- ifelse(df3$prop.lower.se < 0, 0, df3$prop.lower.se)
df3$prop.upper.se <- with(df3, prop + qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df3$prop.upper.se <- ifelse(df3$prop.upper.se > 1, 1, df3$prop.upper.se)

ggplot(df3, aes(x="Sick", y=prop, fill=MPG.name.top)) + geom_col() +
	scale_fill_manual(values = MPG.colour.hex.top.cut)

df1 <- data.frame(df1 %>% rename(timepoint = healthy.before.cat2))
df2 <- data.frame(df2 %>% rename(timepoint = healthy.after.cat2))
df3 <- data.frame(df3 %>% mutate(timepoint = "Sick"))
df3 <- df3[,c(6,1:5)]

df <- rbind(df1, df3, df2)

# Display error bars for moraxella only

df4 <- log.table.df2 %>%
	filter(well.sick == "Well") %>%
	group_by(MPG.name.top) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(MPG.name.top %in% c("Moraxella MPG", "Streptococcus MPG1", 
														 "Haemophilus MPG1", "Haemophilus MPG2")) %>%
	mutate(MPG.name.top = factor(MPG.name.top, 
															 levels=rev(c("Moraxella MPG", "Streptococcus MPG1", 
															 						 "Haemophilus MPG1", "Haemophilus MPG2"))))
df4 <- df4 %>% rename(prop.healthy.mean = prop) %>%
	mutate(count = NULL)
df <- merge(df, df4, by="MPG.name.top", all=T)

df$facet <- ifelse(as.numeric(df$timepoint) %in% c(1:7), "Before",
										ifelse(as.numeric(df$timepoint) == 8, "Sick", 
													 ifelse(as.numeric(df$timepoint) %in% c(9:13), "After",
													 			 ifelse(as.numeric(df$timepoint) == 14, "Never", NA))))
df$facet <- factor(df$facet, 
										levels=c("Before", "Sick", "After", "Never"))

ggplot(df, aes(x=timepoint, y=prop, fill=MPG.name.top)) + geom_col() +
	geom_errorbar(aes(ymin=prop.lower.se, ymax=prop.upper.se), width=0.1) +
	scale_fill_manual(values = MPG.colour.hex.top.cut) + 
	geom_hline(aes(yintercept = prop.healthy.mean), linetype = 2) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
	facet_grid(MPG.name.top~facet, scales="free_x",space="free_x")

df0 <- df
df0$prop.lower.se[which(df0$MPG.name.top != "Moraxella MPG")] <- NA
df0$prop.upper.se[which(df0$MPG.name.top != "Moraxella MPG")] <- NA


ggplot(df0, aes(x=timepoint, y=prop, fill=MPG.name.top)) + geom_col() +
	geom_errorbar(aes(ymin=prop.lower.se, ymax=prop.upper.se), width=0.1) +
	scale_fill_manual(values = MPG.colour.hex.top.cut) + 
	geom_hline(yintercept = 
						 	df4$prop.healthy.mean[which(df4$MPG.name.top == "Moraxella MPG")], 
						 linetype = 2) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
	facet_grid(.~facet, scales="free_x",space="free_x")

# Virus

df1.vir <- log.table.df2 %>% 
	filter(!is.na(healthy.before.cat2), !is.na(No.vir)) %>%
	group_by(healthy.before.cat2, No.vir) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(No.vir == 0) 

df1.vir$healthy.before.cat2 <- 
	factor(df1.vir$healthy.before.cat2,
				 levels=c("More than 12m before", "6-12m before", "4-6m before", 
				 				 "2-4m before", "1-2m before", "2-4w before", "0-2w before"))
df1.vir$prop.lower.se <- with(df1.vir, prop - qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df1.vir$prop.lower.se <- ifelse(df1.vir$prop.lower.se < 0, 0, df1.vir$prop.lower.se)
df1.vir$prop.upper.se <- with(df1.vir, prop + qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df1.vir$prop.upper.se <- ifelse(df1.vir$prop.upper.se > 1, 1, df1.vir$prop.upper.se)

df2.vir <- log.table.df2 %>% 
	filter(!is.na(healthy.after.cat2), !is.na(No.vir)) %>%
	group_by(healthy.after.cat2, No.vir) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(No.vir == 0) 

df2.vir$healthy.after.cat2 <- 
	factor(df2.vir$healthy.after.cat2,
				 levels=c("1-2m after", "2-4m after", "4-6m after", 
				 				 "6-12m after", "More than 12m after", "No prior sickness"))
df2.vir$prop.lower.se <- with(df2.vir, prop - qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df2.vir$prop.lower.se <- ifelse(df2.vir$prop.lower.se < 0, 0, df2.vir$prop.lower.se)
df2.vir$prop.upper.se <- with(df2.vir, prop + qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df2.vir$prop.upper.se <- ifelse(df2.vir$prop.upper.se > 1, 1, df2.vir$prop.upper.se)

df3.vir <- log.table.df2 %>%
	filter(well.sick == "Sick", !is.na(No.vir)) %>%
	group_by(No.vir) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(No.vir == 0) 

df3.vir$prop.lower.se <- with(df3.vir, prop - qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df3.vir$prop.lower.se <- ifelse(df3.vir$prop.lower.se < 0, 0, df3.vir$prop.lower.se)
df3.vir$prop.upper.se <- with(df3.vir, prop + qnorm(0.8413)*sqrt((prop * (1-prop))/count))
df3.vir$prop.upper.se <- ifelse(df3.vir$prop.upper.se > 1, 1, df3.vir$prop.upper.se)

df1.vir <- data.frame(df1.vir %>% rename(timepoint = healthy.before.cat2))
df2.vir <- data.frame(df2.vir %>% rename(timepoint = healthy.after.cat2))
df3.vir <- data.frame(df3.vir %>% mutate(timepoint = "Sick"))
df3.vir <- df3.vir[,c(6,1:5)]

df.vir <- rbind(df1.vir, df3.vir, df2.vir)

df4.vir <- log.table.df2 %>%
	filter(well.sick == "Well", !is.na(No.vir)) %>%
	group_by(No.vir) %>%
	summarize(count = n()) %>%
	mutate(prop = count/sum(count)) %>%
	filter(No.vir == 0)

df4.vir <- df4.vir %>% rename(prop.healthy.mean = prop) %>%
	mutate(count = NULL)
df.vir <- merge(df.vir, df4.vir, by="No.vir", all=T)

df.vir$facet <- ifelse(as.numeric(df.vir$timepoint) %in% c(1:7), "Before",
									 ifelse(as.numeric(df.vir$timepoint) == 8, "Sick", 
									 			 ifelse(as.numeric(df.vir$timepoint) %in% c(9:13), "After",
									 			 			 ifelse(as.numeric(df.vir$timepoint) == 14, "Never", NA))))
df.vir$facet <- factor(df.vir$facet, 
									 levels=c("Before", "Sick", "After", "Never"))

ggplot(data=df0, aes(x=timepoint, y=prop)) + geom_col(data=df0, aes(fill=MPG.name.top)) +
	geom_errorbar(data=df0, aes(ymin=prop.lower.se, ymax=prop.upper.se), width=0.1) +
	scale_fill_manual(values = MPG.colour.hex.top.cut) + 
	geom_hline(data=df0, yintercept = 
						 	df4$prop.healthy.mean[which(df4$MPG.name.top == "Moraxella MPG")], 
						 linetype = 2) + 
	geom_point(data=df.vir, aes(x=timepoint, y=prop, group=facet, color="Virus"), shape=4) + 
	geom_line(data=df.vir, aes(x=timepoint, y=prop, group=facet, color="Virus")) +
	geom_errorbar(data=df.vir, aes(ymin=prop.lower.se, ymax=prop.upper.se, 
								color="Virus"), width=0.1) +
	scale_color_manual(values = "#bfbfbf") + 
	geom_hline(data=df.vir, color="#bfbfbf", yintercept = df4.vir$prop.healthy.mean, 
						 linetype = 2) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
	ylim(0,1) +
	facet_grid(.~facet, scales="free_x",space="free_x")
