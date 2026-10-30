# ==============================================================================
# 2018-02-20 - Microbiome distributions and association analyses (GEEs)
# VIRUSES
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

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("output/2018-02-20_micro_assoc_analysis.Rdata")


# Ggplot and graphical settings
# ---------------------------------------
theme_set(theme_bw())

seqPalette <- colorRampPalette(c("#fef0d9","#fdd49e","#fdbb84",
																 "#fc8d59","#e34a33","#b30000"))
quaPalette <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
								'#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
# quaPalette <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
# 								'#ff7f00','#ffff33','#a65628','#f781bf',
# 								'#999999')
# quaPalette <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4',
#								'#fed9a6','#ffffcc','#e5d8bd','#fddaec',
#								'#f2f2f2')
divPalette <- colorRampPalette(rev(c("#d73027","#fc8d59","#fee090",
																		 "#e0f3f8","#91bfdb","#4575b4")))
# ------------------------------------------------------------------------------
# VIRUS
# ------------------------------------------------------------------------------

# Create new variables
# -------------------------------------

# First, get source lab of samples

virus_y1_typed.df$source <- "Johnston"
virus_y2_typed.df$source <- "Gern"
virus_y3_typed.df$source <- "Gern"

virus_y1_typed.df2 <- virus_y1_typed.df %>%
	rename(LogNo = LogNoS, MPV = HMPV)

virus_y2_typed.df$Adeno <- with(virus_y2_typed.df, 
																 as.numeric(AdVB + AdVC + AdVE > 0))
virus_y2_typed.df$Corona <- with(virus_y2_typed.df, 
																 as.numeric(CoVOC43 + CoVNL63 > 0))
virus_y2_typed.df$Paraflu <- with(virus_y2_typed.df, 
																 as.numeric(PIV1 + PIV2 + PIV3 + PIV4b > 0))
virus_y3_typed.df$Adeno <- with(virus_y3_typed.df, 
																 as.numeric(AdVB + AdVC + AdVE > 0))
virus_y3_typed.df$Corona <- with(virus_y3_typed.df, 
																	as.numeric(CoVOC43 + CoVNL63 > 0))
virus_y3_typed.df$Paraflu <- with(virus_y3_typed.df, 
																	 as.numeric(PIV1 + PIV2 + PIV3 + PIV4b > 0))

virus_y2_typed.df2 <- virus_y2_typed.df %>%
	rename(LogNo = LogNoG)
virus_y3_typed.df2 <- virus_y3_typed.df %>%
	rename(LogNo = LogNoG)

virus.all.typed.df <- merge(virus_y1_typed.df2, virus_y2_typed.df2,
			by = c("LogNo", "subjid", "Category", "age", "Adeno", "Corona", "HRV",
						 "HRVA", "HRVB", "HRVC", "Influ", "MPV", "Paraflu", "RSV", "source"),
			all=T)

virus.all.typed.df <- merge(virus.all.typed.df, virus_y3_typed.df2,
			by = c("LogNo", "subjid", "Category", "age", 
						 "Adeno", "AdVB", "AdVC", "AdVE", 
						 "BoV", "Corona", "CoVOC43", "CoVNL63", "EV", 
						 "HRV", "HRVA", "HRVB", "HRVC", 
						 "Influ", "FluA", "FluB", "MPV", 
						 "Paraflu", "PIV1", "PIV2", "PIV3", "PIV4b",
						 "RSV", "RSVA", "RSVB",
						 "source"),
			all=T)

virus.all.typed.df[which(duplicated(virus.all.typed.df$LogNo) |
														duplicated(virus.all.typed.df$LogNo, fromLast = T)),]

# Find and fix duplicated log numbers
# ----
# dup_LogNo <- virus.all.typed.df$LogNo[duplicated(virus.all.typed.df$LogNo)==1]

# There are none:
# # 358   619   902  1832  1885   543  1933  2487 10671
# #
# # Interrogating what they are can be done by:
# virus.all.merged.pre.dup.fix.df <- 
# 	virus.all.merged.df[virus.all.merged.df$LogNo %in% dup_LogNo,]
# 
# # or:
# # virus_y1_typed.df[virus_y1_typed.df$LogNoS %in% dup_LogNo,]
# # virus_y2_typed.df[virus_y2_typed.df$LogNoG %in% dup_LogNo,]
# # virus_y3_typed.df[virus_y3_typed.df$LogNoG %in% dup_LogNo,]
# # ari_y1.df[ari_y1.df$LogNoS %in% dup_LogNo,]
# # ari_y1.df[ari_y1.df$LogNoG %in% dup_LogNo,]
# # ari_y1.df[ari_y1.df$LogNoG %in% dup_LogNo,]
# 
# # These are semi-automatically fixed as follows:
# useful_cols <- which(colnames(virus.all.merged.df) %in%
# 										 	c("age", 
# 										 		"Adeno", "AdVB", "AdVC", "AdVE", 
# 										 		"BoV", "Corona", "CoVOC43", "CoVNL63", "EV", 
# 										 		"HRV", "HRVA", "HRVB", "HRVC", 
# 										 		"Influ", "FluA", "FluB", "MPV", 
# 										 		"Paraflu", "PIV1", "PIV2", "PIV3", "PIV4b",
# 										 		"RSV", "RSVA", "RSVB",
# 										 		"sex",
# 										 		"uri", "lri", "flri", "wlri", "wflri"))
# 
# for (i in dup_LogNo) {
# 	# Identify row index of each duplicated LogNo
# 	x <- which(virus.all.merged.df$LogNo==i)
# 	
# 	# If more than 1 duplicate, break loop
# 	if(length(x) > 2) {
# 		paste("Error: need to resolve more than 1 duplicate for", i)
# 		break
# 	}
# 	
# 	# Subset out rows with that LogNo
# 	freq_var_test <- virus.all.merged.df[x,]
# 	
# 	# For each column/variable related to pathogen/infection/age:
# 	for (j in colnames(virus.all.merged.df)) {
# 		# if both rows have NAs, value is NA
# 		# else, if at least one NA or zero, summate and treat NAs as zeroes
# 		if (all(is.na(virus.all.merged.df[x,j]))) {
# 			freq_var_test[,j] <- NA
# 		} else if (any(is.na(virus.all.merged.df[x,j]) |
# 									 virus.all.merged.df[x,j]==0)) {
# 			freq_var_test[,j] <- sum(virus.all.merged.df[x,j], na.rm=T)
# 		}
# 		# else, keep original value
# 	}
# 	
# 	# If infection reported --- need to change Category to Infectious
# 	if (with(freq_var_test, sum(uri, lri, wlri, flri, na.rm=T)>0)) {
# 		freq_var_test$Category <- "Infectious"
# 	} else {
# 		freq_var_test$Category <- "control"
# 	}
# 	virus.all.merged.df[x,] <- freq_var_test
# }

virus.all.typed.df <- unique(virus.all.typed.df)
levels(virus.all.typed.df$Category) <- c("control", "infectious", "control")
virus.all.typed.df <- virus.all.typed.df %>%
	rename(SUBJ_ID = subjid,
				 LogNo_old = LogNo)

log.table.df2 <- 
	merge(log.table.df, virus.all.typed.df,
				by=c("SUBJ_ID", "LogNo_old", "age", "Category"),
				all=T)

table(log.table.df2$Category, useNA="ifany")

# We need to analyse HRVA,B,C subtypes separately from other cohorts
# In particular, for Year 1 samples, only look for subtypes in those samples with positive HRV
# Separate analyses for Year 1 and 2+3

df <- log.table.df2[which(log.table.df2$source == "Johnston"),]
with(df, table(anypcr, 
							 Adeno + Corona + HRV + MPV + Influ + Paraflu + RSV + Pico +
							 	Chlam + Myc != 0, 
							 useNA = "ifany"))
with(df, table(anypcr, 
							 Adeno + Corona + HRV + MPV + Influ + Paraflu + RSV + Pico != 0, 
							 useNA = "ifany"))
df[which(df$anypcr == 1 & (df$Adeno + df$Corona + df$HRV + df$MPV + df$Influ +
													 	df$Paraflu + df$RSV + df$Pico + 
													 	df$Chlam + df$Myc != 0)), 752:782]
df[which(df$anypcr == 1 & (df$Adeno + df$Corona + df$HRV + df$MPV + df$Influ + 
				 	df$Paraflu + df$RSV + df$Pico == 0)), 752:782]
df[which(df$anypcr == 1 & 
				 	is.na(df$Adeno + df$Corona + df$HRV + df$MPV + df$Influ +
				 					df$Paraflu + df$RSV + df$Pico)), 752:782]

# In Year 1 / Johnston analyses, anypcr represents any *viral* pcr, not Chlamydophila
# or Mycoplasma

log.table.df2$No.vir <- NA
log.table.df2$No.vir[which(log.table.df2$source == "Johnston")] <- 
	as.numeric(log.table.df2$anypcr[which(log.table.df2$source == "Johnston")] == 0)
			 
log.table.df2$No.vir[which(log.table.df2$source == "Gern")] <- 
	with(log.table.df2[which(log.table.df2$source == "Gern"),],
			 ifelse(
			 	ifelse(is.na(Adeno), 0, Adeno) +
			 		ifelse(is.na(BoV), 0, BoV) +
			 		ifelse(is.na(Corona), 0, Corona) +
			 		ifelse(is.na(EV), 0, EV) +
			 		ifelse(is.na(HRV), 0, HRV) +
			 		ifelse(is.na(Influ), 0, Influ) +
			 		ifelse(is.na(MPV), 0, MPV) +
			 		ifelse(is.na(Paraflu), 0, Paraflu) +
			 		ifelse(is.na(RSV), 0, RSV) == 0, 1, 0))

# Get top viruses
# ---------------------------------------
vir.count.table.y1.df <- with(log.table.df2[which(log.table.df2$source == "Johnston"),], 
													 data.frame(LogNo,
													 					 HRV, HRVA, HRVB, HRVC, RSV, Influ, Adeno,
													 					 Corona, MPV, Paraflu, Chlam, Myc, Pico,
													 					 No.vir))
vir.count.table.y1.df.m <- melt(vir.count.table.y1.df, id.var = c("LogNo"),
														 variable.name="virus.name",
														 value.name="count")
vir.count.table.y1.df.m2 <- with(vir.count.table.y1.df.m, 
															as.data.frame.matrix(table(virus.name, count == 1)))
colnames(vir.count.table.y1.df.m2) <- c("Absent", "Present")

vir.count.table.y1.df.m2$Total <- with(vir.count.table.y1.df.m2,
																		Absent + Present)
write.csv(vir.count.table.y1.df.m2, "output/2018-02-12_virus_count_y1.csv", quote=F,
					row.names=T)


vir.count.table.y23.df <- with(log.table.df2[which(log.table.df2$source == "Gern"),], 
													 data.frame(LogNo,
													 					 HRV, HRVA, HRVB, HRVC, RSV, Influ, Adeno,
													 					 Corona, MPV, Paraflu, EV, BoV,
													 					 No.vir))
vir.count.table.y23.df.m <- melt(vir.count.table.y23.df, id.var = c("LogNo"),
														 variable.name="virus.name",
														 value.name="count")
vir.count.table.y23.df.m2 <- with(vir.count.table.y23.df.m, 
															as.data.frame.matrix(table(virus.name, count == 1)))
colnames(vir.count.table.y23.df.m2) <- c("Absent", "Present")

vir.count.table.y23.df.m2$Total <- with(vir.count.table.y23.df.m2,
																		Absent + Present)
write.csv(vir.count.table.y23.df.m2, "output/2018-02-12_virus_count_y23.csv", quote=F,
					row.names=T)

# maximum age for johnston samples
with(log.table.df2[which(log.table.df2$source == "Johnston"),], max(age, na.rm=T))
with(log.table.df2[which(log.table.df2$source == "Johnston" &
												 	log.table.df2$age != 392),], max(age, na.rm=T))

# Minimum age for gern samples
with(log.table.df2[which(log.table.df2$source == "Gern"),], min(age, na.rm=T))

# Top viruses:
# RV.new, RSV, Parainfluenza, Coronavirus, Adenovirus, Bocavirus, Influenza, MPV, Other

virus.real.name.y1 <- c("Adenovirus", "RV", "Picornavirus",
										 "Coronavirus", "Influenza", "RSV", "Parainfluenza", 
										  "MPV", "None")
virus.var.name.y1 <- c("Adeno", "HRV",  "Pico", 
										"Corona", "Influ",  "RSV", "Paraflu",
										"MPV", "No.vir")
virus.simp.name.y1 <- c("adeno", "RV",  "pico",
										 "corona", "influ", "RSV", "para",
										  "MPV", "novir")
virus.colors.y1 <- c("#33a02c", "#ff0000",  "#cd5b45", 
										 "#ff7f00", "#1f78b4", "#a6cee3", "#00008B", 
										 "#6a3d9a", "#808080")

# "RV-A", "RV-B", "RV-C", 
# "HRVA", "HRVB", "HRVC",
# "RVA", "RVB", "RVC",
# "#800000", "#fb9a99", "#e31a1c",

virus.real.name.y23 <- c("Adenovirus", "Bocavirus", "RV", "Enterovirus",
												"Coronavirus", "Influenza", "RSV", "Parainfluenza", 
												"MPV", "None")
virus.var.name.y23 <- c("Adeno", "BoV", "HRV",   "EV",
											 "Corona", "Influ",  "RSV", "Paraflu",
											 "MPV", "No.vir")
virus.simp.name.y23 <- c("adeno", "boca", "RV", "EV",
												"corona", "influ", "RSV", "para",
												"MPV", "novir")
virus.colors.y23 <- c("#33a02c", "#b2df8a", "#ff0000", "#cd5b45", 
										 "#ff7f00", "#1f78b4", "#a6cee3", "#00008B", 
										 "#6a3d9a", "#808080")

# "RV-A", "RV-B", "RV-C", 
# "HRVA", "HRVB", "HRVC",
# "RVA", "RVB", "RVC", 
# "#800000", "#fb9a99", "#e31a1c",

# Graph vs. time (Y1 samples)
# ---------------------------------------

log.table.y1.df2 <- log.table.df2[which(log.table.df2$source == "Johnston"),]
vir.age.distr.y1.df <- log.table.y1.df2[,which(colnames(log.table.y1.df2) %in% 
																					c("LogNo_old", "age", "well.sick", virus.var.name.y1))] %>%
	mutate(age.custom = cut(age, 
													breaks=c(0, 365/12*4, 365/12*9, 365*1.25),
													labels=c("2mth", "6mth", "1y")))

vir.age.distr.y1.df.m <- melt(vir.age.distr.y1.df, 
													 id.var = c("LogNo_old", "age", "age.custom", "well.sick"),
													 variable.name="virus.name",
													 value.name="Count") %>%
	mutate(Present = 1) %>%
	group_by(age.custom, virus.name, well.sick) %>%
	summarize(sum = sum(Count, na.rm=T), 
						sum.all = sum(Present, na.rm=T)) %>%
	mutate(prop = sum/sum.all)
vir.age.distr.y1.df.m$virus.name <- factor(vir.age.distr.y1.df.m$virus.name,
													levels=virus.var.name.y1,
													labels=virus.real.name.y1)

pdf("output/2018-01-31_prop_mild_sick_virus_y1.pdf", width=8, height=5)
vir.age.distr.y1.df.m %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=prop, group=virus.name, color=virus.name)) + geom_point() + 
	geom_line() + 
	scale_color_manual(values = virus.colors.y1, 
										 name = "Virus") +
	scale_x_discrete(name="Timepoint of collection") + 
	scale_y_continuous(name="Proportion of samples") + 
	theme(panel.grid.minor = element_blank()) + 
	facet_grid(well.sick~.)
dev.off()

pdf("output/2018-01-31_count_mild_sick_virus_y1.pdf", width=8, height=5)
vir.age.distr.y1.df.m %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=sum, group=virus.name, color=virus.name)) + geom_point() + 
	geom_line() + 
	scale_color_manual(values = virus.colors.y1, 
										 name = "Virus") +										 
	scale_x_discrete(name="Timepoint of collection") + 
	scale_y_continuous(name="Number of samples") + 
	theme(panel.grid.minor = element_blank()) + 
	facet_grid(well.sick~.)
dev.off()

# Graph vs. time (Y23 samples)
# ---------------------------------------

log.table.y23.df2 <- log.table.df2[which(log.table.df2$source == "Gern"),]
vir.age.distr.y23.df <- log.table.y23.df2[,which(colnames(log.table.y23.df2) %in% 
																c("LogNo_old", "age", "well.sick", virus.var.name.y23))] %>%
	mutate(age.custom = cut(age, 
													breaks=c(365/12*9, 365*1.25, 365*1.75, 365*2.25, 365*2.75, 365*3.25),
													labels=c("1y", "1.5y", "2y", "2.5y", "3y")))

vir.age.distr.y23.df.m <- melt(vir.age.distr.y23.df, 
													 id.var = c("LogNo_old", "age", "age.custom", "well.sick"),
													 variable.name="virus.name",
													 value.name="Count") %>%
	mutate(Present = 1) %>%
	group_by(age.custom, virus.name, well.sick) %>%
	summarize(sum = sum(Count, na.rm=T), 
						sum.all = sum(Present, na.rm=T)) %>%
	mutate(prop = sum/sum.all)
vir.age.distr.y23.df.m$virus.name <- factor(vir.age.distr.y23.df.m$virus.name,
																				levels=virus.var.name.y23,
																				labels=virus.real.name.y23)

pdf("output/2018-01-31_prop_mild_sick_virus_y23.pdf", width=8, height=5)
vir.age.distr.y23.df.m %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=prop, group=virus.name, color=virus.name)) + geom_point() + 
	geom_line() + 
	scale_color_manual(values = virus.colors.y23, 
										 name = "Virus") +
	scale_x_discrete(name="Timepoint of collection") + 
	scale_y_continuous(name="Proportion of samples") + 
	theme(panel.grid.minor = element_blank()) + 
	facet_grid(well.sick~.)
dev.off()

pdf("output/2018-01-31_count_mild_sick_virus_y23.pdf", width=8, height=5)
vir.age.distr.y23.df.m %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=sum, group=virus.name, color=virus.name)) + geom_point() + 
	geom_line() + 
	scale_color_manual(values = virus.colors.y23, 
										 name = "Virus") +										 
	scale_x_discrete(name="Timepoint of collection") + 
	scale_y_continuous(name="Number of samples") + 
	theme(panel.grid.minor = element_blank()) + 
	facet_grid(well.sick~.)
dev.off()

# Virus distribution table (y1 samples)
# ---------------------------------------
for (i in 1:(length(virus.var.name.y1))) {
	df <- cbind(
		data.frame(with(log.table.y1.df2, table(get(virus.var.name.y1[i])))), 
		as.data.frame.matrix(with(log.table.y1.df2, 
															table(get(virus.var.name.y1[i]), well.sick)))
	)
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	colnames(df)[1:2] <- c("Virus", "All")
	assign(paste0("vir.count.", virus.simp.name.y1[i], ".df"), df)
	
	df <- cbind(
		data.frame(with(log.table.y1.df2, prop.table(table(get(virus.var.name.y1[i]))))),
		as.data.frame.matrix(with(log.table.y1.df2, 
															prop.table(table(get(virus.var.name.y1[i]), well.sick),2)))
	)
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	colnames(df)[1:2] <- c("Virus", "All")
	assign(paste0("vir.prop.", virus.simp.name.y1[i], ".df"), df)	
}

vir.count.y1.df <- do.call("rbind", mget(paste0("vir.count.", 
																									 virus.simp.name.y1, ".df")))
vir.count.y1.df$Virus <- NULL
row.names(vir.count.y1.df) <- NULL
colnames(vir.count.y1.df) <- c("All.Freq", "Well.Freq", "Sick.Freq", "Status")

vir.prop.y1.df <- do.call("rbind", mget(paste0("vir.prop.", 
																						 virus.simp.name.y1, ".df")))
vir.prop.y1.df$Virus <- NULL
row.names(vir.prop.y1.df) <- NULL
colnames(vir.prop.y1.df) <- c("All.Prop", "Well.Prop", "Sick.Prop", "Status")

vir.distr.y1.df <- data.frame(vir.name = vir.count.y1.df$Status,
													 All.Freq = vir.count.y1.df$All.Freq,
													 All.Prop = paste0(good.signif(vir.prop.y1.df$All.Prop*100,2), "%"),
													 Well.Freq = vir.count.y1.df$Well.Freq,
													 Well.Prop = paste0(good.signif(vir.prop.y1.df$Well.Prop*100,2), "%"),
													 Sick.Freq = vir.count.y1.df$Sick.Freq,
													 Sick.Prop = paste0(good.signif(vir.prop.y1.df$Sick.Prop*100,2), "%"), 
													 stringsAsFactors = F)
vir.distr.y1.df2 <- data.frame(vir.name = vir.distr.y1.df$vir.name,
														All = paste0(vir.distr.y1.df$All.Freq, " (", 
																				 vir.distr.y1.df$All.Prop, ")"),
														Well = paste0(vir.distr.y1.df$Well.Freq, " (", 
																					vir.distr.y1.df$Well.Prop, ")"),
														Sick = paste0(vir.distr.y1.df$Sick.Freq, " (", 
																					vir.distr.y1.df$Sick.Prop, ")"))

vir.distr.y1.df3 <- vir.distr.y1.df2[seq(2, nrow(vir.distr.y1.df2), 2),]

# Virus distribution table (y23 samples)
# ---------------------------------------
for (i in 1:(length(virus.var.name.y23))) {
	df <- cbind(
		data.frame(with(log.table.y23.df2, table(get(virus.var.name.y23[i])))), 
		as.data.frame.matrix(with(log.table.y23.df2, 
															table(get(virus.var.name.y23[i]), well.sick)))
	)
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	colnames(df)[1:2] <- c("Virus", "All")
	assign(paste0("vir.count.", virus.simp.name.y23[i], ".df"), df)
	
	df <- cbind(
		data.frame(with(log.table.y23.df2, prop.table(table(get(virus.var.name.y23[i]))))),
		as.data.frame.matrix(with(log.table.y23.df2, 
															prop.table(table(get(virus.var.name.y23[i]), well.sick),2)))
	)
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	colnames(df)[1:2] <- c("Virus", "All")
	assign(paste0("vir.prop.", virus.simp.name.y23[i], ".df"), df)	
}

vir.count.y23.df <- do.call("rbind", mget(paste0("vir.count.", 
																								virus.simp.name.y23, ".df")))
vir.count.y23.df$Virus <- NULL
row.names(vir.count.y23.df) <- NULL
colnames(vir.count.y23.df) <- c("All.Freq", "Well.Freq", "Sick.Freq", "Status")

vir.prop.y23.df <- do.call("rbind", mget(paste0("vir.prop.", 
																							 virus.simp.name.y23, ".df")))
vir.prop.y23.df$Virus <- NULL
row.names(vir.prop.y23.df) <- NULL
colnames(vir.prop.y23.df) <- c("All.Prop", "Well.Prop", "Sick.Prop", "Status")

vir.distr.y23.df <- data.frame(vir.name = vir.count.y23.df$Status,
															All.Freq = vir.count.y23.df$All.Freq,
															All.Prop = paste0(good.signif(vir.prop.y23.df$All.Prop*100,2), "%"),
															Well.Freq = vir.count.y23.df$Well.Freq,
															Well.Prop = paste0(good.signif(vir.prop.y23.df$Well.Prop*100,2), "%"),
															Sick.Freq = vir.count.y23.df$Sick.Freq,
															Sick.Prop = paste0(good.signif(vir.prop.y23.df$Sick.Prop*100,2), "%"), 
															stringsAsFactors = F)
vir.distr.y23.df2 <- data.frame(vir.name = vir.distr.y23.df$vir.name,
															 All = paste0(vir.distr.y23.df$All.Freq, " (", 
															 						 vir.distr.y23.df$All.Prop, ")"),
															 Well = paste0(vir.distr.y23.df$Well.Freq, " (", 
															 							vir.distr.y23.df$Well.Prop, ")"),
															 Sick = paste0(vir.distr.y23.df$Sick.Freq, " (", 
															 							vir.distr.y23.df$Sick.Prop, ")"))

vir.distr.y23.df3 <- vir.distr.y23.df2[seq(2, nrow(vir.distr.y23.df2), 2),]

# Graph vs. season (y1 sample)
# ---------------------------------------

for (i in 1:length(virus.var.name.y1)) {
	df <- with(log.table.y1.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y1.df2[which(log.table.y1.df2$well.sick == "Well"),],
						 			 table(season, get(virus.var.name.y1[i]))
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y1[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.well.", virus.simp.name.y1[i], ".df"), df)
	
	df <- with(log.table.y1.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y1.df2[which(log.table.y1.df2$well.sick == "Sick"),],
						 			 table(season, get(virus.var.name.y1[i]))
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y1[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.sick.", virus.simp.name.y1[i], ".df"), df)
}

# Fix well.influ and well.MPV
vir.season.well.influ.df$"1" <- 0
vir.season.well.influ.df <- vir.season.well.influ.df[,c(1,4,2,3)]
vir.season.well.MPV.df$"1" <- 0
vir.season.well.MPV.df <- vir.season.well.MPV.df[,c(1,4,2,3)]

vir.season.well.y1.df <- do.call("rbind", mget(paste0("vir.season.well.", 
																									 virus.simp.name.y1, ".df")))
row.names(vir.season.well.y1.df) <- NULL
vir.season.well.y1.df$well.sick <- "Well"
vir.season.well.y1.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 9)
colnames(vir.season.well.y1.df)[2] <- "Count"

vir.season.sick.y1.df <- do.call("rbind", mget(paste0("vir.season.sick.", 
																											virus.simp.name.y1, ".df")))
row.names(vir.season.sick.y1.df) <- NULL
vir.season.sick.y1.df$well.sick <- "Sick"
vir.season.sick.y1.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 9)
colnames(vir.season.sick.y1.df)[2] <- "Count"

vir.season.well.sick.y1.df <- rbind(vir.season.well.y1.df, vir.season.sick.y1.df)
vir.season.well.sick.y1.df$well.sick <- 
	factor(vir.season.well.sick.y1.df$well.sick,
				 levels=c("Well", "Sick"))
vir.season.well.sick.y1.df$vir.name <- 
	factor(vir.season.well.sick.y1.df$vir.name,
				 levels = virus.real.name.y1)
vir.season.well.sick.y1.df$season <- 
	factor(vir.season.well.sick.y1.df$season,
				 levels = c("Winter", "Spring", "Summer", "Autumn"),
				 labels = c("Winter", "Spring", "Summer", "Fall"))

pdf("output/2018-02-12_vir_distr_season_y1.count.pdf", width=8, height=5)
vir.season.well.sick.y1.df %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=season, y=Count, group=vir.name)) + 
	geom_line(aes(color=vir.name)) + 
	geom_point(aes(color=vir.name)) + 
	scale_color_manual(values = virus.colors.y1) +
	facet_grid(.~well.sick) +
	labs(color = "Virus",
			 x = "Season of collection",
			 y = "Number of samples")
dev.off()

for (i in 1:length(virus.var.name.y1)) {
	df <- with(log.table.y1.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y1.df2[which(log.table.y1.df2$well.sick == "Well"),],
						 			 prop.table(table(season, get(virus.var.name.y1[i])),1)
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y1[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.well.", virus.simp.name.y1[i], ".df"), df)
	
	df <- with(log.table.y1.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y1.df2[which(log.table.y1.df2$well.sick == "Sick"),],
						 			 prop.table(table(season, get(virus.var.name.y1[i])),1)
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y1[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.sick.", virus.simp.name.y1[i], ".df"), df)
}

# Fix well.influ and well.MPV
vir.season.well.influ.df$"1" <- 0
vir.season.well.influ.df <- vir.season.well.influ.df[,c(1,4,2,3)]
vir.season.well.MPV.df$"1" <- 0
vir.season.well.MPV.df <- vir.season.well.MPV.df[,c(1,4,2,3)]

vir.season.well.y1.df <- do.call("rbind", mget(paste0("vir.season.well.", 
																									 virus.simp.name.y1, ".df")))
row.names(vir.season.well.y1.df) <- NULL
vir.season.well.y1.df$well.sick <- "Well"
vir.season.well.y1.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 9)
colnames(vir.season.well.y1.df)[2] <- "Prop"

vir.season.sick.y1.df <- do.call("rbind", mget(paste0("vir.season.sick.", 
																											virus.simp.name.y1, ".df")))
row.names(vir.season.sick.y1.df) <- NULL
vir.season.sick.y1.df$well.sick <- "Sick"
vir.season.sick.y1.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 9)
colnames(vir.season.sick.y1.df)[2] <- "Prop"

vir.season.well.sick.y1.df <- rbind(vir.season.well.y1.df, vir.season.sick.y1.df)
vir.season.well.sick.y1.df$well.sick <- 
	factor(vir.season.well.sick.y1.df$well.sick, 
				 levels=c("Well", "Sick"))
vir.season.well.sick.y1.df$vir.name <- 
	factor(vir.season.well.sick.y1.df$vir.name,
				 levels = virus.real.name.y1)
vir.season.well.sick.y1.df$season <- 
	factor(vir.season.well.sick.y1.df$season,
				 levels = c("Winter", "Spring", "Summer", "Autumn"),
				 labels = c("Winter", "Spring", "Summer", "Fall"))

pdf("output/2018-02-12_vir_distr_season_y1.pdf", width=8, height=5)
vir.season.well.sick.y1.df %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=season, y=Prop, group=vir.name)) + 
	geom_line(aes(color=vir.name)) + 
	geom_point(aes(color=vir.name)) + 
	scale_color_manual(values = virus.colors.y1) +
	facet_grid(.~well.sick) +
	labs(color = "Virus",
			 x = "Season of collection",
			 y = "Proportion")
dev.off()


# Graph vs. season (y23 sample)
# ---------------------------------------

for (i in 1:length(virus.var.name.y23)) {
	df <- with(log.table.y23.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y23.df2[which(log.table.y23.df2$well.sick == "Well"),],
						 			 table(season, get(virus.var.name.y23[i]))
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y23[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.well.", virus.simp.name.y23[i], ".df"), df)
	
	df <- with(log.table.y23.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y23.df2[which(log.table.y23.df2$well.sick == "Sick"),],
						 			 table(season, get(virus.var.name.y23[i]))
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y23[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.sick.", virus.simp.name.y23[i], ".df"), df)
}

vir.season.well.y23.df <- do.call("rbind", mget(paste0("vir.season.well.", 
																											virus.simp.name.y23, ".df")))
row.names(vir.season.well.y23.df) <- NULL
vir.season.well.y23.df$well.sick <- "Well"
vir.season.well.y23.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 10)
colnames(vir.season.well.y23.df)[2] <- "Count"

vir.season.sick.y23.df <- do.call("rbind", mget(paste0("vir.season.sick.", 
																											virus.simp.name.y23, ".df")))
row.names(vir.season.sick.y23.df) <- NULL
vir.season.sick.y23.df$well.sick <- "Sick"
vir.season.sick.y23.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 10)
colnames(vir.season.sick.y23.df)[2] <- "Count"

vir.season.well.sick.y23.df <- rbind(vir.season.well.y23.df, vir.season.sick.y23.df)
vir.season.well.sick.y23.df$well.sick <- 
	factor(vir.season.well.sick.y23.df$well.sick,
				 levels=c("Well", "Sick"))
vir.season.well.sick.y23.df$vir.name <- 
	factor(vir.season.well.sick.y23.df$vir.name,
				 levels = virus.real.name.y23)
vir.season.well.sick.y23.df$season <- 
	factor(vir.season.well.sick.y23.df$season,
				 levels = c("Winter", "Spring", "Summer", "Autumn"),
				 labels = c("Winter", "Spring", "Summer", "Fall"))

pdf("output/2018-02-12_vir_distr_season_y23.count.pdf", width=8, height=5)
vir.season.well.sick.y23.df %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=season, y=Count, group=vir.name)) + 
	geom_line(aes(color=vir.name)) + 
	geom_point(aes(color=vir.name)) + 
	scale_color_manual(values = virus.colors.y23) +
	facet_grid(.~well.sick) +
	labs(color = "Virus",
			 x = "Season of collection",
			 y = "Number of samples")
dev.off()

for (i in 1:length(virus.var.name.y23)) {
	df <- with(log.table.y23.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y23.df2[which(log.table.y23.df2$well.sick == "Well"),],
						 			 prop.table(table(season, get(virus.var.name.y23[i])),1)
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y23[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.well.", virus.simp.name.y23[i], ".df"), df)
	
	df <- with(log.table.y23.df2, 
						 as.data.frame.matrix(
						 	with(log.table.y23.df2[which(log.table.y23.df2$well.sick == "Sick"),],
						 			 prop.table(table(season, get(virus.var.name.y23[i])),1)
						 	)
						 )
	)
	df$vir.name <- virus.real.name.y23[i]
	df$season <- row.names(df)
	row.names(df) <- NULL
	assign(paste0("vir.season.sick.", virus.simp.name.y23[i], ".df"), df)
}

vir.season.well.y23.df <- do.call("rbind", mget(paste0("vir.season.well.", 
																											virus.simp.name.y23, ".df")))
row.names(vir.season.well.y23.df) <- NULL
vir.season.well.y23.df$well.sick <- "Well"
vir.season.well.y23.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 10)
colnames(vir.season.well.y23.df)[2] <- "Prop"

vir.season.sick.y23.df <- do.call("rbind", mget(paste0("vir.season.sick.", 
																											virus.simp.name.y23, ".df")))
row.names(vir.season.sick.y23.df) <- NULL
vir.season.sick.y23.df$well.sick <- "Sick"
vir.season.sick.y23.df$season <- rep(c("Winter", "Spring", "Summer", "Autumn"), 10)
colnames(vir.season.sick.y23.df)[2] <- "Prop"

vir.season.well.sick.y23.df <- rbind(vir.season.well.y23.df, vir.season.sick.y23.df)
vir.season.well.sick.y23.df$well.sick <- 
	factor(vir.season.well.sick.y23.df$well.sick, 
				 levels=c("Well", "Sick"))
vir.season.well.sick.y23.df$vir.name <- 
	factor(vir.season.well.sick.y23.df$vir.name,
				 levels = virus.real.name.y23)
vir.season.well.sick.y23.df$season <- 
	factor(vir.season.well.sick.y23.df$season,
				 levels = c("Winter", "Spring", "Summer", "Autumn"),
				 labels = c("Winter", "Spring", "Summer", "Fall"))

pdf("output/2018-02-12_vir_distr_season_y23.pdf", width=8, height=5)
vir.season.well.sick.y23.df %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=season, y=Prop, group=vir.name)) + 
	geom_line(aes(color=vir.name)) + 
	geom_point(aes(color=vir.name)) + 
	scale_color_manual(values = virus.colors.y23) +
	facet_grid(.~well.sick) +
	labs(color = "Virus",
			 x = "Season of collection",
			 y = "Proportion")
dev.off()

# Well/sick ~ virus (+ gender + age + season) - y1
# ---------------------------------------
vir.list <- virus.var.name.y1

results.df <- data.frame(OR = character(length(vir.list)),
												 p = character(length(vir.list)),
												 OR_alt = character(length(vir.list)),
												 p_alt = character(length(vir.list)),
												 stringsAsFactors=F)

for(i in 1:length(vir.list)){
	gee.out <- tryCatch(gee(as.formula(paste0('I(well.sick == "Sick") ~ ', 
																	 vir.list[i], 
																	 '+ sex + age')),
								 id = SUBJ_ID, 
								 data = log.table.y1.df2, 
								 family= binomial,
								 corstr = "unstructured"),
								 error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=4, ncol=5))
	} else{
		s  <- summary(gee.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}
vir.distr.y1.df3$OR <- results.df$OR
vir.distr.y1.df3$p <- results.df$p

for(i in 1:length(vir.list)){
	gee.out <- tryCatch(gee(as.formula(paste0('I(well.sick == "Sick") ~ ', 
																	 vir.list[i], 
																	 '+ sex + age + I(season2=="Winter.Spring")')),
								 id = SUBJ_ID, 
								 data = log.table.y1.df2, 
								 family= binomial,
								 corstr = "unstructured"),
								 error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=4, ncol=5))
	} else{
		s  <- summary(gee.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR_alt[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p_alt[i] <- paste0(p)[2]
}
vir.distr.y1.df3$OR_alt <- results.df$OR_alt
vir.distr.y1.df3$p_alt <- results.df$p_alt
write.csv(vir.distr.y1.df3, "output/2018-02-13_virus.distr_y1.csv", 
					quote=F, row.names=F)


# Well/sick ~ virus (+ gender + age + season) - y23
# ---------------------------------------
vir.list <- virus.var.name.y23

results.df <- data.frame(OR = character(length(vir.list)),
												 p = character(length(vir.list)),
												 OR_alt = character(length(vir.list)),
												 p_alt = character(length(vir.list)),
												 stringsAsFactors=F)

for(i in 1:length(vir.list)){
	gee.out <- tryCatch(gee(as.formula(paste0('I(well.sick == "Sick") ~ ', 
																						vir.list[i], 
																						'+ sex + age')),
													id = SUBJ_ID, 
													data = log.table.y23.df2, 
													family= binomial,
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=4, ncol=5))
	} else{
		s  <- summary(gee.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}
vir.distr.y23.df3$OR <- results.df$OR
vir.distr.y23.df3$p <- results.df$p

for(i in 1:length(vir.list)){
	gee.out <- tryCatch(gee(as.formula(paste0('I(well.sick == "Sick") ~ ', 
																						vir.list[i], 
																						'+ sex + age + I(season2=="Winter.Spring")')),
													id = SUBJ_ID, 
													data = log.table.y23.df2, 
													family= binomial,
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=4, ncol=5))
	} else{
		s  <- summary(gee.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR_alt[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p_alt[i] <- paste0(p)[2]
}
vir.distr.y23.df3$OR_alt <- results.df$OR_alt
vir.distr.y23.df3$p_alt <- results.df$p_alt
write.csv(vir.distr.y23.df3, "output/2018-02-13_virus.distr_y23.csv", 
					quote=F, row.names=F)

# Virus distribution table ~ season, well/sick - y1
# ---------------------------------------

for (i in 1:length(virus.var.name.y1)) {
	df <- as.data.frame.matrix(with(log.table.y1.df2, 
															table(get(virus.var.name.y1[i]), season)))
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.season.", virus.simp.name.y1[i], ".df"), df)
	
	df <- as.data.frame.matrix(with(log.table.y1.df2, 
															prop.table(table(get(virus.var.name.y1[i]), season),2)))
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.season.", virus.simp.name.y1[i], ".df"), df)	
}

vir.count.season.y1.df <- do.call("rbind", mget(paste0("vir.count.season.", 
																						 virus.simp.name.y1, ".df")))
row.names(vir.count.season.y1.df) <- NULL
colnames(vir.count.season.y1.df) <- c("Winter.Freq", "Spring.Freq", "Summer.Freq", 
																	 "Fall.Freq", "Status")

vir.prop.season.y1.df <- do.call("rbind", mget(paste0("vir.prop.season.", 
																											virus.simp.name.y1, ".df")))
row.names(vir.prop.season.y1.df) <- NULL
colnames(vir.prop.season.y1.df) <- c("Winter.Prop", "Spring.Prop", "Summer.Prop", 
																	"Fall.Prop", "Status")

vir.season.distr.y1.df <- data.frame(vir.name = vir.count.season.y1.df$Status,
													 Winter.Freq = vir.count.season.y1.df$Winter.Freq,
													 Winter.Prop = paste0(good.signif(vir.prop.season.y1.df$Winter.Prop*100,2), "%"),
													 Spring.Freq = vir.count.season.y1.df$Spring.Freq,
													 Spring.Prop = paste0(good.signif(vir.prop.season.y1.df$Spring.Prop*100,2), "%"),
													 Summer.Freq = vir.count.season.y1.df$Summer.Freq,
													 Summer.Prop = paste0(good.signif(vir.prop.season.y1.df$Summer.Prop*100,2), "%"),
													 Fall.Freq = vir.count.season.y1.df$Fall.Freq,
													 Fall.Prop = paste0(good.signif(vir.prop.season.y1.df$Fall.Prop*100,2), "%"), 
													 stringsAsFactors = F)
vir.season.distr.y1.df2 <- data.frame(vir.name = vir.season.distr.y1.df$vir.name,
														Winter = paste0(vir.season.distr.y1.df$Winter.Freq, " (", 
																						vir.season.distr.y1.df$Winter.Prop, ")"),
														Spring = paste0(vir.season.distr.y1.df$Spring.Freq, " (", 
																						vir.season.distr.y1.df$Spring.Prop, ")"),
														Summer = paste0(vir.season.distr.y1.df$Summer.Freq, " (", 
																						vir.season.distr.y1.df$Summer.Prop, ")"),
														Fall = paste0(vir.season.distr.y1.df$Fall.Freq, " (", 
																					vir.season.distr.y1.df$Fall.Prop, ")"))

vir.season.distr.y1.df3 <- vir.season.distr.y1.df2[seq(2, nrow(vir.season.distr.y1.df2), 2),]


# Virus distribution table ~ season, well/sick - y23
# ---------------------------------------

for (i in 1:length(virus.var.name.y23)) {
	df <- as.data.frame.matrix(with(log.table.y23.df2, 
																	table(get(virus.var.name.y23[i]), season)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.season.", virus.simp.name.y23[i], ".df"), df)
	
	df <- as.data.frame.matrix(with(log.table.y23.df2, 
																	prop.table(table(get(virus.var.name.y23[i]), season),2)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.season.", virus.simp.name.y23[i], ".df"), df)	
}

vir.count.season.y23.df <- do.call("rbind", mget(paste0("vir.count.season.", 
																											 virus.simp.name.y23, ".df")))
row.names(vir.count.season.y23.df) <- NULL
colnames(vir.count.season.y23.df) <- c("Winter.Freq", "Spring.Freq", "Summer.Freq", 
																			"Fall.Freq", "Status")

vir.prop.season.y23.df <- do.call("rbind", mget(paste0("vir.prop.season.", 
																											virus.simp.name.y23, ".df")))
row.names(vir.prop.season.y23.df) <- NULL
colnames(vir.prop.season.y23.df) <- c("Winter.Prop", "Spring.Prop", "Summer.Prop", 
																		 "Fall.Prop", "Status")

vir.season.distr.y23.df <- data.frame(vir.name = vir.count.season.y23.df$Status,
																		 Winter.Freq = vir.count.season.y23.df$Winter.Freq,
																		 Winter.Prop = paste0(good.signif(vir.prop.season.y23.df$Winter.Prop*100,2), "%"),
																		 Spring.Freq = vir.count.season.y23.df$Spring.Freq,
																		 Spring.Prop = paste0(good.signif(vir.prop.season.y23.df$Spring.Prop*100,2), "%"),
																		 Summer.Freq = vir.count.season.y23.df$Summer.Freq,
																		 Summer.Prop = paste0(good.signif(vir.prop.season.y23.df$Summer.Prop*100,2), "%"),
																		 Fall.Freq = vir.count.season.y23.df$Fall.Freq,
																		 Fall.Prop = paste0(good.signif(vir.prop.season.y23.df$Fall.Prop*100,2), "%"), 
																		 stringsAsFactors = F)
vir.season.distr.y23.df2 <- data.frame(vir.name = vir.season.distr.y23.df$vir.name,
																			Winter = paste0(vir.season.distr.y23.df$Winter.Freq, " (", 
																											vir.season.distr.y23.df$Winter.Prop, ")"),
																			Spring = paste0(vir.season.distr.y23.df$Spring.Freq, " (", 
																											vir.season.distr.y23.df$Spring.Prop, ")"),
																			Summer = paste0(vir.season.distr.y23.df$Summer.Freq, " (", 
																											vir.season.distr.y23.df$Summer.Prop, ")"),
																			Fall = paste0(vir.season.distr.y23.df$Fall.Freq, " (", 
																										vir.season.distr.y23.df$Fall.Prop, ")"))

vir.season.distr.y23.df3 <- vir.season.distr.y23.df2[seq(2, nrow(vir.season.distr.y23.df2), 2),]

# Season ~ virus (+ gender + age + well/sick) - y1
# ---------------------------------------
vir.list <- virus.var.name.y1

results.df <- data.frame(OR = character(length(vir.list)),
												 p = character(length(vir.list)),
												 OR_alt = character(length(vir.list)),
												 p_alt = character(length(vir.list)),
												 stringsAsFactors=F)
for(i in 1:length(vir.list)){
	gee.out <- tryCatch(gee(as.formula(paste0('I(season2=="Winter.Spring") ~ ', 
																	 vir.list[i], 
																	 '+ sex + age + I(well.sick == "Sick")')),
								 id = SUBJ_ID, 
								 data = log.table.y1.df2, 
								 family= binomial,
								 corstr = "unstructured"),
								 error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=4, ncol=5))
	} else{
		s  <- summary(gee.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}
vir.season.distr.y1.df3$OR <- results.df$OR
vir.season.distr.y1.df3$p <- results.df$p

write.csv(vir.season.distr.y1.df3, "output/2018-02-12_vir.season.distr.y1.csv", quote=F,
					row.names=F)

# Virus distribution table ~ season, well/sick - y23
# ---------------------------------------

for (i in 1:length(virus.var.name.y23)) {
	df <- as.data.frame.matrix(with(log.table.y23.df2, 
																	table(get(virus.var.name.y23[i]), season)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.season.", virus.simp.name.y23[i], ".df"), df)
	
	df <- as.data.frame.matrix(with(log.table.y23.df2, 
																	prop.table(table(get(virus.var.name.y23[i]), season),2)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.season.", virus.simp.name.y23[i], ".df"), df)	
}

vir.count.season.y23.df <- do.call("rbind", mget(paste0("vir.count.season.", 
																												virus.simp.name.y23, ".df")))
row.names(vir.count.season.y23.df) <- NULL
colnames(vir.count.season.y23.df) <- c("Winter.Freq", "Spring.Freq", "Summer.Freq", 
																			 "Fall.Freq", "Status")

vir.prop.season.y23.df <- do.call("rbind", mget(paste0("vir.prop.season.", 
																											 virus.simp.name.y23, ".df")))
row.names(vir.prop.season.y23.df) <- NULL
colnames(vir.prop.season.y23.df) <- c("Winter.Prop", "Spring.Prop", "Summer.Prop", 
																			"Fall.Prop", "Status")

vir.season.distr.y23.df <- data.frame(vir.name = vir.count.season.y23.df$Status,
																			Winter.Freq = vir.count.season.y23.df$Winter.Freq,
																			Winter.Prop = paste0(good.signif(vir.prop.season.y23.df$Winter.Prop*100,2), "%"),
																			Spring.Freq = vir.count.season.y23.df$Spring.Freq,
																			Spring.Prop = paste0(good.signif(vir.prop.season.y23.df$Spring.Prop*100,2), "%"),
																			Summer.Freq = vir.count.season.y23.df$Summer.Freq,
																			Summer.Prop = paste0(good.signif(vir.prop.season.y23.df$Summer.Prop*100,2), "%"),
																			Fall.Freq = vir.count.season.y23.df$Fall.Freq,
																			Fall.Prop = paste0(good.signif(vir.prop.season.y23.df$Fall.Prop*100,2), "%"), 
																			stringsAsFactors = F)
vir.season.distr.y23.df2 <- data.frame(vir.name = vir.season.distr.y23.df$vir.name,
																			 Winter = paste0(vir.season.distr.y23.df$Winter.Freq, " (", 
																			 								vir.season.distr.y23.df$Winter.Prop, ")"),
																			 Spring = paste0(vir.season.distr.y23.df$Spring.Freq, " (", 
																			 								vir.season.distr.y23.df$Spring.Prop, ")"),
																			 Summer = paste0(vir.season.distr.y23.df$Summer.Freq, " (", 
																			 								vir.season.distr.y23.df$Summer.Prop, ")"),
																			 Fall = paste0(vir.season.distr.y23.df$Fall.Freq, " (", 
																			 							vir.season.distr.y23.df$Fall.Prop, ")"))

vir.season.distr.y23.df3 <- vir.season.distr.y23.df2[seq(2, nrow(vir.season.distr.y23.df2), 2),]

# Season ~ virus (+ gender + age + well/sick) - y23
# ---------------------------------------
vir.list <- virus.var.name.y23

results.df <- data.frame(OR = character(length(vir.list)),
												 p = character(length(vir.list)),
												 OR_alt = character(length(vir.list)),
												 p_alt = character(length(vir.list)),
												 stringsAsFactors=F)
for(i in 1:length(vir.list)){
	gee.out <- tryCatch(gee(as.formula(paste0('I(season2=="Winter.Spring") ~ ', 
																						vir.list[i], 
																						'+ sex + age + I(well.sick == "Sick")')),
													id = SUBJ_ID, 
													data = log.table.y23.df2, 
													family= binomial,
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=4, ncol=5))
	} else{
		s  <- summary(gee.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}
vir.season.distr.y23.df3$OR <- results.df$OR
vir.season.distr.y23.df3$p <- results.df$p

write.csv(vir.season.distr.y23.df3, "output/2018-02-12_vir.season.distr.y23.csv", quote=F,
					row.names=F)

# MPG by virus heatmap (facet by season, well/sick) - y1
# --------------------------------------

# Winter/Spring, Sick
for (i in 1:length(virus.var.name.y1)) {
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Winter.Spring" & 
			log.table.y1.df2$well.sick == "Sick"),], 
																	table(get(virus.var.name.y1[i]), MPG.name.top)))
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y1[i], ".df1"), df)
	
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Winter.Spring" & 
			log.table.y1.df2$well.sick == "Sick"),], 
																	prop.table(table(get(virus.var.name.y1[i]), MPG.name.top),1)))
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y1[i], ".df1"), df)	
}

vir.count.MPG.y1.df1 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																										 virus.simp.name.y1, ".df1")))
row.names(vir.count.MPG.y1.df1) <- NULL
vir.count.MPG.y1.df1 <- vir.count.MPG.y1.df1[seq(2, nrow(vir.count.MPG.y1.df1), 2),]
vir.count.MPG.y1.df1$well.sick <- "Sick"
vir.count.MPG.y1.df1$season2 <- "Winter.Spring"

vir.prop.MPG.y1.df1 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																										virus.simp.name.y1, ".df1")))
row.names(vir.prop.MPG.y1.df1) <- NULL
vir.prop.MPG.y1.df1 <- vir.prop.MPG.y1.df1[seq(2, nrow(vir.prop.MPG.y1.df1), 2),]
vir.prop.MPG.y1.df1$well.sick <- "Sick"
vir.prop.MPG.y1.df1$season2 <- "Winter.Spring"

# Winter/Spring, Well
for (i in 1:length(virus.var.name.y1)) {
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Winter.Spring" & 
			log.table.y1.df2$well.sick == "Well"),], 
		table(get(virus.var.name.y1[i]), MPG.name.top)))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y1[i], ".df2"), df)
	
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Winter.Spring" & 
			log.table.y1.df2$well.sick == "Well"),], 
		prop.table(table(get(virus.var.name.y1[i]), MPG.name.top),1)))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y1[i], ".df2"), df)	
}

vir.count.MPG.y1.df2 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																									virus.simp.name.y1, ".df2")))
row.names(vir.count.MPG.y1.df2) <- NULL
vir.count.MPG.y1.df2 <- vir.count.MPG.y1.df2[seq(2, nrow(vir.count.MPG.y1.df2), 2),]
vir.count.MPG.y1.df2$well.sick <- "Well"
vir.count.MPG.y1.df2$season2 <- "Winter.Spring"

vir.prop.MPG.y1.df2 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																										virus.simp.name.y1, ".df2")))
row.names(vir.prop.MPG.y1.df2) <- NULL
vir.prop.MPG.y1.df2 <- vir.prop.MPG.y1.df2[seq(2, nrow(vir.prop.MPG.y1.df2), 2),]
vir.prop.MPG.y1.df2$well.sick <- "Well"
vir.prop.MPG.y1.df2$season2 <- "Winter.Spring"

# Summer/Fall, Sick
for (i in 1:length(virus.var.name.y1)) {
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Summer.Fall" & 
			log.table.y1.df2$well.sick == "Sick"),], 
		table(get(virus.var.name.y1[i]), MPG.name.top)))
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y1[i], ".df3"), df)
	
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Summer.Fall" & 
			log.table.y1.df2$well.sick == "Sick"),], 
		prop.table(table(get(virus.var.name.y1[i]), MPG.name.top),1)))
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]), 
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y1[i], ".df3"), df)	
}

vir.count.MPG.y1.df3 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																										 virus.simp.name.y1, ".df3")))
row.names(vir.count.MPG.y1.df3) <- NULL
vir.count.MPG.y1.df3 <- vir.count.MPG.y1.df3[seq(2, nrow(vir.count.MPG.y1.df3), 2),]
vir.count.MPG.y1.df3$well.sick <- "Sick"
vir.count.MPG.y1.df3$season2 <- "Summer.Fall"

vir.prop.MPG.y1.df3 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																								 virus.simp.name.y1, ".df3")))
row.names(vir.prop.MPG.y1.df3) <- NULL
vir.prop.MPG.y1.df3 <- vir.prop.MPG.y1.df3[seq(2, nrow(vir.prop.MPG.y1.df3), 2),]
vir.prop.MPG.y1.df3$well.sick <- "Sick"
vir.prop.MPG.y1.df3$season2 <- "Summer.Fall"

# Summer/Fall, Well
for (i in 1:length(virus.var.name.y1)) {
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Summer.Fall" & 
			log.table.y1.df2$well.sick == "Well"),], 
		table(get(virus.var.name.y1[i]), MPG.name.top), useNA="always"))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]),
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y1[i], ".df4"), df)
	
	df <- as.data.frame.matrix(with(log.table.y1.df2[which(
		log.table.y1.df2$season2 == "Summer.Fall" & 
			log.table.y1.df2$well.sick == "Well"),], 
		prop.table(table(get(virus.var.name.y1[i]), MPG.name.top),1)))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y1[i]),
										 paste0(virus.real.name.y1[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y1[i], ".df4"), df)	
}

vir.count.MPG.y1.df4 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																									virus.simp.name.y1, ".df4")))
row.names(vir.count.MPG.y1.df4) <- NULL
vir.count.MPG.y1.df4 <- vir.count.MPG.y1.df4[seq(2, nrow(vir.count.MPG.y1.df4), 2),]
vir.count.MPG.y1.df4$well.sick <- "Well"
vir.count.MPG.y1.df4$season2 <- "Summer.Fall"

vir.prop.MPG.y1.df4 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																								 virus.simp.name.y1, ".df4")))
row.names(vir.prop.MPG.y1.df4) <- NULL
vir.prop.MPG.y1.df4 <- vir.prop.MPG.y1.df4[seq(2, nrow(vir.prop.MPG.y1.df4), 2),]
vir.prop.MPG.y1.df4$well.sick <- "Well"
vir.prop.MPG.y1.df4$season2 <- "Summer.Fall"

vir.count.MPG.y1.df <- do.call("rbind", list(vir.count.MPG.y1.df1, 
																					vir.count.MPG.y1.df2,
																					vir.count.MPG.y1.df3,
																					vir.count.MPG.y1.df4))
vir.count.MPG.y1.df.m <- melt(vir.count.MPG.y1.df, 
													 id.vars = c("Status", "well.sick", "season2"),
													 variable.name = "MPG.name",
													 value.name = "Count")
vir.count.MPG.y1.df.m$Vir.name <- 
	factor(vir.count.MPG.y1.df.m$Status,
				 levels=virus.real.name.y1)
vir.count.MPG.y1.df.m$Status <- NULL

vir.count.MPG.y1.df.m$season2 <- 
	factor(vir.count.MPG.y1.df.m$season2, 
				 levels=c("Winter.Spring", "Summer.Fall"))

vir.count.MPG.y1.df.m$well.sick <- 
	factor(vir.count.MPG.y1.df.m$well.sick, 
				 levels=c("Well", "Sick"))

vir.prop.MPG.y1.df <- do.call("rbind", list(vir.prop.MPG.y1.df1, 
																					vir.prop.MPG.y1.df2,
																					vir.prop.MPG.y1.df3,
																					vir.prop.MPG.y1.df4))
vir.prop.MPG.y1.df.m <- melt(vir.prop.MPG.y1.df, 
													 id.vars = c("Status", "well.sick", "season2"),
													 variable.name = "MPG.name",
													 value.name = "Prop")
vir.prop.MPG.y1.df.m$Vir.name <- 
	factor(vir.prop.MPG.y1.df.m$Status,
				 levels=virus.real.name.y1)
vir.prop.MPG.y1.df.m$Status <- NULL

vir.prop.MPG.y1.df.m$season2 <- 
	factor(vir.prop.MPG.y1.df.m$season2, 
				 levels=c("Winter.Spring", "Summer.Fall"))

vir.prop.MPG.y1.df.m$well.sick <- 
	factor(vir.prop.MPG.y1.df.m$well.sick, 
				 levels=c("Well", "Sick"))

vir.prop.MPG.y1.df.m$Count <- vir.count.MPG.y1.df.m$Count
vir.count.MPG.y1.df.m$Prop <- vir.prop.MPG.y1.df.m$Prop

pdf("output/2018-02-12_MPG.vir.season.well.sick.count.y1.pdf", width=12, height=7.5)
ggplot(vir.count.MPG.y1.df.m, aes(x=Vir.name, y=MPG.name)) + 
	geom_tile(aes(fill = Count)) + 
	geom_text(aes(label = Count)) +
	facet_grid(factor(season2, labels=c("Winter/Spring", "Summer/Fall")) ~ 
						 	factor(well.sick, labels=c("Well", "Sick"))) + 
	scale_fill_gradientn(colours = seqPalette(100)) + 
	labs(y = "MPG", x = "Virus") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("output/2018-02-12_MPG.vir.season.well.sick.prop.y1.pdf", width=12, height=7.5)
ggplot(vir.prop.MPG.y1.df.m, aes(x=Vir.name, y=MPG.name)) + 
	geom_tile(aes(fill = Prop)) + 
	geom_text(aes(label = paste0(good.signif(Prop*100,2), "%"))) +
	facet_grid(factor(season2, labels=c("Winter/Spring", "Summer/Fall"))  ~ 
						 	factor(well.sick, labels=c("Well", "Unwell"))) + 
	scale_fill_gradientn(colours = seqPalette(100), name = "Proportion \n of columns") + 
	labs(y = "MPG", x = "Virus") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("output/2018-02-12_MPG.vir.season.well.sick.prop.countlabel.y1.pdf", width=16, height=10)
ggplot(vir.prop.MPG.y1.df.m, aes(x=Vir.name, y=MPG.name)) + 
	geom_tile(aes(fill = Prop)) + 
	geom_text(aes(label = Count)) +
	facet_grid(season2 ~ well.sick) + 
	scale_fill_gradientn(colours = seqPalette(100))
dev.off()

# MPG by virus heatmap (facet by season, well/sick) - y23
# --------------------------------------

# Winter/Spring, Sick
for (i in 1:length(virus.var.name.y23)) {
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Winter.Spring" & 
			log.table.y23.df2$well.sick == "Sick"),], 
		table(get(virus.var.name.y23[i]), MPG.name.top)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y23[i], ".df1"), df)
	
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Winter.Spring" & 
			log.table.y23.df2$well.sick == "Sick"),], 
		prop.table(table(get(virus.var.name.y23[i]), MPG.name.top),1)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y23[i], ".df1"), df)	
}

vir.count.MPG.y23.df1 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																										 virus.simp.name.y23, ".df1")))
row.names(vir.count.MPG.y23.df1) <- NULL
vir.count.MPG.y23.df1 <- vir.count.MPG.y23.df1[seq(2, nrow(vir.count.MPG.y23.df1), 2),]
vir.count.MPG.y23.df1$well.sick <- "Sick"
vir.count.MPG.y23.df1$season2 <- "Winter.Spring"

vir.prop.MPG.y23.df1 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																										virus.simp.name.y23, ".df1")))
row.names(vir.prop.MPG.y23.df1) <- NULL
vir.prop.MPG.y23.df1 <- vir.prop.MPG.y23.df1[seq(2, nrow(vir.prop.MPG.y23.df1), 2),]
vir.prop.MPG.y23.df1$well.sick <- "Sick"
vir.prop.MPG.y23.df1$season2 <- "Winter.Spring"

# Winter/Spring, Well
for (i in 1:length(virus.var.name.y23)) {
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Winter.Spring" & 
			log.table.y23.df2$well.sick == "Well"),], 
		table(get(virus.var.name.y23[i]), MPG.name.top)))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y23[i], ".df2"), df)
	
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Winter.Spring" & 
			log.table.y23.df2$well.sick == "Well"),], 
		prop.table(table(get(virus.var.name.y23[i]), MPG.name.top),1)))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y23[i], ".df2"), df)	
}

vir.count.MPG.y23.df2 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																										 virus.simp.name.y23, ".df2")))
row.names(vir.count.MPG.y23.df2) <- NULL
vir.count.MPG.y23.df2 <- vir.count.MPG.y23.df2[seq(2, nrow(vir.count.MPG.y23.df2), 2),]
vir.count.MPG.y23.df2$well.sick <- "Well"
vir.count.MPG.y23.df2$season2 <- "Winter.Spring"

vir.prop.MPG.y23.df2 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																										virus.simp.name.y23, ".df2")))
row.names(vir.prop.MPG.y23.df2) <- NULL
vir.prop.MPG.y23.df2 <- vir.prop.MPG.y23.df2[seq(2, nrow(vir.prop.MPG.y23.df2), 2),]
vir.prop.MPG.y23.df2$well.sick <- "Well"
vir.prop.MPG.y23.df2$season2 <- "Winter.Spring"

# Summer/Fall, Sick
for (i in 1:length(virus.var.name.y23)) {
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Summer.Fall" & 
			log.table.y23.df2$well.sick == "Sick"),], 
		table(get(virus.var.name.y23[i]), MPG.name.top)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y23[i], ".df3"), df)
	
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Summer.Fall" & 
			log.table.y23.df2$well.sick == "Sick"),], 
		prop.table(table(get(virus.var.name.y23[i]), MPG.name.top),1)))
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]), 
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y23[i], ".df3"), df)	
}

vir.count.MPG.y23.df3 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																										 virus.simp.name.y23, ".df3")))
row.names(vir.count.MPG.y23.df3) <- NULL
vir.count.MPG.y23.df3 <- vir.count.MPG.y23.df3[seq(2, nrow(vir.count.MPG.y23.df3), 2),]
vir.count.MPG.y23.df3$well.sick <- "Sick"
vir.count.MPG.y23.df3$season2 <- "Summer.Fall"

vir.prop.MPG.y23.df3 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																										virus.simp.name.y23, ".df3")))
row.names(vir.prop.MPG.y23.df3) <- NULL
vir.prop.MPG.y23.df3 <- vir.prop.MPG.y23.df3[seq(2, nrow(vir.prop.MPG.y23.df3), 2),]
vir.prop.MPG.y23.df3$well.sick <- "Sick"
vir.prop.MPG.y23.df3$season2 <- "Summer.Fall"

# Summer/Fall, Well
for (i in 1:length(virus.var.name.y23)) {
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Summer.Fall" & 
			log.table.y23.df2$well.sick == "Well"),], 
		table(get(virus.var.name.y23[i]), MPG.name.top), useNA="always"))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]),
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.count.MPG.", virus.simp.name.y23[i], ".df4"), df)
	
	df <- as.data.frame.matrix(with(log.table.y23.df2[which(
		log.table.y23.df2$season2 == "Summer.Fall" & 
			log.table.y23.df2$well.sick == "Well"),], 
		prop.table(table(get(virus.var.name.y23[i]), MPG.name.top),1)))
	if (nrow(df) == 1 && row.names(df) == 0) {
		df <- rbind(df, rep(0, ncol(df)))
	}
	if (nrow(df) == 1 && row.names(df) == 1) {
		df <- rbind(rep(0, ncol(df)), df)
	}
	row.names(df) <- c(paste0("No ", virus.real.name.y23[i]),
										 paste0(virus.real.name.y23[i]))
	df$Status <- row.names(df)
	assign(paste0("vir.prop.MPG.", virus.simp.name.y23[i], ".df4"), df)	
}

vir.count.MPG.y23.df4 <- do.call("rbind", mget(paste0("vir.count.MPG.", 
																										 virus.simp.name.y23, ".df4")))
row.names(vir.count.MPG.y23.df4) <- NULL
vir.count.MPG.y23.df4 <- vir.count.MPG.y23.df4[seq(2, nrow(vir.count.MPG.y23.df4), 2),]
vir.count.MPG.y23.df4$well.sick <- "Well"
vir.count.MPG.y23.df4$season2 <- "Summer.Fall"

vir.prop.MPG.y23.df4 <- do.call("rbind", mget(paste0("vir.prop.MPG.", 
																										virus.simp.name.y23, ".df4")))
row.names(vir.prop.MPG.y23.df4) <- NULL
vir.prop.MPG.y23.df4 <- vir.prop.MPG.y23.df4[seq(2, nrow(vir.prop.MPG.y23.df4), 2),]
vir.prop.MPG.y23.df4$well.sick <- "Well"
vir.prop.MPG.y23.df4$season2 <- "Summer.Fall"

vir.count.MPG.y23.df <- do.call("rbind", list(vir.count.MPG.y23.df1, 
																						 vir.count.MPG.y23.df2,
																						 vir.count.MPG.y23.df3,
																						 vir.count.MPG.y23.df4))
vir.count.MPG.y23.df.m <- melt(vir.count.MPG.y23.df, 
															id.vars = c("Status", "well.sick", "season2"),
															variable.name = "MPG.name",
															value.name = "Count")
vir.count.MPG.y23.df.m$Vir.name <- 
	factor(vir.count.MPG.y23.df.m$Status,
				 levels=virus.real.name.y23)
vir.count.MPG.y23.df.m$Status <- NULL

vir.count.MPG.y23.df.m$season2 <- 
	factor(vir.count.MPG.y23.df.m$season2, 
				 levels=c("Winter.Spring", "Summer.Fall"))

vir.count.MPG.y23.df.m$well.sick <- 
	factor(vir.count.MPG.y23.df.m$well.sick, 
				 levels=c("Well", "Sick"))

vir.prop.MPG.y23.df <- do.call("rbind", list(vir.prop.MPG.y23.df1, 
																						vir.prop.MPG.y23.df2,
																						vir.prop.MPG.y23.df3,
																						vir.prop.MPG.y23.df4))
vir.prop.MPG.y23.df.m <- melt(vir.prop.MPG.y23.df, 
														 id.vars = c("Status", "well.sick", "season2"),
														 variable.name = "MPG.name",
														 value.name = "Prop")
vir.prop.MPG.y23.df.m$Vir.name <- 
	factor(vir.prop.MPG.y23.df.m$Status,
				 levels=virus.real.name.y23)
vir.prop.MPG.y23.df.m$Status <- NULL

vir.prop.MPG.y23.df.m$season2 <- 
	factor(vir.prop.MPG.y23.df.m$season2, 
				 levels=c("Winter.Spring", "Summer.Fall"))

vir.prop.MPG.y23.df.m$well.sick <- 
	factor(vir.prop.MPG.y23.df.m$well.sick, 
				 levels=c("Well", "Sick"))

vir.prop.MPG.y23.df.m$Count <- vir.count.MPG.y23.df.m$Count
vir.count.MPG.y23.df.m$Prop <- vir.prop.MPG.y23.df.m$Prop

pdf("output/2018-02-12_MPG.vir.season.well.sick.count.y23.pdf", width=12, height=7.5)
ggplot(vir.count.MPG.y23.df.m, aes(x=Vir.name, y=MPG.name)) + 
	geom_tile(aes(fill = Count)) + 
	geom_text(aes(label = Count)) +
	facet_grid(factor(season2, labels=c("Winter/Spring", "Summer/Fall")) ~ 
						 	factor(well.sick, labels=c("Well", "Sick"))) + 
	scale_fill_gradientn(colours = seqPalette(100)) + 
	labs(y = "MPG", x = "Virus") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("output/2018-02-12_MPG.vir.season.well.sick.prop.y23.pdf", width=12, height=7.5)
ggplot(vir.prop.MPG.y23.df.m, aes(x=Vir.name, y=MPG.name)) + 
	geom_tile(aes(fill = Prop)) + 
	geom_text(aes(label = paste0(good.signif(Prop*100,2), "%"))) +
	facet_grid(factor(season2, labels=c("Winter/Spring", "Summer/Fall"))  ~ 
						 	factor(well.sick, labels=c("Well", "Unwell"))) + 
	scale_fill_gradientn(colours = seqPalette(100), name = "Proportion \n of columns") + 
	labs(y = "MPG", x = "Virus") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("output/2018-02-12_MPG.vir.season.well.sick.prop.countlabel.y23.pdf", width=16, height=10)
ggplot(vir.prop.MPG.y23.df.m, aes(x=Vir.name, y=MPG.name)) + 
	geom_tile(aes(fill = Prop)) + 
	geom_text(aes(label = Count)) +
	facet_grid(season2 ~ well.sick) + 
	scale_fill_gradientn(colours = seqPalette(100))
dev.off()

save.image("output/2018-02-20_micro_assoc_analysis_vir.Rdata")

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# MPG by virus heatmap (facet by season, well/sick)
# --------------------------------------

MPG.order <- c("Haemophilus MPG2", "Streptococcus MPG1", "Haemophilus MPG1", 
							 "Moraxella MPG", "Streptococcus MPG2", "Corynebacterium MPG",
							 "Alloiococcus MPG", "Staphylococcus MPG", "Other MPG")
	
for (i in 1:length(virus.var.name)) {
	df <- as.data.frame.matrix(with(log.table.df2, 
		table(MPG.name.top, get(virus.var.name[i]))))
	colnames(df) <- c("No virus", "Virus")
	df <- df[rev(row.names(df)),]
	assign(paste0("mpg.count.", virus.simp.name[i], ".df"), df)
	
	df2 <- as.data.frame.matrix(with(log.table.df2, 
		prop.table(table(MPG.name.top, get(virus.var.name[i])),2)))
	colnames(df2) <- c("No virus", "Virus")
	df2 <- df2[rev(row.names(df2)),]
	assign(paste0("mpg.prop.", virus.simp.name[i], ".df"), df2)	
	
	df3 <- data.frame(MPG.name = row.names(df),
									 No.vir.Freq = df[, "No virus"],
										No.vir.Prop = paste0(good.signif(df2[, "No virus"]*100,2), "%"),
										Vir.Freq = df[, "Virus"],
										Vir.Prop = paste0(good.signif(df2[, "Virus"]*100,2), "%"),
										stringsAsFactors = F)
	
	df4 <- data.frame(MPG.name = df3$MPG.name,
										No.vir = paste0(df3$No.vir.Freq, " (", df3$No.vir.Prop, ")"),
										Vir = paste0(df3$Vir.Freq, " (", df3$Vir.Prop, ")"),
										stringsAsFactors = F)
	
	MPG.list <- df4$MPG.name
	
	results.df <- data.frame(OR = character(length(MPG.list)),
													 p = character(length(MPG.list)),
													 OR.ws = character(length(MPG.list)),
													 p.ws = character(length(MPG.list)),
													 OR.season = character(length(MPG.list)),
													 p.season = character(length(MPG.list)),
													 OR.ws.season = character(length(MPG.list)),
													 p.ws.season = character(length(MPG.list)),
													 stringsAsFactors=F)
	for(j in 1:length(MPG.list)){
		gee.out <- tryCatch(gee(I(get(virus.var.name[i]) == 1) ~
									 	I(MPG.name == MPG.list[j]) + gender + age.day,
									 id = subj, 
									 data = log.table.df2, 
									 family= binomial,
									 corstr = "unstructured"),
									 error = function(err) NA)
		if(is.na(gee.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=4, ncol=5))
		} else{
			s  <- summary(gee.out)
		}
		OR  <-  exp(s$coefficients[,1])
		lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
		uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
		p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
		OR <- good.signif(OR,2)
		lci <- good.signif(lci,2)
		uci <- good.signif(uci,2)
		p <- signif(p,2)
		results.df$OR[j] <- paste0(OR, " (", lci, "-", uci, ")")[2]
		results.df$p[j] <- paste0(p)[2]
	}
	df4$OR <- results.df$OR
	df4$p <- results.df$p
	
	for(j in 1:length(MPG.list)){
		gee.out <- tryCatch(gee(I(get(virus.var.name[i]) == 1) ~
															I(MPG.name == MPG.list[j]) + gender + age.day + well.sick,
														id = subj, 
														data = log.table.df2, 
														family= binomial,
														corstr = "unstructured"),
												error = function(err) NA)
		if(is.na(gee.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=5, ncol=5))
		} else{
			s  <- summary(gee.out)
		}
		OR  <-  exp(s$coefficients[,1])
		lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
		uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
		p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
		OR <- good.signif(OR,2)
		lci <- good.signif(lci,2)
		uci <- good.signif(uci,2)
		p <- signif(p,2)
		results.df$OR.ws[j] <- paste0(OR, " (", lci, "-", uci, ")")[2]
		results.df$p.ws[j] <- paste0(p)[2]
	}
	df4$OR.ws <- results.df$OR.ws
	df4$p.ws <- results.df$p.ws

	for(j in 1:length(MPG.list)){
		gee.out <- tryCatch(gee(I(get(virus.var.name[i]) == 1) ~
															I(MPG.name == MPG.list[j]) + gender + age.day + season,
														id = subj, 
														data = log.table.df2, 
														family= binomial,
														corstr = "unstructured"),
												error = function(err) NA)
		if(is.na(gee.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=5, ncol=5))
		} else{
			s  <- summary(gee.out)
		}
		OR  <-  exp(s$coefficients[,1])
		lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
		uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
		p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
		OR <- good.signif(OR,2)
		lci <- good.signif(lci,2)
		uci <- good.signif(uci,2)
		p <- signif(p,2)
		results.df$OR.season[j] <- paste0(OR, " (", lci, "-", uci, ")")[2]
		results.df$p.season[j] <- paste0(p)[2]
	}
	df4$OR.season <- results.df$OR.season
	df4$p.season <- results.df$p.season
	
	for(j in 1:length(MPG.list)){
		gee.out <- tryCatch(gee(I(get(virus.var.name[i]) == 1) ~
															I(MPG.name == MPG.list[j]) + gender + age.day + 
															well.sick + season,
														id = subj, 
														data = log.table.df2, 
														family= binomial,
														corstr = "unstructured"),
												error = function(err) NA)
		if(is.na(gee.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(gee.out)
		}
		OR  <-  exp(s$coefficients[,1])
		lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
		uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
		p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
		OR <- good.signif(OR,2)
		lci <- good.signif(lci,2)
		uci <- good.signif(uci,2)
		p <- signif(p,2)
		results.df$OR.ws.season[j] <- paste0(OR, " (", lci, "-", uci, ")")[2]
		results.df$p.ws.season[j] <- paste0(p)[2]
	}
	df4$OR.ws.season <- results.df$OR.ws.season
	df4$p.ws.season <- results.df$p.ws.season
	
	df4 <- df4[match(MPG.order, df4$MPG.name),]
	
	write.csv(df4, file=paste0("output/2018-02-13_mpg.", virus.simp.name[i], 
														 ".csv"), quote=F, row.names=F)
}

# Well/sick ~ virus * MPG * season + gender + age
# ---------------------------------------
risk.MPG <- c("Haemophilus MPG2", "Streptococcus MPG1", "Haemophilus MPG1", 
							"Moraxella MPG")

results.df <- data.frame(Variable = c("Any virus", "Any risk MPG", 
																			"Season = Winter/Spring", "Gender = Male",
																			"Age", "Virus x MPG", "Virus x Season", 
																			"MPG x Season", "Virus x MPG x Season"),
												 OR = character(9),
												 p = character(9))
												 
gee.out <- gee(I(well.sick == "Sick") ~
							 	I(virus.any == 1) * I(MPG.name.top %in% risk.MPG) * 
							 	I(season2 == "Winter.Spring") + gender + age.day,
							 id = subj, 
							 data = log.table.df2, 
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.df$OR <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.df$p <- paste0(p)[-1]

write.csv(results.df, "output/2018-02-13_virus.risk.mpg.season.interact.csv", 
					quote=F, row.names=F)

protective.MPG <- c("Staphylococcus MPG", "Streptococcus MPG2", "Alloiococcus MPG", 
							"Corynebacterium MPG")

results.df <- data.frame(Variable = c("Any virus", "Any protective MPG", 
																			"Season = Winter/Spring", "Gender = Male",
																			"Age", "Virus x MPG", "Virus x Season", 
																			"MPG x Season", "Virus x MPG x Season"),
												 OR = character(9),
												 p = character(9))

gee.out <- gee(I(well.sick == "Sick") ~
							 	I(virus.any == 1) * I(MPG.name.top %in% protective.MPG) * 
							 	I(season2 == "Winter.Spring") + gender + age.day,
							 id = subj, 
							 data = log.table.df2, 
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.df$OR <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.df$p <- paste0(p)[-1]

write.csv(results.df, "output/2018-02-13_virus.protect.mpg.season.interact.csv", 
					quote=F, row.names=F)


# Well/sick ~ rhinovirus * MPG * season + gender + age
# ---------------------------------------
results.df <- data.frame(Variable = c("RV", "Any risk MPG", 
																			"Season = Winter/Spring", "Gender = Male",
																			"Age", "Virus x MPG", "Virus x Season", 
																			"MPG x Season", "Virus x MPG x Season"),
												 OR = character(9),
												 p = character(9))

gee.out <- gee(I(well.sick == "Sick") ~
							 	I(RV.new == 1) * I(MPG.name.top %in% risk.MPG) * 
							 	I(season2 == "Winter.Spring") + gender + age.day,
							 id = subj, 
							 data = log.table.df2, 
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.df$OR <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.df$p <- paste0(p)[-1]

write.csv(results.df, "output/2018-02-13_RV.risk.mpg.season.interact.csv", 
					quote=F, row.names=F)

results.df <- data.frame(Variable = c("RV", "Any protective MPG", 
																			"Season = Winter/Spring", "Gender = Male",
																			"Age", "Virus x MPG", "Virus x Season", 
																			"MPG x Season", "Virus x MPG x Season"),
												 OR = character(9),
												 p = character(9))

gee.out <- gee(I(well.sick == "Sick") ~
							 	I(RV.new == 1) * I(MPG.name.top %in% protective.MPG) * 
							 	I(season2 == "Winter.Spring") + gender + age.day,
							 id = subj, 
							 data = log.table.df2, 
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.df$OR <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.df$p <- paste0(p)[-1]

write.csv(results.df, "output/2018-02-13_RV.protect.mpg.season.interact.csv", 
					quote=F, row.names=F)


results.df <- data.frame(Variable = c("Any virus", "Any protective MPG", 
																			"Season = Winter/Spring", 
																			"Any risk MPG",
																			"Gender = Male",
																			"Age", 
																			"Virus x protective MPG", 
																			"Virus x Season", 
																			"Protective MPG x Season", 
																			"Virus x risk MPG",
																			"Season x risk MPG", 
																			"Virus x protective MPG x Season",
																			"Virus x risk MPG x Season"),
												 OR = character(13),
												 p = character(13))

gee.out <- gee(I(well.sick == "Sick") ~
							 	I(virus.any == 1) * I(MPG.name.top %in% protective.MPG) * 
							 	I(season2 == "Winter.Spring") + 
							 	I(virus.any == 1) * I(MPG.name.top %in% risk.MPG) * 
							 	I(season2 == "Winter.Spring") + 
							 	gender + age.day,
							 id = subj, 
							 data = log.table.df2, 
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.df$OR <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.df$p <- paste0(p)[-1]

write.csv(results.df, "output/2018-02-13_virus.all.mpg.season.interact.csv", 
					quote=F, row.names=F)


results.df <- data.frame(Variable = c("RV", "Any protective MPG", 
																			"Season = Winter/Spring", 
																			"Any risk MPG",
																			"Gender = Male",
																			"Age", 
																			"Virus x protective MPG", 
																			"Virus x Season", 
																			"Protective MPG x Season", 
																			"Virus x risk MPG",
																			"Season x risk MPG", 
																			"Virus x protective MPG x Season",
																			"Virus x risk MPG x Season"),
												 OR = character(13),
												 p = character(13))

gee.out <- gee(I(well.sick == "Sick") ~
							 	I(RV.new == 1) * I(MPG.name.top %in% protective.MPG) * 
							 	I(season2 == "Winter.Spring") + 
							 	I(RV.new == 1) * I(MPG.name.top %in% risk.MPG) * 
							 	I(season2 == "Winter.Spring") + 
							 	gender + age.day,
							 id = subj, 
							 data = log.table.df2, 
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.df$OR <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.df$p <- paste0(p)[-1]

write.csv(results.df, "output/2018-02-13_RV.all.mpg.season.interact.csv", 
					quote=F, row.names=F)


# All
# ---------------------------------------

results.df <- data.frame(
	Variable = c("RV", "RSV", "Influenza", "Adenovirus", "Coronavirus", "Parainfluenza",
							 "Corynebacterium MPG", "Alloiococcus MPG", "Staphylococcus MPG", 
							 "Streptococcus MPG2", "Streptococcus MPG1", "Moraxella MPG", 
							 "Haemophilus MPG2", "Haemophilus MPG1", "Season = Winter/Spring", 
							 "Gender = Male", "Age"))
gee.out <- gee(I(well.sick == "Sick") ~
							 	RV.new + RSV + Influenza + Adenovirus + 
							 	Coronavirus + Parainfluenza + MPG.name.top + 
							 	I(season2 == "Winter.Spring") + gender + age.day,
							 id = subj,
							 data = log.table.df2,
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.df$OR <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.df$p <- paste0(p)[-1]

# Winter/Spring only
# ---------------------------------------
results.ws.df <- data.frame(
	Variable = c("RV", "RSV", "Influenza", "Adenovirus", "Coronavirus", "Parainfluenza",
							 "Corynebacterium MPG", "Alloiococcus MPG", "Staphylococcus MPG", 
							 "Streptococcus MPG2", "Streptococcus MPG1", "Moraxella MPG", 
							 "Haemophilus MPG2", "Haemophilus MPG1", "Gender = Male", "Age"))
gee.out <- gee(I(well.sick == "Sick") ~
							 	RV.new + RSV + Influenza + Adenovirus +  
							 	Coronavirus + Parainfluenza + MPG.name.top + 
							 	gender + age.day,
							 id = subj,
							 data = log.table.df2[which(log.table.df2$season2 == "Winter.Spring"),],
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.ws.df$OR.ws <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.ws.df$p.ws <- paste0(p)[-1]


# Summer/Fall only
# ---------------------------------------
results.sf.df <- data.frame(
	Variable = c("RV", "Influenza", "Adenovirus", "Coronavirus", "Parainfluenza",
							 "Corynebacterium MPG", "Alloiococcus MPG", "Staphylococcus MPG", 
							 "Streptococcus MPG2", "Streptococcus MPG1", "Moraxella MPG", 
							 "Haemophilus MPG2", "Haemophilus MPG1", "Gender = Male", "Age"))
gee.out <- gee(I(well.sick == "Sick") ~
							 	RV.new + Influenza + Adenovirus +  
							 	Coronavirus + Parainfluenza + MPG.name.top + 
							 	gender + age.day,
							 id = subj,
							 data = log.table.df2[which(log.table.df2$season2 == "Summer.Fall"),],
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.sf.df$OR.sf <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.sf.df$p.sf <- paste0(p)[-1]

# Absent viruses
# ---------------------------------------

results.no.maj.vir.df <- data.frame(
	Variable = c("Corynebacterium MPG", "Alloiococcus MPG", 
							 "Streptococcus MPG2", "Streptococcus MPG1", "Moraxella MPG", 
							 "Haemophilus MPG2", "Haemophilus MPG1", 
							 "Season = Winter/Spring", "Gender = Male", "Age"))
df <- log.table.df2[which(log.table.df2$virus.any == 0 &
												 	log.table.df2$MPG.name.top != "Staphylococcus MPG"),]
levels(df$MPG.name.top) <- c("Other MPG", "Corynebacterium MPG", "Alloiococcus MPG", NA, 
														 "Streptococcus MPG2", "Streptococcus MPG1", "Moraxella MPG", "Haemophilus MPG2",
														 "Haemophilus MPG1")
gee.out <- gee(I(well.sick == "Sick") ~
							 	MPG.name.top +
							 	(season2 == "Winter.Spring") + gender + age.day,
							 id = subj,
							 data = df,
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.no.maj.vir.df$OR.no.maj.vir <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.no.maj.vir.df$p.no.maj.vir <- paste0(p)[-1]

# Present viruses
# ---------------------------------------
results.maj.vir.df <- data.frame(
	Variable = c("Corynebacterium MPG", "Alloiococcus MPG", "Staphylococcus MPG", 
							 "Streptococcus MPG2", "Streptococcus MPG1", "Moraxella MPG", 
							 "Haemophilus MPG2", "Haemophilus MPG1", 
							 "Season = Winter/Spring", "Gender = Male", "Age"))

gee.out <- gee(I(well.sick == "Sick") ~
							 	MPG.name.top +
							 	(season2 == "Winter.Spring") + gender + age.day,
							 id = subj,
							 data =  log.table.df2[which(log.table.df2$virus.any == 1),],
							 family= binomial,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
results.maj.vir.df$OR.maj.vir <- paste0(OR, " (", lci, "-", uci, ")")[-1]
results.maj.vir.df$p.maj.vir <- paste0(p)[-1]

results.all.df <- merge(results.df, results.ws.df, by="Variable", all.x = T)
results.all.df <- merge(results.all.df, results.sf.df, by="Variable", all.x = T)
results.all.df <- merge(results.all.df, results.no.maj.vir.df, by="Variable", all.x = T)
results.all.df <- merge(results.all.df, results.maj.vir.df, by="Variable", all.x = T)

Variable.order <- c("RV", "RSV", "Influenza", "Adenovirus", "Coronavirus", "Parainfluenza",
										"Haemophilus MPG2", "Streptococcus MPG1", "Haemophilus MPG1", "Moraxella MPG",
										"Streptococcus MPG2", "Alloiococcus MPG", "Corynebacterium MPG", "Staphylococcus MPG",
										"Season = Winter/Spring", "Gender = Male", "Age")

results.all.df <- results.all.df[match(Variable.order, results.all.df$Variable),]

write.csv(results.all.df, "output/2018-02-13_well.sick.gee.results.csv", quote=F,
					row.names=F)


# Correlation matrix of viruses
# =====

x.df <- with(log.table.df2, 
						 data.frame(RV.new, RSV, Parainfluenza, Coronavirus, Adenovirus, 
						 					 Bocavirus, Influenza, MPV, Enterovirus, Echovirus, 
						 					 Coxsackie, TTV, HSV))
vir.cor.mat <- cor(x.df, method="spearman")
vir.cor.mat[upper.tri(vir.cor.mat)] <- NA
vir.cor.m.df <- melt(vir.cor.mat)

vir.cor.test.mat <- corr.test(x.df, method="spearman", adjust = "none")$p
vir.cor.test.mat[upper.tri(vir.cor.test.mat)] <- NA
vir.cor.test.mat.m.df <- melt(vir.cor.test.mat)
colnames(vir.cor.test.mat.m.df) <- c("Var1", "Var2", "p")

vir.cor.m.df <- merge(vir.cor.m.df, vir.cor.test.mat.m.df,
											by=c("Var1", "Var2"))
vir.cor.m.df$p.star <- ifelse(vir.cor.m.df$p < 0.0001, 
															"***", ifelse(vir.cor.m.df$p < 0.001, 
																						"**", ifelse(vir.cor.m.df$p < 0.05,
																												 "*", "")))
vir.cor.m.df$p.factor <- factor(vir.cor.m.df$p.star, 
																levels=c("", "*", "**", "***"),
																labels=c("0", "1", "2", "3"))
levels(vir.cor.m.df$p.factor) <- c("NS", "p<0.05", "p<0.05", "p<0.05")

levels(vir.cor.m.df$Var1)[which(levels(vir.cor.m.df$Var1) == "RV.new")] <- "RV"
levels(vir.cor.m.df$Var2)[which(levels(vir.cor.m.df$Var2) == "RV.new")] <- "RV"

pdf("output/2018-02-14_virus_corr.pdf", width=12, height=7.5)
ggplot(vir.cor.m.df, aes(x=Var1, y=Var2)) + 
	geom_tile(aes(fill = value)) + 
	scale_fill_gradientn(colours = divPalette(100), 
											 limit = c(-1,1)) + 
	geom_text(aes(label = signif(value,2), color=p.factor)) +
	scale_color_manual(values = c("black", "red")) +
	labs(y = "Virus", x = "Virus",
			 fill = "Rho \n(Spearman correlation)",
			 color = "Significance") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# For each column, count number of rows in each other column that also contain 1
virus.var.list <- c("RV.new", "RSV", "Parainfluenza", "Coronavirus", "Adenovirus", 
										"Bocavirus", "Influenza", "MPV", "Enterovirus", "Echovirus", 
										"Coxsackie", "TTV", "HSV")
for (i in 1:length(virus.var.list)) {
	vir.target <- virus.var.list[i]
	y.df <- x.df %>%
		filter(get(vir.target) == 1) %>%
		melt(variable.names = "virus",
				 value.name = "count") %>%
		table() %>%
		as.data.frame.matrix()
	colnames(y.df) <- c("No", paste0(vir.target))
	y.df$virus <- row.names(y.df)
	y.df$No <- NULL
	assign(paste0(vir.target, ".count.df"), y.df)
}

vir.all.count.df <- merge(RV.new.count.df, RSV.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Parainfluenza.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Coronavirus.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Adenovirus.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Bocavirus.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Influenza.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, MPV.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Enterovirus.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Echovirus.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, Coxsackie.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, TTV.count.df, by="virus", all=T)
vir.all.count.df <- merge(vir.all.count.df, HSV.count.df, by="virus", all=T)

vir.all.count.df <- vir.all.count.df[match(virus.var.list, vir.all.count.df$virus),]
vir.all.prop.df <- vir.all.count.df
vir.all.prop.df[, -1] <- vir.all.count.df[, -1]/colSums(vir.all.count.df[, -1])

vir.all.count.mat <- as.matrix(vir.all.count.df[,-1])
vir.all.count.mat[upper.tri(vir.all.count.mat)] <- NA
vir.all.count.df[,-1] <- vir.all.count.mat

vir.all.count.df.m <- melt(vir.all.count.df, id.vars = "virus")
vir.all.count.df.m$virus <- factor(vir.all.count.df.m$virus, 
																	 levels = virus.var.list)
vir.all.count.df.m$variable <- factor(vir.all.count.df.m$variable, 
																	 levels = virus.var.list)
levels(vir.all.count.df.m$virus)[
	which(levels(vir.all.count.df.m$virus) == "RV.new")] <- "RV"
levels(vir.all.count.df.m$variable)[
	which(levels(vir.all.count.df.m$variable) == "RV.new")] <- "RV"

pdf("output/2018-02-14_virus_overlap_count.pdf", width=12, height=7.5)
ggplot(vir.all.count.df.m, aes(x=virus, y=variable)) + 
	geom_tile(aes(fill = value)) + 
	geom_text(aes(label = value)) +
	scale_fill_gradientn(colours = seqPalette(100)) + 
	labs(y = "Virus", x = "Virus") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

vir.all.prop.df.m <- melt(vir.all.prop.df, id.vars = "virus")
vir.all.prop.df.m$virus <- factor(vir.all.prop.df.m$virus, 
																	 levels = virus.var.list)
vir.all.prop.df.m$variable <- factor(vir.all.prop.df.m$variable, 
																			levels = virus.var.list)
levels(vir.all.prop.df.m$virus)[
	which(levels(vir.all.prop.df.m$virus) == "RV.new")] <- "RV"
levels(vir.all.prop.df.m$variable)[
	which(levels(vir.all.prop.df.m$variable) == "RV.new")] <- "RV"

pdf("output/2018-02-14_virus_overlap_prop.pdf", width=12, height=7.5)
ggplot(vir.all.prop.df.m, aes(x=virus, y=variable)) + 
	geom_tile(aes(fill = value)) + 
	geom_text(aes(label = paste0(signif(value*100,2), "%"))) +
	scale_fill_gradientn(colours = seqPalette(100),
											 limit = c(0,1)) + 
	labs(y = "Virus", x = "Virus",
			 fill = "Proportion \n of column") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# Michael Evans analysis check
# =====
with(log.table.df2[which(log.table.df2$age.day <= 365 &
													log.table.df2$RV.new == 1),],
		 table(visittype %in% c("Well", "Well/Sick"), visittype2 == "Sick"))

# In first year, among RV-positive samples,
# 39 of 348+39 = 387 (10%) of routine samples were actually sick (score > 5)
# 39 does not match 25 HRV MSIs in scheduled visits given in Lee 2012
# 260+39 = 299 does not match 214 HRV MSIs in total

log.table.df2.origin <- sample_data(coast.biom)
log.table.df2.origin$RV <- as.numeric(grepl("Rhino", log.table.df2.origin$virus))
log.table.df2.origin$RV.new <- 
	as.numeric(log.table.df2.origin$RV == 1 | !is.na(log.table.df2.origin$RV.type))
with(log.table.df2.origin[which(log.table.df2.origin$age.day <= 365 &
															 	log.table.df2.origin$RV.new == 1),],
		 table(visittype %in% c("Well", "Well/Sick"), visittype2 == "Sick"))

# Does not match either

with(log.table.df2,
		 table(visittype %in% c("Well", "Well/Sick"), visittype2 == "Sick",
		 season, RV.new))

with(log.table.df2[which(log.table.df2$age.day <= 365),],
		 table(visittype %in% c("Well", "Well/Sick"), visittype2 == "Sick",
		 			season, RV.new))

with(log.table.df2,
		 table(visittype, visittype2))

with(log.table.df2,
		 table(visittype %in% c("Well", "Well/Sick"), visittype2 == "Sick"))


