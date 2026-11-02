# ==============================================================================
# 2018-03-07 - Create plink phenotype files for repeatabel
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

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
load("output/2018-03-07_coast_data_for_ober.Rdata")

# ------------------------------------------------------------------------------
# Generate phenotype files
# ------------------------------------------------------------------------------
# phenotype data in genabel format
# data.frame format: id, sex, age, y

# Load sex phenotype file
# ---------------------------------------
sex.df <- with(log.table.df.trunc, data.frame(id = subj, sex = as.character(gender),
																				stringsAsFactors = F))
sex.df <- unique(sex.df)
sex.df$sex[which(is.na(sex.df$sex))] <- NA
# Female: Genabel 0, plink 2
sex.df$sex[which(sex.df$sex == "Female")] <- 0
# Male: Genabel 1, Plink 1
sex.df$sex[which(sex.df$sex == "Male")] <- 1

write.table(sex.df, file="output/coast_sex.phe", quote=F, 
						col.names=T, row.names=F)

# major otu phenotypes
# ---------------------------------------

# Alloiococcus..dd2e15140d20c84712456e4c12ece08b
df0 <- with(log.table.df.trunc, data.frame(id = subj, 
																		 LogNo, age = age.day))
df1 <- merge(df0, sex.df, by="id")
df2 <- with(log.esv.df, data.frame(LogNo, Alloiococcus..dd2e15140d20c84712456e4c12ece08b))
df3 <- merge(df1, df2, by="LogNo", all.x=T)
df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/alloiococcus.phe", quote=F, 
						col.names=T, row.names=F)

# Corynebacterium..cb50e79d177c169ec5b9d5e7d308874eb
df2 <- with(log.esv.df, data.frame(LogNo, Corynebacterium..cb50e79d177c169ec5b9d5e7d308874e))
df3 <- merge(df1, df2, by="LogNo", all.x=T)
df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/corynebacterium.phe", quote=F, 
						col.names=T, row.names=F)

# Haemophilus.influenzae..bc0d904a0469d05c65b270998488ee91
df2 <- with(log.esv.df, 
						data.frame(LogNo,
											 Haemophilus.influenzae..bc0d904a0469d05c65b270998488ee91))
df3 <- merge(df1, df2, by="LogNo", all.x=T)
df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/haemophilus.bc0d.phe", quote=F, 
						col.names=T, row.names=F)

# Haemophilus..f579ddd8d7c6a3d812ca33d72110f0d6
df2 <- with(log.esv.df, 
						data.frame(LogNo,
											 Haemophilus..f579ddd8d7c6a3d812ca33d72110f0d6))
df3 <- merge(df1, df2, by="LogNo", all.x=T)
df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/haemophilus.f579.phe", quote=F, 
						col.names=T, row.names=F)

# Moraxella..d253ca966efcd811d057dc2aa6430774
df2 <- with(log.esv.df, data.frame(LogNo, Moraxella..d253ca966efcd811d057dc2aa6430774))
df3 <- merge(df1, df2, by="LogNo", all.x=T)

df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/moraxella.phe", quote=F, 
						col.names=T, row.names=F)

# Staphylococcus..29eb8e14ce61e9c4549de0f9cf200c19
df2 <- with(log.esv.df, data.frame(LogNo, Staphylococcus..29eb8e14ce61e9c4549de0f9cf200c19))
df3 <- merge(df1, df2, by="LogNo", all.x=T)
df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/staphylococcus.phe", quote=F, 
						col.names=T, row.names=F)

# Streptococcus..357551f644064bb2b49614b1caab2fa5
df2 <- with(log.esv.df, data.frame(LogNo, Streptococcus..357551f644064bb2b49614b1caab2fa5))
df3 <- merge(df1, df2, by="LogNo", all.x=T)
df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/streptococcus.3575.phe", quote=F, 
						col.names=T, row.names=F)

# Streptococcus..4060107ffdf52e44f9d72343b7332609
df2 <- with(log.esv.df, data.frame(LogNo, Streptococcus..4060107ffdf52e44f9d72343b7332609))
df3 <- merge(df1, df2, by="LogNo", all.x=T)
df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL
colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]
df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/streptococcus.4060.phe", quote=F, 
						col.names=T, row.names=F)

# Other outputs
# ---------------------------------------
write.csv(log.esv.df, "output/coast_otu_table_major_otus.csv", quote=F,
					row.names=F)

# major MPG phenotypes - proportions
# ---------------------------------------
df <- log.table.df.trunc %>% 
	filter(!is.na(age.year)) %>%
	group_by(age.year, subj, MPG.name) %>%
	summarize(count = n()) %>%
	mutate (prop = count/sum(count))

# Alloiococcus
df1 <- df %>% filter(MPG.name == "Alloiococcus MPG") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/alloiococcus.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Corynebacterium
df1 <- df %>% filter(MPG.name == "Corynebacterium MPG") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/corynebacterium.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Haemophilus MPG1
df1 <- df %>% filter(MPG.name == "Haemophilus MPG1") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/haemophilus.mpg1.phe", quote=F, 
						col.names=T, row.names=F)

# Haemophilus MPG2
df1 <- df %>% filter(MPG.name == "Haemophilus MPG2") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/haemophilus.mpg2.phe", quote=F, 
						col.names=T, row.names=F)

# Moraxella MPG
df1 <- df %>% filter(MPG.name == "Moraxella MPG") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/moraxella.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Staphylococcus MPG
df1 <- df %>% filter(MPG.name == "Staphylococcus MPG") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/staphylococcus.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Streptococcus MPG1
df1 <- df %>% filter(MPG.name == "Streptococcus MPG1") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/streptococcus.mpg1.phe", quote=F, 
						col.names=T, row.names=F)

# Streptococcus MPG2 / Mixed MPG
df1 <- df %>% filter(MPG.name == "Streptococcus MPG2") %>%
	rename(id = subj,
				 age = age.year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]
write.table(df3, file="output/streptococcus.mpg2.phe", quote=F, 
						col.names=T, row.names=F)
