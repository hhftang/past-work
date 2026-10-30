# Create plink phenotype files for repeatabel

# data.frame: id, sex, age, y
rm(list=ls())
load("output/2018-02-20_micro_assoc_analysis_vir.Rdata")
	
	example_pheno.df <- read.table(file="data/example.phe")
	colnames(example_pheno.df) <- c("FamID", "IndID", "Phenotype")
	premaster_pheno.df <- premaster.df[which(premaster.df$SUBJ_ID %in% 
																					 	example_pheno.df$FamID),]
	master_pheno.df <- master.df[which(master.df$SUBJ_ID %in%
																		 	example_pheno.df$FamID),]

# phenotype data for genabel format
# ----

sex.df <- with(premaster_pheno.df, data.frame(id = SUBJ_ID,
																							sex = sex))
sex.df$sex[which(is.na(sex.df$sex))] <- NA
# Female: CAS dataset 1, Genabel 0, plink 2
sex.df$sex[which(sex.df$sex == 1)] <- 0
# Male: CAS dataset 2, Genabel 1, Plink 1
sex.df$sex[which(sex.df$sex == 2)] <- 1
write.table(sex.df, file="output/cas_imp_3_sex.phe", quote=F, 
						col.names=T, row.names=F)

# major otu phenotypes
# ----

# Alloiococcus..dd2e15140d20c84712456e4c12ece08b
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
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
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
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
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
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
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
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
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
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
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
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
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
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
df0 <- with(log.table.df, data.frame(id = SUBJ_ID, 
																		 LogNo, age))
df1 <- merge(df0, sex.df, by="id", all.y=T)
df2 <- with(log.esv.df, data.frame(LogNo, Streptococcus..4060107ffdf52e44f9d72343b7332609))

df3 <- merge(df1, df2, by="LogNo", all.x=T)

df3$sex[which(is.na(df3$sex))] <- NA
df3$LogNo <- NULL

colnames(df3)[4] <- "y"
df3 <- df3[order(df3$id, df3$age),]

df4 <- with(df3, data.frame(id, sex, age, y))
write.table(df4, file="output/streptococcus.4060.phe", quote=F, 
						col.names=T, row.names=F)

# major MPG phenotypes - proportions
# ----
df <- log.table.df %>% 
	filter(!is.na(year)) %>%
	group_by(year, SUBJ_ID, MPG.name) %>%
	summarize(count = n()) %>%
	mutate (prop = count/sum(count))

# Alloiococcus
df1 <- df %>% filter(MPG.name == "Alloiococcus MPG") %>%
	rename(id = SUBJ_ID,
			 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/alloiococcus.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Corynebacterium
df1 <- df %>% filter(MPG.name == "Corynebacterium MPG") %>%
	rename(id = SUBJ_ID,
				 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/corynebacterium.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Haemophilus MPG1
df1 <- df %>% filter(MPG.name == "Haemophilus MPG1") %>%
	rename(id = SUBJ_ID,
				 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/haemophilus.mpg1.phe", quote=F, 
						col.names=T, row.names=F)

# Haemophilus MPG2
df1 <- df %>% filter(MPG.name == "Haemophilus MPG2") %>%
	rename(id = SUBJ_ID,
				 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/haemophilus.mpg2.phe", quote=F, 
						col.names=T, row.names=F)

# Moraxella MPG
df1 <- df %>% filter(MPG.name == "Moraxella MPG") %>%
	rename(id = SUBJ_ID,
				 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/moraxella.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Staphylococcus MPG
df1 <- df %>% filter(MPG.name == "Staphylococcus MPG") %>%
	rename(id = SUBJ_ID,
				 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/staphylococcus.mpg.phe", quote=F, 
						col.names=T, row.names=F)

# Streptococcus MPG1
df1 <- df %>% filter(MPG.name == "Streptococcus MPG1") %>%
	rename(id = SUBJ_ID,
				 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/streptococcus.mpg1.phe", quote=F, 
						col.names=T, row.names=F)

# Streptococcus MPG2 / Mixed MPG
df1 <- df %>% filter(MPG.name == "Mixed MPG") %>%
	rename(id = SUBJ_ID,
				 age = year)
df2 <- merge(df1, sex.df, by="id", all.y=T)
df3 <- with(df2, data.frame(id, sex, age, y = prop))
df3 <- df3[order(df3$id, df3$age),]

write.table(df3, file="output/mixed.mpg.phe", quote=F, 
						col.names=T, row.names=F)

