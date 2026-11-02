# Create phenotype file
# =====

load("data/2016-12-05_npem_global_heatmaps.Rdata")
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
sex.df$sex[which(sex.df$sex == 1)] <- 0
sex.df$sex[which(sex.df$sex == 2)] <- 1
write.table(sex.df, file="output/cas_imp_3_sex.phe", quote=F, 
						col.names=T, row.names=F)

# wlri serial test
# ----
library(reshape2)
wlri_num.df <- with(premaster_pheno.df, data.frame(id = SUBJ_ID,
																									sex = sex,
																									wlri_1 = wlri_1,
																									wlri_2 = wlri_2,
																									wlri_3 = wlri_3))
wlri_num.df$sex[which(is.na(wlri_num.df$sex))] <- NA
wlri_num.df$sex[which(wlri_num.df$sex == 1)] <- 0
wlri_num.df$sex[which(wlri_num.df$sex == 2)] <- 1
wlri_num_melt.df <- melt(wlri_num.df, id.vars = c("id", "sex"), 
												 variable.name = "age",
												 value.name = "y")
wlri_num_melt.df$age <- factor(wlri_num_melt.df$age, 
															 labels = c(1,2,3),
															 levels = c("wlri_1", "wlri_2", "wlri_3"))
write.table(wlri_num_melt.df, file="output/wlri_num_test.phe", quote=F, 
						col.names=T, row.names=F)
						
# hrva_lri serial test
hrva_lri_num.df <- with(premaster_pheno.df, data.frame(id = SUBJ_ID,
																									 sex = sex,
																									 hrva_lri_1 = hrva_lri_1,
																									 hrva_lri_2 = hrva_lri_2,
																									 hrva_lri_3 = hrva_lri_3))
hrva_lri_num.df$sex[which(is.na(hrva_lri_num.df$sex))] <- NA
hrva_lri_num.df$sex[which(hrva_lri_num.df$sex == 1)] <- 0
hrva_lri_num.df$sex[which(hrva_lri_num.df$sex == 2)] <- 1
hrva_lri_num_melt.df <- melt(hrva_lri_num.df, id.vars = c("id", "sex"), 
												 variable.name = "age",
												 value.name = "y")
hrva_lri_num_melt.df$age <- factor(hrva_lri_num_melt.df$age, 
															 labels = c(1,2,3),
															 levels = c("hrva_lri_1", "hrva_lri_2", "hrva_lri_3"))
write.table(hrva_lri_num_melt.df, file="output/hrva_lri_num_test.phe", quote=F, 
						col.names=T, row.names=F)

# lri serial test
lri_num.df <- with(premaster_pheno.df, data.frame(id = SUBJ_ID,
																											 sex = sex,
																											 lri_1 = lri_1,
																											 lri_2 = lri_2,
																											 lri_3 = lri_3))
lri_num.df$sex[which(is.na(lri_num.df$sex))] <- NA
lri_num.df$sex[which(lri_num.df$sex == 1)] <- 0
lri_num.df$sex[which(lri_num.df$sex == 2)] <- 1
lri_num_melt.df <- melt(lri_num.df, id.vars = c("id", "sex"), 
														 variable.name = "age",
														 value.name = "y")
lri_num_melt.df$age <- factor(lri_num_melt.df$age, 
																	 labels = c(1,2,3),
																	 levels = c("lri_1", "lri_2", "lri_3"))
write.table(lri_num_melt.df, file="output/lri_num_test.phe", quote=F, 
						col.names=T, row.names=F)

# crwz serial test
crwz_bin.df <- with(premaster_pheno.df, data.frame(id = SUBJ_ID,
																									sex = sex,
																									crwz_1 = crwz1,
																									crwz_2 = crwz2,
																									crwz_3 = crwz3,
																									crwz_4 = crwz4,
																									crwz_5 = crwz5))
crwz_bin.df$sex[which(is.na(crwz_bin.df$sex))] <- NA
crwz_bin.df$sex[which(crwz_bin.df$sex == 1)] <- 0
crwz_bin.df$sex[which(crwz_bin.df$sex == 2)] <- 1
crwz_bin_melt.df <- melt(crwz_bin.df, id.vars = c("id", "sex"), 
												variable.name = "age",
												value.name = "y")
crwz_bin_melt.df$age <- factor(crwz_bin_melt.df$age, 
															labels = c(1,2,3,4,5),
															levels = c("crwz_1", "crwz_2", "crwz_3", "crwz_4", 
																				 "crwz_5"))
write.table(crwz_bin_melt.df, file="output/crwz_bin_test.phe", quote=F, 
						col.names=T, row.names=F)