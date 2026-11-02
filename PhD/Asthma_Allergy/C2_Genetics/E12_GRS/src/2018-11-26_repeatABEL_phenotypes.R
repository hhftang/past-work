# Get phenotype files for repeatABEL
# ==============================================================================
# id, sex, age, y

rm(list=ls())

# npem cluster data
load("data/2017-07-06_npem_global_heatmaps.Rdata")

# Microbiome data
load("data/2018-06-28_MPG_illness_assoc_cas.Rdata")

cas.phe.df <- read.table("data/cas_imp_3_sex.phe", header=T, stringsAsFactors = F)
cas.phe.df

library(dplyr)
library(reshape2)

# Get crwz
crwz.df <- master.df %>% select(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, crwz5) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(1,2,3,4,5)))))
																							
crwz.phe.df <- merge(cas.phe.df, crwz.df, 
										 by="id", all.x=T)

write.table(crwz.phe.df, "output/pheno/crwz.phe", row.names=F, quote=F)

# Get cras
cras.df <- master.df %>% select(SUBJ_ID, cras_3, cras_4, cras_5) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(3,4,5)))))

cras.phe.df <- merge(cas.phe.df, cras.df, 
										 by="id", all.x=T)

write.table(cras.phe.df, "output/pheno/cras.phe", row.names=F, quote=F)


# Get lri
lri.df <- master.df %>% select(SUBJ_ID, lri_1, lri_2, lri_3) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(1,2,3)))))

lri.phe.df <- merge(cas.phe.df, lri.df, 
										 by="id", all.x=T)

write.table(lri.phe.df, "output/pheno/lri.phe", row.names=F, quote=F)

# Get flri
flri.df <- master.df %>% select(SUBJ_ID, flri_1, flri_2, flri_3) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(1,2,3)))))

flri.phe.df <- merge(cas.phe.df, flri.df, 
										by="id", all.x=T)

write.table(flri.phe.df, "output/pheno/flri.phe", row.names=F, quote=F)

# Get wlri
wlri.df <- master.df %>% select(SUBJ_ID, wlri_1, wlri_2, wlri_3) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(1,2,3)))))

wlri.phe.df <- merge(cas.phe.df, wlri.df, 
										 by="id", all.x=T)

write.table(wlri.phe.df, "output/pheno/wlri.phe", row.names=F, quote=F)

# Get hrva_lri
hrva_lri.df <- master.df %>% select(SUBJ_ID, hrva_lri_1, hrva_lri_2, hrva_lri_3) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(1,2,3)))))

hrva_lri.phe.df <- merge(cas.phe.df, hrva_lri.df, 
										 by="id", all.x=T)

write.table(hrva_lri.phe.df, "output/pheno/hrva_lri.phe", row.names=F, quote=F)

# Get hrvc_lri
hrvc_lri.df <- master.df %>% select(SUBJ_ID, hrvc_lri_1, hrvc_lri_2, hrvc_lri_3) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(1,2,3)))))

hrvc_lri.phe.df <- merge(cas.phe.df, hrvc_lri.df, 
												 by="id", all.x=T)

write.table(hrvc_lri.phe.df, "output/pheno/hrvc_lri.phe", row.names=F, quote=F)

# Get totige
totige.df <- master.df %>% select(SUBJ_ID, totige6m, totige1,
																	totige2, totige3, totige4, totige5) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(0.5,1,2,3,4,5)))))

totige.phe.df <- merge(cas.phe.df, totige.df, 
											 by="id", all.x=T)

write.table(totige.phe.df, "output/pheno/totige.phe", row.names=F, quote=F)

# Get hdmige
hdmige.df <- master.df %>% select(SUBJ_ID, hdmige6m, hdmige1,
																	hdmige2, hdmige3, hdmige4, hdmige5) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(0.5,1,2,3,4,5)))))

hdmige.phe.df <- merge(cas.phe.df, hdmige.df, 
												 by="id", all.x=T)

write.table(hdmige.phe.df, "output/pheno/hdmige.phe", row.names=F, quote=F)

# Get pntige
pntige.df <- master.df %>% select(SUBJ_ID, pntige6m, pntige1,
																	pntige2, pntige3, pntige4, pntige5) %>%
	melt(id.vars = "SUBJ_ID",
			 variable.name = "age",
			 value.name = "y") %>%
	rename(id = SUBJ_ID) %>%
	mutate(age = as.numeric(as.character(
		factor(age, labels=c(0.5,1,2,3,4,5)))))

pntige.phe.df <- merge(cas.phe.df, pntige.df, 
											 by="id", all.x=T)

write.table(pntige.phe.df, "output/pheno/pntige.phe", row.names=F, quote=F)


