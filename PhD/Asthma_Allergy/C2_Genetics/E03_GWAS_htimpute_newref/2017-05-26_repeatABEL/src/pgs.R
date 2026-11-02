rm(list=ls())
# install.packages("RCurl",repos='https://cran.ms.unimelb.edu.au/')
# install.packages("httr",repos='https://cran.ms.unimelb.edu.au/')
library(RCurl)
library(httr)
library(devtools)
# set_config( config( ssl_verifypeer = 0L ) )
# install_github("YinanZheng/PGS",
#							 dependencies=TRUE)
library(PGS)
# BJdata()
# 
# # Blood micro RNA (miRNA) expression data of 166 miRNAs and lung function (FEV1) 
# # data of 120 individuals with two measurements in Beijing Air Pollution Study.
# 
# BJlung$gender = factor(BJlung$gender)
# BJlung$heat = factor(BJlung$heat)
# BJlung$cigwear = factor(BJlung$cigwear)
# 
# # Merge and order by subject id (SID) and examination day (WD)
# BJ_dat <- merge(BJmirna, BJlung, by=c("SID","WD"))
# BJ_dat <- BJ_dat[with(BJ_dat, order(SID, WD)), ]
# 
# y.vect <- BJ_dat$FEV1       ### Extract dependent variable
# id.vect <- BJ_dat$SID       ### Extract subject ID variable to indicate repeated measures
# M <- BJ_dat[,3:168]         ### Extract microRNA data matrix (row=samples, column=miRNA probes)
# COV <- BJ_dat[,170:179]     ### Extract covariate data matrix
# 
# colnames(BJ_dat)[1:2] # SID, WD
# colnames(BJ_dat)[3:168] # snps/microRNA
# colnames(BJ_dat)[169] # y
# colnames(BJ_dat)[170:179] # covariates (PCs, age)
# dim(BJ_dat)
# 
# sis_LMM = sis(y.vect, id.vect, M, COV)
# sis_GEE = sis(y.vect, id.vect, M, COV, method = "GEE")
# 
# coef(sis_LMM)
# plot(sis_LMM)
# plot(sis_GEE)
# 
# PGSfit = pgsfit(y.vect, id.vect, M, COV, 
# 								sis.obj = sis_LMM, 
# 								lambda.lim = c(3,5), 
# 								pm.n = 14, 
# 								pm.max = 140, 
# 								seed = 1)
# 
# plot(PGSfit)

# Use LongGWAS test data

source('src/longGWAS.R')

# First load the SNP and phenotype files
snps = as.matrix(read.table('data/snps.tab')) # 100 x 500 matrix of SNP by individual
phenos = as.matrix(read.table('data/phenos.tab')) # 500 x 2 matrix of Individual by time point

snps.t <- t(snps)
snps.df <- data.frame(snps.t)
snps.df$id = seq(1:ncol(snps))
dim(snps.df)

phenos.df <- data.frame(id = seq(1:ncol(snps)),
												sex = rep(1, ncol(snps)),
												phenos)
colnames(phenos.df) <- c("id", "sex", 1,2)
phenos.m.df <- melt(phenos.df, 
										id.vars = c("id", "sex"),
										variable.name = "age",
										value.name = "y")
dim(phenos.m.df)
test.df <- merge(snps.df, phenos.m.df, by="id", all=T)
test.df

start.time.1 <- Sys.time()
# Get PCs
pc <- prcomp(snps, rank=7)
plot(pc$rotation[,1], pc$rotation[,2])
pc.df <- data.frame(pc$rotation)
pc.df$id <- seq(1:ncol(snps))
dim(pc.df)

test.df <- merge(test.df, pc.df, by="id", all=T)
colnames(test.df)

y.vect <- test.df$y       ### Extract dependent variable
id.vect <- test.df$id       ### Extract subject ID variable to indicate repeated measures
M <- test.df[,2:101]        ### Extract microRNA data matrix (row=samples, column=miRNA probes)
COV <- test.df[,c(103,105:111)]     ### Extract covariate data matrix

end.time.1 <- Sys.time()

# LMM
start.time.2 <- Sys.time()
sis_LMM <- sis(y.vect, id.vect, M, COV)
lmm_results.df <- data.frame(name = sis_LMM$name, 
														 p = sis_LMM$p.value)
lmm_results.df <- lmm_results.df[match(colnames(snps.df), lmm_results.df$name),]
lmm_results.df <- lmm_results.df[which(!is.na(lmm_results.df$name)),]
end.time.2 <- Sys.time()

pdf("output/pgs_sis_lmm_man_plot.pdf")
plot(-log10(lmm_results.df$p))													
dev.off()

pdf("output/pgs_sis_lmm_qq_plot.pdf")
qq(lmm_results.df$p)
dev.off()

chisq <- qchisq(1-lmm_results.df$p,1)
lambda <- median(chisq)/qchisq(0.5,1)
write.table(lambda, "output/pgs_sis_lmm_lambda.txt")

total.time <- (end.time.1 - start.time.1) + 
	(end.time.2 - start.time.2)
write.table(total.time, "output/pgs_sis_lmm_compute_time.txt")

results.df <- data.frame(id = seq(1:nrow(lmm_results.df)),
												 p = lmm_results.df$p)
write.table(results.df, "output/pgs_sis_lmm_results.txt",
						row.names=F, quote=F)

# GEE
start.time.3 <- Sys.time()
sis_GEE <- sis(y.vect, id.vect, M, COV, method = "GEE")
gee_results.df <- data.frame(name = sis_GEE$name, 
														 p = sis_GEE$p.value)
gee_results.df <- gee_results.df[match(colnames(snps.df), gee_results.df$name),]
gee_results.df <- gee_results.df[which(!is.na(gee_results.df$name)),]
end.time.3 <- Sys.time()

pdf("output/pgs_sis_gee_man_plot.pdf")
plot(-log10(gee_results.df$p))													
dev.off()

pdf("output/pgs_sis_gee_qq_plot.pdf")
qq(gee_results.df$p)
dev.off()

chisq <- qchisq(1-gee_results.df$p,1)
lambda <- median(chisq)/qchisq(0.5,1)
write.table(lambda, "output/pgs_sis_gee_lambda.txt")


total.time <- (end.time.1 - start.time.1) + 
	(end.time.3 - start.time.3)
write.table(total.time, "output/pgs_sis_gee_compute_time.txt")

results.df <- data.frame(id = seq(1:nrow(gee_results.df)),
												 p = gee_results.df$p)
write.table(results.df, "output/pgs_sis_gee_results.txt",
						row.names=F, quote=F)