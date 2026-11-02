desired.packages <- c("lme4", "ggplot2", "reshape2")
new.packages <- desired.packages[!(desired.packages %in% 
																	 	installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages, 
																						repos='https://cran.ms.unimelb.edu.au/')
for (i in 1:length(desired.packages)){
	library(desired.packages[i], character.only=T)
}

# lme4 test data (sleep deprivation study)

(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
summary(fm1)# (with its own print method; see class?merMod % ./merMod-class.Rd

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

# Get PCs
pc <- prcomp(snps, rank=7)
plot(pc$rotation[,1], pc$rotation[,2])
pc.df <- data.frame(pc$rotation)
pc.df$id <- seq(1:ncol(snps))
dim(pc.df)

test.df <- merge(test.df, pc.df, by="id", all=T)
colnames(test.df)
test.df$id <- factor(test.df$id)

snp.var <- colnames(test.df)[2:101] 
covar.var <- colnames(test.df)[105:107]
f <- as.formula(paste0("y ~ age + ",
						paste(snp.var, collapse=" + "), " + ",
						paste(covar.var, collapse=" + "), " + (1 | id)"))
lmm <- lmer(f, data=test.df)
lmm
