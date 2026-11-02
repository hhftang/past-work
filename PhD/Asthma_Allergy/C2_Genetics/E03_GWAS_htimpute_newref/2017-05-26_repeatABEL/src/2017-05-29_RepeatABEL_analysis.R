library(GenABEL)
library(RepeatABEL)
convert.snp.tped("cas_imp_3.tped", "cas_imp_3.tfam", outfile="cas_imp_3.raw", strand="u", bcast=10000)
gen.data <- load.gwaa.data(phenofile = "pheno/cas_imp_3_sex.phe", genofile = "cas_imp_3.raw")
wlri_phen.data <- read.table("pheno/wlri_num_test.phe", header=T)
wlri_gwas <- rGLS(y ~ age + sex, genabel.data = gen.data, 
							 phenotype.data = wlri_phen.data)
# Get top SNPs
sig_snp_index <- which(wlri_gwas@results$P1df < 5*10^(-8))

wlri_gwas_total.df <- cbind(wlri_gwas@results, wlri_gwas@annotation) 
wlri_gwas_total.df$SNP <- rownames(wlri_gwas_total.df)
colnames(wlri_gwas_total.df)[which(colnames(wlri_gwas_total.df) == "Chromosome")] <- "CHR"
colnames(wlri_gwas_total.df)[which(colnames(wlri_gwas_total.df) == "Position")] <- "BP"
colnames(wlri_gwas_total.df)[which(colnames(wlri_gwas_total.df) == "P1df")] <- "P"
wlri_gwas_total.df$CHR <- as.numeric(as.character(wlri_gwas_total.df$CHR))

# Plot manhattan plot and qqplot
library(qqman)

png("output/wlri_test_manplot.png", width=1440, height=720)
manhattan(wlri_gwas_total.df, ylim=c(0, 11),
	main = "Manhattan Plot - wlri", cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

png("output/wlri_test_qqplot.png", width=720,height=720)
qq(wlri_gwas_total.df$P, ylim=c(0,11),
	main = paste0("QQ plot - wlri"), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

sink("output/wlri_lambda.txt")
estlambda(wlri_gwas_total.df$P,plot=F)
sink()

hrva_lri_phen.data <- read.table("pheno/hrva_lri_num_test.phe", header=T)
hrva_lri_gwas <- rGLS(y ~ age + sex, genabel.data = gen.data, 
							 phenotype.data = hrva_lri_phen.data)
# Get top SNPs
sig_snp_index <- which(hrva_lri_gwas@results$P1df < 5*10^(-8))

hrva_lri_gwas_total.df <- cbind(hrva_lri_gwas@results, hrva_lri_gwas@annotation) 
hrva_lri_gwas_total.df$SNP <- rownames(hrva_lri_gwas_total.df)
colnames(hrva_lri_gwas_total.df)[which(colnames(hrva_lri_gwas_total.df) == "Chromosome")] <- "CHR"
colnames(hrva_lri_gwas_total.df)[which(colnames(hrva_lri_gwas_total.df) == "Position")] <- "BP"
colnames(hrva_lri_gwas_total.df)[which(colnames(hrva_lri_gwas_total.df) == "P1df")] <- "P"
hrva_lri_gwas_total.df$CHR <- as.numeric(as.character(hrva_lri_gwas_total.df$CHR))

# Plot manhattan plot and qqplot

png("output/hrva_lri_test_manplot.png", width=1440, height=720)
manhattan(hrva_lri_gwas_total.df, ylim=c(0, 11),
	main = "Manhattan Plot - hrva_lri", cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

png("output/hrva_lri_test_qqplot.png", width=720,height=720)
qq(hrva_lri_gwas_total.df$P, ylim=c(0,11),
	main = paste0("QQ plot - hrva_lri"), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

sink("output/hrva_lri_lambda.txt")
estlambda(hrva_lri_gwas_total.df$P,plot=F)
sink()

# ----

lri_phen.data <- read.table("pheno/lri_num_test.phe", header=T)
lri_gwas <- rGLS(y ~ age + sex, genabel.data = gen.data, 
							 phenotype.data = lri_phen.data)
# Get top SNPs
sig_snp_index <- which(lri_gwas@results$P1df < 5*10^(-8))

lri_gwas_total.df <- cbind(lri_gwas@results, lri_gwas@annotation) 
lri_gwas_total.df$SNP <- rownames(lri_gwas_total.df)
colnames(lri_gwas_total.df)[which(colnames(lri_gwas_total.df) == "Chromosome")] <- "CHR"
colnames(lri_gwas_total.df)[which(colnames(lri_gwas_total.df) == "Position")] <- "BP"
colnames(lri_gwas_total.df)[which(colnames(lri_gwas_total.df) == "P1df")] <- "P"
lri_gwas_total.df$CHR <- as.numeric(as.character(lri_gwas_total.df$CHR))

# Plot manhattan plot and qqplot

png("output/lri_test_manplot.png", width=1440, height=720)
manhattan(lri_gwas_total.df, ylim=c(0, 11),
	main = "Manhattan Plot - lri", cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

png("output/lri_test_qqplot.png", width=720,height=720)
qq(lri_gwas_total.df$P, ylim=c(0,11),
	main = paste0("QQ plot - lri"), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

sink("output/lri_lambda.txt")
estlambda(lri_gwas_total.df$P,plot=F)
sink()

# ----

crwz_phen.data <- read.table("pheno/crwz_bin_test.phe", header=T)
crwz_gwas <- rGLS(y ~ age + sex, genabel.data = gen.data, 
							 phenotype.data = crwz_phen.data)
# Get top SNPs
sig_snp_index <- which(crwz_gwas@results$P1df < 5*10^(-8))

crwz_gwas_total.df <- cbind(crwz_gwas@results, crwz_gwas@annotation) 
crwz_gwas_total.df$SNP <- rownames(crwz_gwas_total.df)
colnames(crwz_gwas_total.df)[which(colnames(crwz_gwas_total.df) == "Chromosome")] <- "CHR"
colnames(crwz_gwas_total.df)[which(colnames(crwz_gwas_total.df) == "Position")] <- "BP"
colnames(crwz_gwas_total.df)[which(colnames(crwz_gwas_total.df) == "P1df")] <- "P"
crwz_gwas_total.df$CHR <- as.numeric(as.character(crwz_gwas_total.df$CHR))

# Plot manhattan plot and qqplot

png("output/crwz_test_manplot.png", width=1440, height=720)
manhattan(crwz_gwas_total.df, ylim=c(0, 11),
	main = "Manhattan Plot - crwz", cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

png("output/crwz_test_qqplot.png", width=720,height=720)
qq(crwz_gwas_total.df$P, ylim=c(0,11),
	main = paste0("QQ plot - crwz"), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

sink("output/crwz_lambda.txt")
estlambda(crwz_gwas_total.df$P,plot=F)
sink()
