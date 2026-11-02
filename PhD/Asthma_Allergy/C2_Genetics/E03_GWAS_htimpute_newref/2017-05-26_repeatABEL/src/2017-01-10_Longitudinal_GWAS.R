# 2017-01-10 - Longitudinal and repeated-measures GWAS
# ========

# This script produces results from different R packages that perform 
# longitudinal GWASs via generalized estimating equations (GEE) or 
# generalized linear mixed models (GLMM)

# Load libraries
# ------------------------------------------------------------------------------

desired.packages <- c("GenABEL", "RepeatABEL", "qqman")
new.packages <- desired.packages[!(desired.packages %in% 
	installed.packages()[,"Package"])]

if(length(new.packages)>0) install.packages(new.packages, 
	repos='https://cran.ms.unimelb.edu.au/')
for (i in 1:length(desired.packages)){
	library(desired.packages[i], character.only=T)
}

# RepeatABEL (LMM)
# ------------------------------------------------------------------------------

# Load data to genABEL-readable format
# ---------------------------------------
convert.snp.tped("cas_imp_3.tped", "cas_imp_3.tfam", 
	outfile="cas_imp_3.raw", strand="u", bcast=10000)

gen.data <- load.gwaa.data(phenofile = "pheno/cas_imp_3_sex.phe", 
	genofile = "cas_imp_3.raw")

# Test wlri (y1-3) phenotype
# ---------------------------------------
wlri.phen.data <- read.table("pheno/wlri_num_test.phe", header=T)

wlri.gwas <- rGLS(y ~ age + sex, genabel.data = gen.data, 
							 phenotype.data = wlri.phen.data)

wlri.gwas.total.df <- cbind(wlri.gwas@results, wlri.gwas@annotation) 
wlri.gwas.total.df$SNP <- rownames(wlri.gwas.total.df)
colnames(wlri.gwas.total.df)[which(colnames(wlri.gwas.total.df) == "Chromosome")] <- "CHR"
colnames(wlri.gwas.total.df)[which(colnames(wlri.gwas.total.df) == "Position")] <- "BP"
colnames(wlri.gwas.total.df)[which(colnames(wlri.gwas.total.df) == "P1df")] <- "P"
wlri.gwas.total.df$CHR <- as.numeric(as.character(wlri.gwas.total.df$CHR))
wlri.gwas.total.df <- wlri.gwas.total.df[order(wlri.gwas.total.df$P, na.last=T),]

# Write table
write.table(wlri.gwas.total.df, "output/wlri_p_value.txt", row.names=F)

# Plot manhattan plot and qqplot
png("output/wlri_test_manplot.png", width=1440, height=720)
manhattan(wlri.gwas.total.df, ylim=c(0, 11),
	main = "Manhattan Plot - wlri", cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

png("output/wlri_test_qqplot.png", width=720,height=720)
qq(wlri.gwas.total.df$P, ylim=c(0,11),
	main = paste0("QQ plot - wlri"), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

sink("output/wlri_lambda.txt")
estlambda(wlri.gwas.total.df$P,plot=F)
sink()

# LongGWAS (Furlotte) (LMM)
# ------------------------------------------------------------------------------

# Load data from variant-major additive component format
# ---------------------------------------
gen_ad.mat <- as.matrix(read.table("cas_imp_3_ad.mat", header=T, stringsAsFactors=F))
gen_ad.mat <- gen_ad.mat/2

# Load LongGWAS.R
source('src/longGWAS.R')

# get kinship matrix
K <- myFastLMM.kinship(gen_ad.mat)