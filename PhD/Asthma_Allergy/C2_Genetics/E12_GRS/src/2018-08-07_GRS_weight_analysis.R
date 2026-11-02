# =============================================================================
# 2018-08-07 - GRS weight analysis
# =============================================================================

rm(list=ls())

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------

# Bioconductor-specific packages
# ---------------------------------------
desired.biocon.packages <- c("biomformat",
														 "metagenomeSeq",
														 "microbiome")

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
desired.packages <- c("cluster",
											"dplyr", 
											"ecodist",
											"gee",
											"gplots",
											"ggplot2", 
											"GMD", 
											"gtools",
											"meta",
											"RColorBrewer",
											"reshape2",
											"riverplot")

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
# Custom functions
# ------------------------------------------------------------------------------

# A better signif
# ---------------------------------------
good.signif <- function(x, n) {
	trimws(formatC(signif(x,digits=n), digits=n,format="fg"))
}
# Note that this converts number to character string


# Inverse rank-based transform
# ---------------------------------------

iNT <- function(x, lower, upper) {
	c = 3/8  
	qnorm( 
		(rank(x, na.last="keep") - c) / (length(which(!is.na(x))) - 2*c + 1) 
	) 
}


# ------------------------------------------------------------------------------
# Custom settings
# ------------------------------------------------------------------------------

# ggplot and graphical settings
# ---------------------------------------
theme_set(theme_bw())

# Sequential Palette - light tan to dark red
seqPalette <- colorRampPalette(c("#fef0d9","#fdd49e","#fdbb84",
																 "#fc8d59","#e34a33","#b30000"))

# Quantitative Palette - red, blue, green, purple,
# orange,yellow, brown, pink, grey
quaPalette <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
								'#ff7f00','#ffff33','#a65628','#f781bf',
								'#999999')
# light version of above
quaPalette.light <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4',
											'#fed9a6','#ffffcc','#e5d8bd','#fddaec',
											'#f2f2f2')

# for npEM groups
quaPalette.npem <- c("#8da0cb","#a6d854","#fc8d62","#e78ac3","#66c2a5","#ffd92f")

# Divergent Palette - blue to yellow-tan to red
divPalette <- colorRampPalette(rev(c("#d73027","#fc8d59","#fee090",
																		 "#e0f3f8","#91bfdb","#4575b4")))

# Utility function to show colours of hex codes
# grid::grid.raster(seqPalette(6), interpolate = FALSE)
# grid::grid.raster(quaPalette, interpolate = FALSE)
# grid::grid.raster(divPalette(6), interpolate = FALSE)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# GRS data
load("output/2018-08-07_GRS_associations.Rdata")

# Load summary statistics for each GRS
allergic.grs.weights.df <- read.table("data/sumstats/allergic_grs_weights.txt",
																			header=F, stringsAsFactors = F)
asthma.grs.weights.df <- read.table("data/sumstats/asthma_grs_weights.txt",
																			header=F, stringsAsFactors = F)

colnames(allergic.grs.weights.df) <- c("SNP", "Effect_allele_allergic", "Beta_allergic", "Position_allergic", "Reference_allele_allergic")
colnames(asthma.grs.weights.df) <- c("SNP", "Effect_allele_asthma", "Beta_asthma", "Beta_asthma_2", "Position_asthma", "Reference_allele_asthma")

allergic.asthma.grs.weights.df <- merge(allergic.grs.weights.df, asthma.grs.weights.df,
																				by="SNP", all=T)
allergic.asthma.grs.weights.df.backup <- allergic.asthma.grs.weights.df
# allergic.asthma.grs.weights.df <- allergic.asthma.grs.weights.df.backup

dim(allergic.asthma.grs.weights.df) # 1.8 million SNPs
dim(allergic.grs.weights.df) # 1.8 million
dim(asthma.grs.weights.df) # 0.5 million

# For each SNP, with unequal effect alleles, multiply beta_asthmas by -1
problem_snps <- which(allergic.asthma.grs.weights.df$Effect_allele_allergic !=
																		 	allergic.asthma.grs.weights.df$Effect_allele_asthma)
allergic.asthma.grs.weights.df$Beta_asthma[problem_snps] <- 
	-1*allergic.asthma.grs.weights.df$Beta_asthma[problem_snps]
allergic.asthma.grs.weights.df$Beta_asthma_2[problem_snps] <- 
	-1*allergic.asthma.grs.weights.df$Beta_asthma_2[problem_snps]

extra_problem_snps <- with(allergic.asthma.grs.weights.df[problem_snps, ], 
													SNP[Effect_allele_allergic != Reference_allele_asthma])

# Only one snp remains problematic

allergic.asthma.grs.weights.df[allergic.asthma.grs.weights.df$SNP == extra_problem_snps,]
# Remove problematic SNP - also has inconsistent position
allergic.asthma.grs.weights.df <- allergic.asthma.grs.weights.df[-which(allergic.asthma.grs.weights.df$SNP == extra_problem_snps),]

allergic.asthma.grs.weights.df2 <- allergic.asthma.grs.weights.df %>%
	filter(!is.na(Beta_asthma), !is.na(Beta_allergic))

allergic.asthma.grs.weights.df2 %>%
	ggplot(aes(x = Beta_allergic, y = Beta_asthma)) + 
	geom_point()

allergic.asthma.grs.weights.df2 %>%
	ggplot() + geom_histogram(aes(x = Beta_allergic))
allergic.asthma.grs.weights.df2 %>%
	ggplot() + geom_histogram(aes(x = Beta_asthma))

# Find signals that are in the same direction / have similar signal
# Use arctan vs. Use dot product on y=x
# Get inverse rank transform of values, and find SNPs >2sd from mean
# sd ~=2 for all iNT'ed variables

allergic.asthma.grs.weights.df2 <- allergic.asthma.grs.weights.df2 %>% 
	mutate(allergic_iNT = iNT(Beta_allergic),
				 asthma_iNT = iNT(Beta_asthma),
				 magnitude = (allergic_iNT^2 + asthma_iNT^2)^0.5,
				 angle = atan2(y = allergic_iNT, x = asthma_iNT),
				 dotproduct = magnitude * cos(angle - atan2(y=1,x=1)),
				 dotproduct_iNT = iNT(dotproduct))

allergic.asthma.grs.weights.df2 %>%
	ggplot(aes(x = Beta_allergic, y = allergic_iNT)) + geom_point()
allergic.asthma.grs.weights.df2 %>%
	ggplot(aes(x = Beta_asthma, y = asthma_iNT)) + geom_point()

allergic.asthma.grs.weights.df2 %>%
	ggplot() + geom_histogram(aes(x = dotproduct))
allergic.asthma.grs.weights.df2 %>%
	ggplot() + geom_histogram(aes(x = dotproduct_iNT)) 
allergic.asthma.grs.weights.df2 %>%
	ggplot() + geom_histogram(aes(x = allergic_iNT))
allergic.asthma.grs.weights.df2 %>%
	ggplot() + geom_histogram(aes(x = asthma_iNT)) 

allergic.asthma.grs.weights.df2 <- allergic.asthma.grs.weights.df2 %>%
	mutate(sig.dp = as.numeric(dotproduct_iNT >= 2 | dotproduct_iNT <= -2),
				 sig.allergic.grs = as.numeric(allergic_iNT >= 2 | allergic_iNT <= -2),
				 sig.asthma.grs = as.numeric(asthma_iNT >= 2 | asthma_iNT <= -2),
				 sig.both.grs = as.numeric(sig.allergic.grs == 1 & sig.asthma.grs == 1),
				 sig.all = as.numeric(sig.both.grs == 1 & sig.dp == 1),
				 sig.all.dir = as.numeric(sig.all == 1 & 
				 												 	((allergic_iNT > 0 & asthma_iNT > 0) | 
				 												 	 	(allergic_iNT < 0 & asthma_iNT < 0))),
				 sig.both.grs.dir = as.numeric(sig.both.grs == 1 & 
				 																((allergic_iNT > 0 & asthma_iNT > 0) | 
				 																 	(allergic_iNT < 0 & asthma_iNT < 0))))

# allergic.asthma.grs.weights.df2 %>%
# 	ggplot(aes(y = dotproduct, x = Beta_allergic)) +
# 	geom_point(aes(color = factor(sig.dp))) +
# 	scale_color_manual(values=c("black", "red"))
# allergic.asthma.grs.weights.df2 %>%
# 	ggplot(aes(y = dotproduct, x = Beta_allergic)) +
# 	geom_point(aes(color = factor(sig.both.grs))) + 
# 	scale_color_manual(values=c("black", "red"))

# allergic.asthma.grs.weights.df2 %>%
# 	ggplot(aes(x = Beta_allergic, y = Beta_asthma)) + 
# 	geom_point(aes(color = factor(sig.all.dir))) + 
# 	scale_color_manual(values=c("black", "red"))
allergic.asthma.grs.weights.df2 %>%
	ggplot(aes(x = allergic_iNT, y = asthma_iNT)) + 
	geom_point(aes(color = factor(sig.both.grs.dir))) + 
	scale_color_manual(values=c("black", "red"))

allergic.asthma.grs.weights.df2 %>%
	ggplot(aes(x = allergic_iNT, y = asthma_iNT)) + 
	geom_point(aes(color = factor(sig.dp))) + 
	scale_color_manual(values=c("black", "red")) 

allergic.asthma.grs.weights.df2 %>%
	ggplot(aes(x = Beta_allergic, y = Beta_asthma)) + 
	geom_point(aes(color = factor(sig.both.grs.dir))) + 
	scale_color_manual(values=c("black", "red"))

allergic.asthma.grs.weights.df2 %>%
	ggplot(aes(x = Beta_allergic, y = Beta_asthma)) + 
	geom_point(aes(color = factor(sig.dp))) + 
	scale_color_manual(values=c("black", "red"))

# consistent: 
# with(allergic.asthma.grs.weights.df2, 
#		 table(sig.all.dir, sig.both.grs.dir))

# Get concurrent SNPS
# SNP.list <- allergic.asthma.grs.weights.df2$SNP[which(
#  	allergic.asthma.grs.weights.df2$sig.both.grs.dir == 1)]
# length(SNP.list) #2900

SNP.list <- allergic.asthma.grs.weights.df2$SNP[which(
	allergic.asthma.grs.weights.df2$sig.dp == 1)]
length(SNP.list) #22308

# Get allergic and asthma scores
allergic.grs.weights.strong.df <- allergic.grs.weights.df[which(allergic.grs.weights.df$SNP %in% SNP.list),]
asthma.grs.weights.strong.df <- asthma.grs.weights.df[which(asthma.grs.weights.df$SNP %in% SNP.list),]

write.table(allergic.grs.weights.strong.df, "output/allergic_grs_weights_strong.txt", row.names=F, quote=F, col.names = F)
write.table(asthma.grs.weights.strong.df, "output/asthma_grs_weights_strong.txt", row.names=F, quote=F, col.names = F)

save.image("output/2018-08-07_GRS_weight_analysis.Rdata")
