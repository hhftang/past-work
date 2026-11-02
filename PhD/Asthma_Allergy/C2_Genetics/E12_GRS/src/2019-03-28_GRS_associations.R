# =============================================================================
# 2019-03-28 - GRS analysis
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
											"riverplot",
											"psych", 
											"pROC", 
											"PRROC", 
											"NMF", 
											"mediation", 
											"ggfortify")

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

# Helpful package with functions for 
# multiple rarefaction
# OTU abundance averaging following Compositional Data Analysis workflow
# Phylogenetic diversity estimation, diversity profiles (Hill numbers!)
# Pairwise dissimilarity plots 
# Prevalence plots (total OTU abundance vs. OTU prevalence)

if(!("metagMisc" %in% installed.packages()[,"Package"])){
	library(devtools)
	install_github("vmikk/metagMisc")
}
library(metagMisc)

# ------------------------------------------------------------------------------
# Custom functions
# ------------------------------------------------------------------------------

# A better signif
# ---------------------------------------
good.signif <- function(x, n) {
	trimws(formatC(signif(x,digits=n), digits=n,format="fg"))
}
# Note that this converts number to character string

# Shumei's custom functions
# ---------------------------------------
source("src/Functions.R")

# "relative_abundance" = calculates relative abundance for rows = samples or OTUs 
# "multiplot" =  multi-panel plotting
# "ggsurvplot" = KM survival curves
# "resetPar" = resets par to default
# "adj2edge.edit" = converts adjacency matrix to edge list for import into cytoscape
# "plotTimeSeries.new" = re-written version of metagenomeSeq function, plots 
# difference in abundance for particular bacteria
# "GetColorHexAndDecimal" = retrieves hex code from R colour terms
# "put.fig.letter" = inserts figure letters in figures

# Custom colour ramp function for heatmap.3
# ---------------------------------------
# Color.FUN accepts length of vector as argument 
color.custom.function <- function(n) {
	c(col2hex("lightgrey"), 
		colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E4E4CB")))(n-1))
}

# ------------------------------------------------------------------------------
# Custom settings
# ------------------------------------------------------------------------------

# ggplot and graphical settings
# ---------------------------------------
theme_set(theme_bw())

# Sequential Palette - light tan to dark red
# ---------------------------------------
seqPalette <- colorRampPalette(c("#fef0d9","#fdd49e","#fdbb84",
																 "#fc8d59","#e34a33","#b30000"))

# Alternate sequential palette - darker tan to dark red
# ---------------------------------------
seqPalette.alt <- colorRampPalette(c("#fdbb84",
																		 "#fc8d59","#e34a33","#b30000"))

# Quantitative Palette - red, blue, green, purple,
# orange,yellow, brown, pink, grey
# ---------------------------------------
quaPalette <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
								'#ff7f00','#ffff33','#a65628','#f781bf',
								'#999999')
# light version of above
# ---------------------------------------
quaPalette.light <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4',
											'#fed9a6','#ffffcc','#e5d8bd','#fddaec',
											'#f2f2f2')

# for npEM groups
# ---------------------------------------
quaPalette.npem <- c("#8da0cb","#a6d854","#fc8d62","#e78ac3","#66c2a5","#ffd92f")

# Divergent Palette - blue to yellow-tan to red
# ---------------------------------------
divPalette <- colorRampPalette(rev(c("#d73027","#fc8d59","#fee090",
																		 "#e0f3f8","#91bfdb","#4575b4")))
divPalette.alt <- colorRampPalette(rev(c("#d73027","#FFFFFF", "#4575b4")))
divPalette.alt2 <- colorRampPalette(rev(c("#d73027","#4575b4")))


# Utility function to show colours of hex codes
# ---------------------------------------
# grid::grid.raster(seqPalette(6), interpolate = FALSE)
# grid::grid.raster(quaPalette, interpolate = FALSE)
# grid::grid.raster(divPalette(6), interpolate = FALSE)


# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

# npem cluster data
load("data/2017-07-06_npem_global_heatmaps.Rdata")

# Microbiome data
load("data/2018-06-28_MPG_illness_assoc_cas.Rdata")

# ------------------------------------------------------------------------------
# Load GRS files
# ------------------------------------------------------------------------------
# Plink snp scoring file (.profile) format
# CNT = Number of non-missing SNPs used for scoring
# CNT2 = The number of named alleles

# List of phenotypes
# ------------------------------------------------------------------------------

grs_phenotype <- c("allergic",
									 "as",
									 "asthma",
									 "aa",
									 "ca",
									 "ar",
									 "nar",
									 "COPD", 
									 "SpiroMeta_FEV1",
									 "UKB_FEV1",
									 "Meta_FEV1",
									 "SpiroMeta_FVC",
									 "UKB_FVC",
									 "Meta_FVC",
									 "SpiroMeta_RATIO",
									 "UKB_RATIO",
									 "Meta_RATIO",
									 "SpiroMeta_PEF",
									 "UKB_PEF",
									 "Meta_PEF")
grs_phenotype_name <- c("Allergic Disease (Ferreira 2017)",
												"Allergic Sens. (Waage 2018)",
												"Asthma (Demenais 2018)",
												"Adult Asthma (Pividori 2018)",
												"Child Asthma (Pividori 2018)",
												"Allergic Rhinitis (Waage 2018)",
												"N-A Rhinitis (Waage 2018)",
												"COPD (Neale 2018)",
												"FEV1 SM (Shrine 2019)",
												"FEV1 UKB (Shrine 2019, UKB)",
												"FEV1 Meta (Shrine 2019)",
												"FVC SM (Shrine 2019)",
												"FVC UKB (Shrine 2019, UKB)",
												"FVC Meta (Shrine 2019)",
												"FEV1:FVC SM (Shrine 2019)",
												"FEV1:FVC UKB (Shrine 2019)",
												"FEV1:FVC Meta (Shrine 2019)",
												"PEF SM (Shrine 2019)",
												"PEF UKB (Shrine 2019)",
												"PEF Meta (Shrine 2019)")
												
# allergic = Allergic disease - Ferreira et al 2017
	# https://www.nature.com/articles/ng.3985
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics/Allergic_disease/SHARE-without23andMe.LDSCORE-GC.SE-META.v0
	# UKB r2=0.9 LD-thinned set
	# Meta-analysis across 13 datasets (including UK biobank and 23andMe)
	# All of European Ancestry
	# Phenotype = asthma and/or hay fever and/or eczema
# as = Allergic Sensitisation (inhalant allergens) - Waage et al 2018 (European ancestry)
	# https://www.nature.com/articles/s41588-018-0157-1
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics/hayfever/
	# UKB r2=0.9 LD-thinned set
	# (includes ALSPAC, COPSAC, RAINE and others; does not include 23andme)
	# Criterion: IgE >0.7 IU/L for inhalant allergens, or SPT positive
# asthma = Demenais et al 2018 (European subpopulation only)
	# https://www.nature.com/articles/s41588-017-0014-7
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics/Asthma/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv
	# UKB r2=0.9 LD-thinned set
	# Using beta from european only fixed effects model
# aa = Adult Asthma - Pividori et al 2018 (UKB)
	# https://www.biorxiv.org/content/early/2018/09/26/427427
	# UKB r2=0.9 LD-thinned set
	# 21564 cases, 318237 controls
# ca = Childhood Asthma - Pividori et al 2018 (UKB)
	# https://www.biorxiv.org/content/early/2018/09/26/427427
	# UKB r2=0.9 LD-thinned set
	# 9433 cases, 318237 controls
# ar = Allergic Rhinitis - Waage et al 2018 (European ancestry)
	# https://www.nature.com/articles/s41588-018-0157-1
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics/hayfever/
	# UKB r2=0.9 LD-thinned set
	# (includes 23andMe dataset, COPSAC, RAINE and others)
	# NOTE: From methods:
	# "Cases were defined as individuals who had ever received a diagnosis or 
	# experienced symptoms of allergic rhinitis, on the basis of available 
	# phenotype definitions in the included studies (Supplementary Table 3 and 
	# cohort-recruitment details in Supplementary Note). ... (i.e. questionnaire-based)
	# To maximize numbers and optimize statistical power, we did not require 
	# doctor-diagnosed allergic rhinitis or verification by allergic sensitization. 
	# This approach was confirmed by a sensitivity analysis in 23andMe on the basis 
	# of association with known risk loci for allergic rhinitis (data not shown). 
	# Controls were defined as individuals who had never received a diagnosis or 
	# experienced symptoms of allergic rhinitis."
# nar = Non-allergic rhinitis - Waage et al 2018  summary statistics (European ancestry)
	# https://www.nature.com/articles/s41588-018-0157-1
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics/hayfever/
	# UKB r2=0.9 LDthinned set
	# (includes 23andMe dataset, COPSAC, RAINE and others)
	# According to Waage 2018:
	# "Nonallergic rhinitis, defined as having rhinitis symptoms without evidence of 
	# allergic sensitization, is a common but poorly understood disease entity. 
	# We performed what is, to our knowledge, the first GWAS on this phenotype, 
	# with the goal of identifying specific rhinitis mechanisms. 
	# The analysis included 2,028 cases and 9,606 controls from nine studies but 
	# did not identify any risk loci at the genome-wide-significant level. 
	# Comparison with the allergic rhinitis results suggested some overlap in 
	# susceptibility loci (Supplementary Note and Supplementary Table 23)."
# COPD = COPD - UK Biobank 
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/raw_SummaryStatistics/UKB_summary_statistics/
	# 20002_1112.assoc.tsv.gz	
	# from Ben Neale, https://www.dropbox.com/s/hi1zl3asxkty320/20002_1112.assoc.tsv.gz?dl=0 
	# self-reported: chronic obstructive airways disease/copd		
	# UKB r2=0.9 LD-thinned set
# Spirometa/UKB/Meta_FEV1 = Forced expiratory volume in 1 second (Shrine 2019)
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/COPD_spirometa_ukb
	# and https://www.amazon.co.uk/clouddrive/share/slxYFgFImzTyIQWfkqR8EvGRQmlVC1Sso3ABkRiRuPh?_encoding=UTF8&%2AVersion%2A=1&%2Aentries%2A=0&mgh=1
	# paper: https://www.nature.com/articles/s41588-018-0321-7
	# UKB r2=0.9 LD-thinned set
# Spirometa/UKB/Meta_FVC = Forced vital capacity (Shrine 2019)
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/COPD_spirometa_ukb
	# and https://www.amazon.co.uk/clouddrive/share/slxYFgFImzTyIQWfkqR8EvGRQmlVC1Sso3ABkRiRuPh?_encoding=UTF8&%2AVersion%2A=1&%2Aentries%2A=0&mgh=1
	# paper: https://www.nature.com/articles/s41588-018-0321-7
	# UKB r2=0.9 LD-thinned set
# Spirometa/UKB/Meta_RATIO = FEV1/FVC ratio (Shrine 2019)
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/COPD_spirometa_ukb
	# and https://www.amazon.co.uk/clouddrive/share/slxYFgFImzTyIQWfkqR8EvGRQmlVC1Sso3ABkRiRuPh?_encoding=UTF8&%2AVersion%2A=1&%2Aentries%2A=0&mgh=1
	# paper: https://www.nature.com/articles/s41588-018-0321-7
	# UKB r2=0.9 LD-thinned set
# Spirometa/UKB/Meta_PEF = Peak expiratory flow (Shrine 2019)
	# Based on summary statistics in 
	# /projects/GenomicPrediction/GRS_resources/COPD_spirometa_ukb
	# and https://www.amazon.co.uk/clouddrive/share/slxYFgFImzTyIQWfkqR8EvGRQmlVC1Sso3ABkRiRuPh?_encoding=UTF8&%2AVersion%2A=1&%2Aentries%2A=0&mgh=1
	# paper: https://www.nature.com/articles/s41588-018-0321-7
	# UKB r2=0.9 LD-thinned set


# Create dataframe for each phenotype
# ------------------------------------------------------------------------------
master_GRS.df <- master_npem_groups.df

for (i in 1:length(grs_phenotype)) {
	df <- read.table(paste0("data/current/", grs_phenotype[i], ".profile"), header=T, stringsAsFactors=F)
	score_grs <- sym(paste0("SCORESUM_", grs_phenotype[i]))
	cnt_grs <- sym(paste0("CNT_", grs_phenotype[i]))
	cnt2_grs <- sym(paste0("CNT2_", grs_phenotype[i]))
	score_grs_std <- sym(paste0("SCORESUM_",grs_phenotype[i],"_std"))
	score_grs_tert <- sym(paste0("SCORESUM_",grs_phenotype[i],"_tert"))
	df <- df %>% rename(SUBJ_ID = FID,
											!!score_grs := SCORESUM,
											!!cnt_grs := CNT,
											!!cnt2_grs := CNT2) %>%	
		mutate(IID = NULL, 
					 PHENO = NULL,
					 !!score_grs_std := as.numeric(scale(!!score_grs)))
	tert <- with(df, quantile(eval(score_grs_std), 
														probs=seq(0,1,1/3), na.rm=T))
	df <- df %>% mutate(!!score_grs_tert := 
												factor(cut(!!score_grs_std, tert), 
															 labels=c("Low", "Med", "High")))
	assign(paste0(grs_phenotype[i], "_GRS.df"), df)
	master_GRS.df <- merge(master_GRS.df, df,
												 by="SUBJ_ID", all=T)
}

CNT2 <- paste0("CNT2_", grs_phenotype)
grs_phenotype_raw <- paste0("SCORESUM_", grs_phenotype)
grs_phenotype_std <- paste0("SCORESUM_", grs_phenotype, "_std")
grs_phenotype_tert <- paste0("SCORESUM_", grs_phenotype, "_tert")

allergic_tert <- quantile(master_GRS.df$SCORESUM_allergic_std,
													probs = seq(0, 1, 1/3), na.rm=T)
master_GRS.df$SCORESUM_allergic_tert <- 
	factor(cut(master_GRS.df$SCORESUM_allergic_std, allergic_tert), 
				 labels=c("Low", "Med", "High"))

# ------------------------------------------------------------------------------
# Descriptives
# ------------------------------------------------------------------------------

# Number of named alleles used in scoring per subject
# ------------------------------------------------------------------------------

pdf("output/2019-03-28_GRS_scored_alleles.pdf", width=5, height=5)
master_GRS.df %>%
	dplyr::select(SUBJ_ID, one_of(CNT2)) %>%
	melt(id.vars = "SUBJ_ID", 
			 variable.name = "GRS",
			 value.name = "Count") %>%
	mutate(GRS = factor(GRS, labels=grs_phenotype_name)) %>%
	ggplot(aes(x = GRS, y = Count)) + geom_boxplot() + 
	labs(y = "Number of alleles scored per subject") + 
	theme(axis.text.x = element_text(angle = 90, hjust=1))
dev.off()

# lapply(master_GRS.df %>%
# 			 	select(one_of(CNT2)), FUN=function(x) {mean(x, na.rm=T)})

pdf("output/2019-03-28_GRS_scored_alleles_scaled.pdf", width=5, height=5)
master_GRS.df %>%
	dplyr::select(SUBJ_ID, one_of(CNT2)) %>%
	mutate_at(vars(CNT2), funs(scale)) %>%
	melt(id.vars = "SUBJ_ID", 
			 variable.name = "GRS",
			 value.name = "Count") %>%
	mutate(GRS = factor(GRS, labels=grs_phenotype_name)) %>%
	ggplot(aes(x = GRS, y = Count)) + geom_boxplot() + 
	labs(y = "Number of alleles scored per subject \n (scaled to standard normal distribution)") + 
	theme(axis.text.x = element_text(angle = 90, hjust=1))
dev.off()

pdf("output/2019-03-28_GRS_scored_alleles_histogram_count.pdf", width=5, height=5)
master_GRS.df %>%
	dplyr::select(SUBJ_ID, one_of(CNT2)) %>%
	melt(id.vars = "SUBJ_ID", 
			 variable.name = "GRS",
			 value.name = "Count") %>%
	mutate(GRS = factor(GRS, labels=grs_phenotype_name)) %>%
	ggplot(aes(x = Count/1000)) + geom_histogram() + 
	facet_wrap(.~GRS, scales = "free_x") +
	labs(x = "Allele count (10^3)", y = "Frequency of allele count")
dev.off()

# No obvious outliers

# ------------------------------------------------------------------------------
# Get variables (BMI, sensitisation)
# ------------------------------------------------------------------------------

# Change crwz5 to factor for ease of graphing/analysis later
# ------------------------------------------------------------------------------
master_GRS.df$crwz5 <- factor(master_GRS.df$crwz5, levels=c(0,1))

# BMI
# ------------------------------------------------------------------------------
master_GRS.df$bmi3 <- master_GRS.df$wt3/(master_GRS.df$ht3/100)^2
master_GRS.df$bmi4 <- master_GRS.df$wt4/(master_GRS.df$ht4/100)^2
master_GRS.df$bmi5 <- master_GRS.df$wt5/(master_GRS.df$ht5/100)^2
master_GRS.df$bmi10 <- master_GRS.df$wt10/(master_GRS.df$ht10)^2

# Get ae-type variables from premaster.df
# ------------------------------------------------------------------------------
ae.df <- with(premaster.df, data.frame(SUBJ_ID, ae6m, ae1, ae2, ae3, 
																			 ae4, ae5, ae10, 
																			 stringsAsFactors = F))
master_GRS.df <- merge(master_GRS.df, ae.df, by="SUBJ_ID")

# ------------------------------------------------------------------------------
# Histogram of score distributions
# ------------------------------------------------------------------------------

pdf("output/2019-03-28_GRS_histograms.pdf", width=10, height=5)
master_GRS.df %>%
	dplyr::select(SUBJ_ID, one_of(grs_phenotype_raw)) %>%
	melt(id.vars = "SUBJ_ID", 
			 variable.name = "GRS",
			 value.name = "Raw") %>%
	mutate(GRS = factor(GRS, labels=grs_phenotype_name)) %>%
	ggplot(aes(x = Raw)) + geom_histogram(bins=20) + 
	facet_wrap(.~ GRS, scales = "free_x") +
	labs(x = "Raw GRS score", y = "Count") + 
	theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

pdf("output/2019-03-28_GRS_histograms_std.pdf", width=10, height=5)
master_GRS.df %>%
	dplyr::select(SUBJ_ID, one_of(grs_phenotype_std)) %>%
	melt(id.vars = "SUBJ_ID", 
			 variable.name = "GRS",
			 value.name = "Std") %>%
	mutate(GRS = factor(GRS, labels=grs_phenotype_name)) %>%
	ggplot(aes(x = Std)) + geom_histogram(fill="grey", bins=20) + 
	facet_wrap(.~ GRS) +
	labs(x = "Standardised GRS score", y = "Count") + 
	stat_function(color = "red",
								fun = function(x, mean, sd, n, bw) {
									dnorm(x = x, mean = mean, sd = sd)*n*bw
								},
								args = c(mean = 0, sd = 1, n = 215, bw = 5/20))
dev.off()

master_GRS.df %>%
	dplyr::select(SUBJ_ID, one_of(grs_phenotype_std)) %>%	
	melt(id.vars = "SUBJ_ID", 
			 variable.name = "GRS",
			 value.name = "Std") %>%
	mutate(GRS = factor(GRS, labels=grs_phenotype_name)) %>%
	group_by(GRS) %>%
	summarize(shapiro.stat = shapiro.test(Std)[[1]],
						shapiro.p = shapiro.test(Std)[[2]])

# p > 0.05 for all except FEV1 Meta and FVC SpiroMeta

# ------------------------------------------------------------------------------
# PCA of scores
# ------------------------------------------------------------------------------

x <- master_GRS.df[, which(colnames(master_GRS.df) %in% grs_phenotype_std)]
pca <- prcomp(na.omit(x))
str(pca)

# PCA plot
# ------------------------------------------------------------------------------
pca.plot <- as.data.frame(pca$x)
pca.plot <- data.frame(SUBJ_ID = 
											 	master_GRS.df$SUBJ_ID[which(
											 		!is.na(master_GRS.df$SCORESUM_asthma_std))], 
											 group = 
											 	master_GRS.df$group[which(
											 		!is.na(master_GRS.df$SCORESUM_asthma_std))], 											 	
											 pca.plot, 
											 stringsAsFactors = F)
pdf("output/2019-03-28_pca_GRS_npem.pdf", height=5, width=6)
pca.plot %>%
	filter(!is.na(group)) %>%
	mutate(group = factor(group, labels=c("CAS1", "CAS2", "CAS3"))) %>%
	ggplot(aes(x = PC1, y = PC2)) + 
	geom_point(aes(color = group)) + 
	scale_color_manual(values = quaPalette.npem) + 
	labs(color = "npEM cluster")
dev.off()

pca.plot <- pca.plot %>%
	mutate(group.color = factor(group, 
															levels=c(1,2,3),
															labels=quaPalette.npem[1:3]))
with(pca.plot, 
		 scatterplot3d::scatterplot3d(x = PC1, y = PC2, z = PC3,
		 														 color = group.color))

# with(pca.plot, 
#  		 rgl::plot3d(x = PC1, y = PC2, z = PC3,
#  		 						col = group.color))

# Loadings
# ------------------------------------------------------------------------------
loadings <- as.data.frame(pca$rotation)
loadings$GRS <- rownames(loadings)

grs_phenotype_colors <- c("#e41a1c", "#377eb8", "#87b1d4", "#214b6e", 
																 "#4daf4a", "#984ea3", "#cba6d1", "#ff7f00",
																 "#b25800", "#b25800", "#b25800", 
																 "#ffa54c", "#ffa54c", "#ffa54c", 
																 "#663200", "#663200", "#663200",
																 "#ffcb99", "#ffcb99", "#ffcb99")

# Colour loadings by phenotype
pdf("output/2019-03-28_pca_GRS_loadings.pdf", height=5, width=16)
loadings %>%
	melt(id.vars = "GRS",
			 variable.name = "PC",
			 value.name = "loading") %>%
	mutate(GRS = factor(GRS, levels=grs_phenotype_std,
											labels=grs_phenotype_name)) %>%
	filter(as.numeric(PC) <= 14) %>%
	ggplot(aes(x = PC, y = loading, group=GRS, fill=GRS)) +
	geom_bar(stat = "identity", position="dodge") + 
	labs(y = "Loading") +
	# geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5)) + 
	scale_fill_manual(values=grs_phenotype_colors)
dev.off()

# loadings plot
pdf("output/2019-03-28_pca_GRS_loadings_v2.pdf", height=5, width=10)
pca.plot %>%
	ggplot(aes(x = PC1, y = PC2)) + 
	geom_point(alpha=0.5, color="grey") + 
	geom_segment(data = loadings, 
							 aes(x = 0, y = 0, xend = (PC1*8), yend=(PC2*8), 
							 		color = factor(GRS, levels=grs_phenotype_std)),
							 arrow = arrow(length=unit(1/2, "picas"))) +
	annotate("text", x = (loadings$PC1*9), y = (loadings$PC2*9), 
					 label=letters[1:length(grs_phenotype_std)]) +
	scale_color_manual(values=grs_phenotype_colors,
										 labels=paste0(letters, " - ", grs_phenotype_name)) + 
	labs(color = "GRS")
dev.off()

# Proportion of variance explained
# ------------------------------------------------------------------------------
pca.importance <- as.data.frame(summary(pca)$importance)
df <- as.data.frame(c(0, t(pca.importance[3,])))
colnames(df) <- "cumvar"
df$PC <- seq(0, (nrow(df)-1))

pdf("output/2019-03-28_pca_GRS_var.pdf", height=5, width=8)
ggplot(df, aes(x = PC, y = cumvar)) + 
	geom_point() + geom_line() + 
	labs(x = "Number of PCs", 
			 y = "Cumulative proportion of variance explained")
dev.off()

# Get PCs and incorporate in table
# ------------------------------------------------------------------------------
master_GRS.df <- merge(master_GRS.df, 
											 with(pca.plot,
											 		 data.frame(SUBJ_ID, PC1, PC2, PC3, PC4, PC5, PC6, PC7,
											 		 					 PC8, PC9, PC10, PC11, PC12, PC13, PC14)), 
											 by="SUBJ_ID", all=T)

# ------------------------------------------------------------------------------
# Global score (not useful at this stage)
# ------------------------------------------------------------------------------
# master_GRS.df$SCORESUM_combined_std <- 
# 	with(master_GRS.df, SCORESUM_allergic_std + SCORESUM_asthma_std + 
# 			 	SCORESUM_adult_asthma_std + SCORESUM_childhood_asthma_std + 
# 			 	SCORESUM_ar_std + SCORESUM_as_std + SCORESUM_nar_std)
# # The sum of normally distributed random variables is normally distributed, with mean = sum(mean),
# # and var = sum(var)
# # https://en.wikipedia.org/wiki/Sum_of_normally_distributed_random_variables
# 

# ------------------------------------------------------------------------------
# Get tertiles for PC1
# ------------------------------------------------------------------------------
master_GRS.df$SCORESUM_PC1_std <- 
 	with(master_GRS.df, as.numeric(scale(master_GRS.df$PC1)))

PC1_tert <- quantile(master_GRS.df$SCORESUM_PC1_std,
										 probs = seq(0, 1, 1/3), na.rm=T)
master_GRS.df$SCORESUM_PC1_tert <- factor(cut(master_GRS.df$SCORESUM_PC1_std, PC1_tert), 
																					labels=c("Low", "Med", "High"))

# ------------------------------------------------------------------------------
# Correlation of scores
# ------------------------------------------------------------------------------
x <- master_GRS.df[, which(colnames(master_GRS.df) %in% grs_phenotype_std)]
cor(x, method="pearson", use="pairwise.complete.obs")
y <- corr.test(x, method="pearson", use="pairwise.complete.obs")

z.r <- y$r
z.r[upper.tri(z.r, diag=T)] <- NA
z.r <- z.r %>%
	melt(value.name="corr.r")

z.p <- y$p
z.p[upper.tri(z.p, diag=T)] <- NA
z.p <- z.p %>%
	melt(value.name="corr.p")

z <- merge(z.r, z.p, by=c("Var1", "Var2"))
z$corr.p.sig <- factor(ifelse(z$corr.p < 0.05, "Sig", "NS"))

pdf("output/2019-03-28_correlation_GRS.pdf", width=20, height=16)
z %>% 
	mutate(Var1 = factor(Var1, labels=grs_phenotype_name),
				 Var2 = factor(Var2, labels=grs_phenotype_name)) %>%
	ggplot(aes(x = Var1, y = Var2)) + 
	geom_tile(aes(fill = corr.r)) +
	geom_text(aes(label = paste0("R = ", good.signif(corr.r,2), 
															 "\n p = ", signif(corr.p, 2)),
								color = corr.p.sig))+
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1,1)) +
	scale_color_manual(breaks=c("NS", "Sig"), values = c("grey", "black")) + 
	labs(x = "GRS", y="GRS", fill="Rho", color="Sig.") +
	theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				axis.text.x = element_text(angle = 45, hjust = 1)) +
	geom_hline(yintercept = c(2.5, 5.5, 7.5, 8.5, 11.5, 14.5, 17.5, 20.5))+
	geom_vline(xintercept = c(2.5, 5.5, 7.5, 8.5, 11.5, 14.5, 17.5, 20.5))
dev.off()

# ------------------------------------------------------------------------------
# Get lung function data.frames
# ------------------------------------------------------------------------------

# Yr 5 lung function
# ------------------------------------------------------------------------------
lungfunction.df <- read.csv("data/CAS yr5 resp data.csv", header=T, stringsAsFactors = F)
lungfunction.df <- lungfunction.df %>%
	rename(Respiratory..postratio = Respiratory..posratio) %>%
	mutate(Respiratory..testdate = as.Date(strptime(Respiratory..testdate, "%m/%d/%Y")),
				 Respiratory..prefvc = as.numeric(as.character(Respiratory..prefvc)),
				 Respiratory..postfvc = as.numeric(as.character(Respiratory..postfvc)),
				 Respiratory..chgfvc_old = as.numeric(as.character(Respiratory..chgfvc))/100,
				 Respiratory..chgfvc = (Respiratory..postfvc - Respiratory..prefvc)/
				 	Respiratory..prefvc,
				 Respiratory..prefev1 = as.numeric(as.character(Respiratory..prefev1)),
				 Respiratory..postfev1 = as.numeric(as.character(Respiratory..postfev1)),
				 Respiratory..chgfev1_old = as.numeric(as.character(Respiratory..chgfev1))/100,
				 Respiratory..chgfev1 = (Respiratory..postfev1 - Respiratory..prefev1)/
				 	Respiratory..prefev1,
				 Respiratory..preratio = as.numeric(as.character(Respiratory..preratio)),
				 Respiratory..postratio = as.numeric(as.character(Respiratory..postratio)),
				 Respiratory..chgratio_old = as.numeric(as.character(Respiratory..chgratio))/100,
				 Respiratory..chgratio = (Respiratory..postratio - Respiratory..preratio)/
				 	Respiratory..preratio,
				 Respiratory..prefef = as.numeric(as.character(Respiratory..prefef)),
				 Respiratory..postfef = as.numeric(as.character(Respiratory..postfef)),
				 Respiratory..chgfef_old = as.numeric(as.character(Respiratory..chgfef))/100,
				 Respiratory..chgfef = (Respiratory..postfef - Respiratory..prefef)/
					Respiratory..prefef,
				 Respiratory..prepef = as.numeric(as.character(Respiratory..prepef)),
				 Respiratory..postpef = as.numeric(as.character(Respiratory..postpef)),
				 Respiratory..chgpef_old = as.numeric(as.character(Respiratory..chgpef))/100,
				 Respiratory..chgpef = (Respiratory..postpef - Respiratory..prepef)/
				 	Respiratory..prepef,
				 FOT..birthdate = as.Date(strptime(FOT..birthdate, "%m/%d/%Y")),
				 FOT..meas.date = as.Date(strptime(FOT..meas.date, "%m/%d/%Y")))

lungfunction.df <- lungfunction.df[which(lungfunction.df$Respiratory..agemonth <= 75),] 

# Yr 10 lung function
# ------------------------------------------------------------------------------

lungfunction_y10.df <- read.csv("data/10yr LFT and SPT.csv", header=T, stringsAsFactors = F)
lungfunction_y10.df$SUBJ_ID_part <- formatC(lungfunction_y10.df$idno, width=4, format="d", flag="0")
lungfunction_y10.df <- lungfunction_y10.df %>%
	mutate(FirstName = NULL, 
				 TypeofMedication = NULL,
				 Frequencyofmedication = NULL,
				 CommAsthma = NULL, 
				 Frequencyofmedication2 = NULL,
				 CommHF = NULL,
				 CommeNO = NULL,
				 TRIGGER = NULL,
				 Hayfever = as.numeric(as.character(
				 	factor(Hayfever, levels=c("No","Yes"), labels=c(0,1)))))

# Age check
lungfunction_y10.df$date_new <- as.Date(strptime("1900-01-01", "%Y-%m-%d")) +
	lungfunction_y10.df$Date - 2
lungfunction_y10.df$DOB_new <- as.Date(strptime("1900-01-01", "%Y-%m-%d")) +
	lungfunction_y10.df$DOB - 2

# need to correct entry 143 in lungfunction_y10.df
# get premaster.df dobs
df <- with(premaster.df, data.frame(SUBJ_ID, dob))
df$SUBJ_ID_part <- substring(df$SUBJ_ID, 0, 4)
df$SUBJ_ID <- NULL
lungfunction_y10.df <- merge(lungfunction_y10.df, df, 
														 by = "SUBJ_ID_part", all=T)
# check
with(lungfunction_y10.df, data.frame(DOB_new, dob, test = DOB_new - dob))

lungfunction_y10.df$age_days <- as.numeric(
	as.character(lungfunction_y10.df$date_new - lungfunction_y10.df$dob))

table(floor(lungfunction_y10.df$age_days/365 * 12))
hist(lungfunction_y10.df$age_days/365)
# most within 6 months of 10 years

# Simplify dataframe terminology for ease of coding
df <- lungfunction_y10.df

# FEV1
# -----
# Find out what best value actually means
df$FEV1best_test <- with(df, pmax(FEV11, FEV12, na.rm=T))
df %>% ggplot(aes(x=FEV1best, y=FEV1best_test)) + geom_point()
# max of two trials

# Find out what FEV1pred means
# FEV1 predicted for age/height/weight/sex?
df %>% ggplot(aes(x=FEV1best, y=FEV1pred)) + geom_point()

# Find out what FEV180 actually means 
# FEV1 80% in dictionary, but value given as numeric, not boolean
# Therefore it does not represent true/false for FEV1<>80%
df %>% ggplot(aes(x=FEV1best, y=FEV180)) + geom_point() # better aligned
df %>% ggplot(aes(x=FEV1pred, y=FEV180)) + geom_point() 
df %>% ggplot(aes(x=FEV1best*0.8, y=FEV180)) + geom_point() # better aligned, roughly linear
df %>% ggplot(aes(x=FEV1pred*0.8, y=FEV180)) + geom_point() 
# rougly equal to FEV1best * 0.8 - why this is needed is not clear (makes more sense if FEV1pred*0.8)

# Find out what the A suffix means
# According to dictionary, FEV1pred_A = FEV1pred%, % of predicted
df %>% ggplot(aes(x=FEV1best, y=FEV1pred_A)) + geom_point()
df %>% ggplot(aes(x=FEV1pred, y=FEV1pred_A)) + geom_point()
df %>% ggplot(aes(x=FEV1best/FEV1pred, y=FEV1pred_A)) + geom_point()
# FEV1pred_A = FEV1best/FEV1pred

# FVC
# -----
df$FVCbest_test <- with(df, pmax(FVC1, FVC2, na.rm=T))
df %>% ggplot(aes(x=FVCbest, y=FVCbest_test)) + geom_point()
# usually max of two trials, but there are some that are not
df[which(df$FVCbest != df$FVCbest_test),]
length(which(df$FVCbest != df$FVCbest_test))
hist(with(df[which(df$FVCbest != df$FVCbest_test),], abs(FVC1 - FVC2)))
# those that are discordant appear to be within 0.5 L of each other (max - min < 0.5)

df[which(df$FVCbest == df$FVCbest_test),]
length(which(df$FVCbest == df$FVCbest_test))
hist(with(df[which(df$FVCbest == df$FVCbest_test),], abs(FVC1 - FVC2)))
# those that are consistent tend to (but not always) have wide differences


df[which(df$FVCbest - df$FVCbest_test < -1),]
x <- which(df$FVCbest - df$FVCbest_test < -1)
# severe outlier (index 235), out by > 1, due to typo
df[x,]$FVCbest <- df[x,]$FVC2

df$FVC_reproducibility <- 
	with(df, (pmax(FVC1, FVC2, na.rm=T) - pmin(FVC1, FVC2, na.rm=T))/
			 	pmax(FVC1, FVC2, na.rm=T))
hist(df$FVC_reproducibility)
df %>% ggplot(aes(x=FVC_reproducibility, y=FVCbest_test - FVCbest)) + geom_point()
# Ideal reproducibility (technically, relative difference rather than similarity) < 0.05 to 0.10
# Does not appear to be any association between reproducibility and choice of FVCbest
# Obviously, if they choose the "wrong" one, the more different (less reproducible), the greater the discordance

# Unlike Year 5 data, we are also given the original raw measurements (for Year 5, 
# we assume that they automatically  chose the best measurement, as suggested
# by FOT..bestFVC and similar variables - although again these are not always
# concordant.)

# decision: replace FVCbest with FVCbest_test
df$FVCbest <- df$FVCbest_test

# FVCpred
df %>% ggplot(aes(x=FVCbest_test, y=FVCpred)) + geom_point()
df %>% ggplot(aes(x=FVCpred, y=FVCpred_A)) + geom_point()
lungfunction_y10.df %>% ggplot(aes(x=FVCbest/FVCpred, y=FVCpred_A)) + geom_point()
df %>% ggplot(aes(x=FVCbest_test/FVCpred, y=FVCpred_A)) + geom_point()
# FVCpred_A = FVCbest/FVCpred

# decision: use new FVCpred_A_test based on FVCbest_test
df$FVCpred_A_test <- with(df, FVCbest_test/FVCpred)
df %>% ggplot(aes(x=FVCbest_test/FVCpred, y=FVCpred_A_test)) + geom_point()
df$FVCpred_A <- df$FVCpred_A_test

# FEF2575
# -----
df$FEF2575best_test <- with(df, pmax(FEF25751, FEF25752, na.rm=T))
df %>% ggplot(aes(x=FEF2575best, y=FEF2575best_test)) + geom_point()
# need to deal with two major outliers
df[which(df$FEF2575best_test > 100),]
# typos - one is clearly wrong, and one is out by a factor of 100
# correct typos
df[which(df$FEF25751 == 999),]$FEF25751 <- NA
df[which(df$FEF25752 == 190),]$FEF25752 <- 1.90

df$FEF2575best_test <- with(df, pmax(FEF25751, FEF25752, na.rm=T))
df %>% ggplot(aes(x=FEF2575best, y=FEF2575best_test)) + geom_point()
# some discrepancies remain
# usually max of two trials, but there are some that are not
df[which(df$FEF2575best != df$FEF2575best_test),]
length(which(df$FEF2575best != df$FEF2575best_test))
hist(with(df[which(df$FEF2575best != df$FEF2575best_test),], abs(FEF25751 - FEF25752)))
# again, not clear why. Those that are discordant tend to be within 0.6 L of each other
df[which(df$FEF2575best == df$FEF2575best_test),]
length(which(df$FEF2575best == df$FEF2575best_test))
hist(with(df[which(df$FEF2575best == df$FEF2575best_test),], abs(FEF25751 - FEF25752)))
# again, those that are consistent tend to (but not always) have wide differences

# decision: replace FEF2575best with FEF2575best_test
df$FEF2575best <- df$FEF2575best_test

# FEF2575pred
df %>% ggplot(aes(x=FEF2575best_test, y=FEF2575pred)) + geom_point()
df %>% ggplot(aes(x=FEF2575pred, y=FEF2575pred_A)) + geom_point()
lungfunction_y10.df %>% ggplot(aes(x=FEF2575best/FEF2575pred, y=FEF2575pred_A)) + geom_point()
df %>% ggplot(aes(x=FEF2575best_test/FEF2575pred, y=FEF2575pred_A)) + geom_point()
# FEF2575pred_A = FEF2575best/FEF2575pred

# decision: use new FEF2575pred_A_test based on FEF2575best_test
df$FEF2575pred_A_test <- with(df, FEF2575best_test/FEF2575pred)
df %>% ggplot(aes(x=FEF2575best_test/FEF2575pred, y=FEF2575pred_A_test)) + geom_point()
df$FEF2575pred_A <- df$FEF2575pred_A_test

# Post-bronchodilator FEV1
# -----
df$FEV1Bbest_test <- with(df, pmax(FEV1B1, FEV1B2, na.rm=T))
df %>% ggplot(aes(x=FEV1Bbest, y=FEV1Bbest_test)) + geom_point()
df[which(df$FEV1Bbest != df$FEV1Bbest_test),]
# Seems like someone had accidentally rounded off FEV1B2 ...
df[which(df$FEV1Bbest == df$FEV1Bbest_test),]
table(df[which(df$FEV1Bbest == df$FEV1Bbest_test),]$FEV1B2)
# All consistent values had FEV1B2 == 0 for some reason.

# Decision - ** replace affected/discordant FEV1B2 with FEV1Bbest **
df[which(df$FEV1Bbest != df$FEV1Bbest_test),]$FEV1B2 <- 
	df[which(df$FEV1Bbest != df$FEV1Bbest_test),]$FEV1Bbest
df$FEV1Bbest_test <- with(df, pmax(FEV1B1, FEV1B2, na.rm=T))
df$FEV1Bbest <- df$FEV1Bbest_test

# FEV1Bpred and FEV1B  = FEV1Bpred_A? (create new variable)
df %>% ggplot(aes(x=FEV1Bbest_test, y=FEV1Bpred)) + geom_point()
df$FEV1Bpred_A <- with(df, FEV1Bbest_test/FEV1Bpred)
df %>% ggplot(aes(x=FEV1B, y=FEV1Bpred_A)) + geom_point()
# Consistent
# But to maintain consistent terminology - rename/create new variable FEV1Bpred_A

# Post-bronchodilator FVC
# -----
df$FVCBbest_test <- with(df, pmax(FVCB1, FVCB2, na.rm=T))
df %>% ggplot(aes(x=FVCBbest, y=FVCBbest_test)) + geom_point()
df[which(df$FVCBbest != df$FVCBbest_test),]
# Again, same error - someone had accidentally rounded off FVCB2 ...
df[which(df$FVCBbest == df$FVCBbest_test),]
table(df[which(df$FVCBbest == df$FVCBbest_test),]$FVCB2)
# Most consistent values had FVCB2 == 0 for some reason.

# Decision - ** replace affected/discordant FVCB2 with FVCBbest **
df[which(df$FVCBbest != df$FVCBbest_test),]$FVCB2 <- 
	df[which(df$FVCBbest != df$FVCBbest_test),]$FVCBbest
df$FVCBbest_test <- with(df, pmax(FVCB1, FVCB2, na.rm=T))
df$FVCBbest <- df$FVCBbest_test

# FVCBpred and FVCBpred_A (create new variable)
df %>% ggplot(aes(x=FVCBbest_test, y=FVCBpred)) + geom_point()
df$FVCBpred_A <- with(df, FVCBbest_test/FVCBpred)
df %>% ggplot(aes(x=FVCB, y=FVCBpred_A)) + geom_point()

# very strong outlier [index 10] - to consider excluding?
df[which(df$FVCB > 200),]
df[which(df$FVCB > 200),]$FVCB <- NA
df[which(df$FVCBpred_A > 2),]$FVCBpred_A <- NA
df[which(df$FVCBbest > 5),]$FVCBbest <- NA
df[which(df$FVCBbest_test > 5),]$FVCBbest_test <- NA

# Post-bronchodilator FEF2575
# -----
df$FEF2575Bbest_test <- with(df, pmax(FEF2575B1, FEF2575B2, na.rm=T))
df %>% ggplot(aes(x=FEF2575Bbest, y=FEF2575Bbest_test)) + geom_point()
df[which(df$FEF2575Bbest != df$FEF2575Bbest_test),]
# Again, same error - someone had accidentally rounded off FEF2575B2 ...
df[which(df$FEF2575Bbest == df$FEF2575Bbest_test),]
table(df[which(df$FEF2575Bbest == df$FEF2575Bbest_test),]$FEF2575B2)
# Almost all consistent values had FEF2575B2 == 0 

# Decision - ** replace affected/discordant FEF2575B2 with FEF2575Bbest **
df[which(df$FEF2575Bbest != df$FEF2575Bbest_test),]$FEF2575B2 <- 
	df[which(df$FEF2575Bbest != df$FEF2575Bbest_test),]$FEF2575Bbest
df$FEF2575Bbest_test <- with(df, pmax(FEF2575B1, FEF2575B2, na.rm=T))
df$FEF2575Bbest <- df$FEF2575Bbest_test

# FEF2575Bpred and FEF2575Bpred_A (create new variable)
df %>% ggplot(aes(x=FEF2575Bbest_test, y=FEF2575Bpred)) + geom_point()
df$FEF2575Bpred_A <- with(df, FEF2575Bbest_test/FEF2575Bpred)
df %>% ggplot(aes(x=FEF2575B, y=FEF2575Bpred_A)) + geom_point()
# Some discrepancies

with(df[which(abs(df$FEF2575B/100 - df$FEF2575Bpred_A) > 0.1),],
		 data.frame(FEF2575Bbest, FEF2575B/100, FEF2575Bpred_A))
with(lungfunction_y10.df[which(abs(df$FEF2575B/100 - df$FEF2575Bpred_A) > 0.1),],
		 data.frame(FEF2575Bbest, FEF2575B/100))
# Discrepancies to do with incorrect (original) FEF2575Bbest and hence FEF2575B...

# Decision - calculate FEF2575Bpred_A accordingly and replace/ignore FEF2575B
df$FEF2575Bpred_A <- with(df, FEF2575Bbest/FEF2575Bpred)
df$FEF2575B <- df$FEF2575Bpred_A*100
hist(df$FEF2575B - lungfunction_y10.df$FEF2575B)

# Sex, gender, height and weight
with(df, table(Sex, gender, useNA="ifany"))
# two individuals of female sex, unspecified gender
# standardise sex notation

df$sex <- as.numeric(if_elsese(df$Sex == "f", 0,
														ifelse(df$Sex == "F", 0,
																	 ifelse(df$Sex == "M", 1, NA))))
with(df, table(sex, Sex, useNA="ifany"))														
df$Sex <- NULL

df$ht10 <- df$Heightcms/100
df$wt10 <- df$Weight

# Overwrite original lungfunction_y10.df with df
lungfunction_y10.df_old <- lungfunction_y10.df
lungfunction_y10.df <- df

df2 <- master_GRS.df %>%
	dplyr::select(SUBJ_ID)
df2$SUBJ_ID_part <- substring(df2$SUBJ_ID, 0, 4)
lungfunction_y10.df <- merge(df2, lungfunction_y10.df, by="SUBJ_ID_part", all=T)

# ------------------------------------------------------------------------------
# Consolidate
# ------------------------------------------------------------------------------

# Define transient and persistent wheeze
# ------------------------------------------------------------------------------
master_GRS.df$transwz <- 
	with(master_GRS.df, as.numeric((crwz1 + crwz2 + crwz3 > 0) & crwz5 == 0))
master_GRS.df$perswz <- 
	with(master_GRS.df, as.numeric((crwz1 + crwz2 + crwz3 > 0) & crwz5 == 1))

# Merge with trajectory data.frame
# ------------------------------------------------------------------------------
traj.df <- read.csv("data/cas.trajectory.csv", header=T)
traj.df <- traj.df %>%
	rename(traj = traj.timepoint)
master_GRS.df <- merge(master_GRS.df, traj.df, by="SUBJ_ID", all.x=T)


# Merge with illness mpg data.frame
# ------------------------------------------------------------------------------
illness.mpg.df <- read.csv("data/cas.illness.mpg.csv", header=T)
illness.mpg.df <- illness.mpg.df  %>%
	rename(illness.mpg.healthy.age2.prop = Prop,
				 illness.mpg.healthy.age2.bin = Prop.bin)
master_GRS.df <- merge(master_GRS.df, illness.mpg.df, by="SUBJ_ID", all.x=T)

# Merge with lung function data.frames
# ------------------------------------------------------------------------------
lungfunction.y5.df <- lungfunction.df %>%
	dplyr::mutate(FVC.pre.y5 = Respiratory..prefvc,
								FVC.post.y5 = Respiratory..postfvc,
								FVC.chg.y5 = Respiratory..chgfvc,
								FEV1.pre.y5 = Respiratory..prefev1,
								FEV1.post.y5 = Respiratory..postfev1,
								FEV1.chg.y5 = Respiratory..chgfev1,
								FEV1.FVC.pre.y5 = Respiratory..preratio/100,
								FEV1.FVC.post.y5 = Respiratory..postratio/100,
								FEV1.FVC.chg.y5 = Respiratory..chgfvc,
								FEF.pre.y5 = Respiratory..prefef,
								FEF.post.y5 = Respiratory..postfef,
								FEF.chg.y5 = Respiratory..chgfef,
								PEF.pre.y5 = Respiratory..prepef,
								PEF.post.y5 = Respiratory..postpef,
								PEF.chg.y5 = Respiratory..chgpef) %>%
	dplyr::select(SUBJ_ID, FVC.pre.y5, FVC.post.y5, FVC.chg.y5,
								FEV1.pre.y5, FEV1.post.y5, FEV1.chg.y5,
								FEV1.FVC.pre.y5, FEV1.FVC.post.y5, FEV1.FVC.chg.y5,
								FEF.pre.y5, FEF.post.y5, FEF.chg.y5,
								PEF.pre.y5, PEF.post.y5, PEF.chg.y5)
master_GRS.df <- merge(master_GRS.df, lungfunction.y5.df, by="SUBJ_ID", all.x=T)

lungfunction.y10.df <- lungfunction_y10.df %>%
	dplyr::mutate(FVC.pre.y10 = FVCbest,
								FVC.post.y10 = FVCBbest,
								FVC.chg.y10 = (FVCBbest - FVCbest)/FEV1best,
								FEV1.pre.y10 = FEV1best,
								FEV1.post.y10 = FEV1Bbest,
								FEV1.chg.y10 = (FEV1Bbest - FEV1best)/FEV1best,
								FEV1.FVC.pre.y10 = FEV1best/FVCbest,
								FEV1.FVC.post.y10 = FEV1Bbest/FVCBbest,
								FEV1.FVC.chg.y10 = (FEV1.FVC.post.y10 - FEV1.FVC.pre.y10)/FEV1.FVC.pre.y10,
								FEF2575.pre.y10 = FEF2575best,
								FEF2575.post.y10 = FEF2575Bbest,
								FEF2575.chg.y10 = (FEF2575.post.y10  - FEF2575.pre.y10)/FEF2575.pre.y10) %>%
	dplyr::select(SUBJ_ID, FVC.pre.y10, FVC.post.y10, FVC.chg.y10,
								FEV1.pre.y10, FEV1.post.y10, FEV1.chg.y10,
								FEV1.FVC.pre.y10, FEV1.FVC.post.y10, FEV1.FVC.chg.y10,
								FEF2575.pre.y10, FEF2575.post.y10, FEF2575.chg.y10)

master_GRS.df <- merge(master_GRS.df, lungfunction.y10.df, by="SUBJ_ID", all.x=T)

# ------------------------------------------------------------------------------
# Violin plots of GRS with certain traits
# ------------------------------------------------------------------------------

# Define new vectors for labels and variables
# ------------------------------------------------------------------------------
# grs_phenotype_tert <- c("SCORESUM_allergic_tert", "SCORESUM_asthma_tert",
# 									"SCORESUM_adult_asthma_tert", "SCORESUM_childhood_asthma_tert",
# 									"SCORESUM_ar_tert", "SCORESUM_as_tert", "SCORESUM_nar_tert",
# 									"SCORESUM_copd_tert",
# 									"SCORESUM_combined_tert", "SCORESUM_PC1_tert")
# tertile_labs <- c("Allergic disease GRS tertile", "Asthma GRS tertile",
# 									"Adult asthma GRS tertile", "Childhood asthma GRS tertile",
# 									"Allergic rhinitis GRS tertile", "Allergic sensitisation GRS tertile",
# 									"Non-allergic rhinitis GRS tertile", "COPD GRS tertile",
# 									"Combined GRS tertile",
# 									"PC1 of all GRS tertile")
# grs_phenotype_std <- c("SCORESUM_allergic_std", "SCORESUM_asthma_std",
# 								"SCORESUM_adult_asthma_std", "SCORESUM_childhood_asthma_std",
# 								"SCORESUM_ar_std", "SCORESUM_as_std", "SCORESUM_nar_std",
# 								"SCORESUM_copd_std",
# 								"SCORESUM_combined_std", "SCORESUM_PC1_std")
# score_labs <- c("Allergic disease GRS", "Asthma GRS",
# 								"Adult asthma GRS", "Childhood asthma GRS",
# 								"Allergic rhinitis GRS", "Allergic sens. GRS",
# 								"Non-allergic rhinitis GRS", "COPD GRS",
# 								"Combined GRS",
# 								"PC1 of all GRS")
grs_phenotype_std <- c(grs_phenotype_std, "SCORESUM_PC1_std")
grs_phenotype_tert <- c(grs_phenotype_tert, "SCORESUM_PC1_tert")
grs_phenotype_name <- c(grs_phenotype_name, "PC1 of all GRS")

outcome_list <- c("crwz5", "crwz10", 
									"cras_5", "cras_10", 
									"arc5", "transwz", "perswz",
									"totige2 >= 2", 
									"totige5 >= 2",
									"hdmige2 >= log10(0.35)",
									"hdmige5 >= log10(0.35)",
									"pntige2 >= log10(0.35)", 
									"ae2", "ae5",
									"lri_1 + lri_2 >= 1", "wlri_1 + wlri_2 >= 1",
									"illness.mpg.healthy.age2.bin == 'High'",
									"group == 3", "group %in% c(2,3)", 
									"traj == 'Traj.A'", "traj == 'Traj.B'", 
									"traj == 'Traj.C1'", "traj == 'Traj.C2'", 
									"traj %in% c('Traj.C1', 'Traj.C2')",
									"smoke2","smoke5")
outcome_name <- c("Wheeze at age 5",
									"Wheeze at age 10",
									"Asthma at age 5",
									"Asthma at age 10",
									"Allergic rhinoconjunctivitis at age 5",
									"Transient wheeze",
									"Persistent wheeze",
									"Total IgE at age 2 >= 100 kU/L",
									"Total IgE at age 5 >= 100 kU/L",
									"HDM IgE at age 2 >= 0.35 kU/L",
									"HDM IgE at age 5 >= 0.35 kU/L",
									"Peanut IgE at age 2 >= 0.35 kU/L",
									"Allergic sensitisation at age 2",
									"Allergic sensitisation at age 5",
									"Any LRI by age 2", "Any wLRI by age 2", 
									"High illness MPG in healthy samples by age 2",
									"npEM cluster 3", "npEM cluster 2/3",
									"Traj.A", "Traj.B", "Traj.C1", "Traj.C2", "Traj.C1/2",
									"Smoke exposure age 2", "Smoke exposure age 10")
outcome_list_cont <- c("totige2", "totige5",
									"hdmige2", "hdmige5", "pntige2", 
									"lri_1 + lri_2", "wlri_1 + wlri_2",
									"illness.mpg.healthy.age2.prop",
									"FEV1.pre.y5", "FEV1.post.y5", "FEV1.chg.y5",
									"FEV1.pre.y10", "FEV1.post.y10", "FEV1.chg.y10",
									"FVC.pre.y5", "FVC.post.y5", "FVC.chg.y5",
									"FVC.pre.y10", "FVC.post.y10", "FVC.chg.y10",
									"FEV1.FVC.pre.y5", "FEV1.FVC.post.y5", "FEV1.FVC.chg.y5",
									"FEV1.FVC.pre.y10", "FEV1.FVC.post.y10", "FEV1.FVC.chg.y10",
									"FEF.pre.y5", "FEF.post.y5", "FEF.chg.y5",
									"FEF2575.pre.y10", "FEF2575.post.y10", "FEF2575.chg.y10",
									"PEF.pre.y5", "PEF.post.y5", "PEF.chg.y5")
outcome_name_cont <- c("Total IgE at age 2",
									"Total IgE at age 5",
									"HDM IgE at age 2",
									"HDM IgE at age 5",
									"Peanut IgE at age 2",
									"Total LRI by age 2",
									"Total wLRI by age 2",
									"Proportion illness MPG in healthy samples by age 2",
									"FEV1 pre-bronch. age 5", "FEV1 post-bronch. age 5", "FEV1 change age 5",
									"FEV1 pre-bronch. age 10", "FEV1 post-bronch. age 10", "FEV1 change age 10",
									"FVC pre-bronch. age 5", "FVC post-bronch. age 5", "FVC change age 5",
									"FVC pre-bronch. age 10", "FVC post-bronch. age 10", "FVC change age 10",
									"FEV1:FVC pre-bronch. age 5", "FEV1:FVC post-bronch. age 5", "FEV1:FVC change age 5",
									"FEV1:FVC pre-bronch. age 10", "FEV1:FVC post-bronch. age 10", "FEV1:FVC change age 10",
									"FEF pre-bronch. age 5", "FEF post-bronch. age 5", "FEF change age 5",
									"FEF2575 pre-bronch. age 10", "FEF2575 post-bronch. age 10", "FEF2575 change age 10",
									"PEF pre-bronch. age 5", "PEF post-bronch. age 5", "PEF change age 5")

# Plot GRS vs. wheeze at age 5
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(crwz5), !is.na(score)) %>%
	mutate(crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(crwz5)) %>%
	group_by(crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(crwz5)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_crwz5.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol=3) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Plot GRS vs. cluster
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, group, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "group"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(group), !is.na(score)) %>%
	mutate(group = factor(group, labels=c("CAS1", "CAS2", "CAS3")))
df2 <- df %>% 
	filter(!is.na(group)) %>%
	group_by(group, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(group)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ group)$statistic,
						kruskal.df = kruskal.test(score ~ group)$parameter,
						kruskal.p = kruskal.test(score ~ group)$p.value,
						kruskal.p.cas3 = kruskal.test(score ~ (group == "CAS3"))$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p.cas3, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_npem.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol = 3) +
	geom_violin(aes(color = group, x = group, y = score)) +
	geom_point(data = df2, aes(color = group, x = group, y = mean)) + 
	geom_errorbar(data = df2, aes(color = group, x = group, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = group, 
														x = group, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("CAS3 p = ", signif(kruskal.p.cas3,2)), 
														x = 2, y = -3.5)) +
	scale_color_manual(name = "npEM cluster", values=quaPalette.npem) +
	labs(x = "npEM cluster", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Interaction with wheeze
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, group, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "group", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(group), !is.na(crwz5), !is.na(score)) %>%
	mutate(group = factor(group, labels=c("CAS1", "CAS2", "CAS3")),
				 crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(group), !is.na(crwz5)) %>%
	group_by(group, crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(group)) %>%
	group_by(GRS, group) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_npem_by_crwz5.pdf", height=40, width=10)
df %>% ggplot() + 
	facet_grid(factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name) ~ group) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Plot GRS vs. ae2
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, ae2, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "ae2"),
			 variable.name = "GRS",
			 value.name = "score") %>%
	filter(!is.na(ae2), !is.na(score)) %>%
	mutate(ae2 = factor(ae2, labels=c("No", 
																		"Yes")))
df2 <- df %>% 
	filter(!is.na(ae2)) %>%
	group_by(ae2, GRS) %>%
	summarize(n = n(),
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(ae2)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ ae2)$statistic,
						kruskal.df = kruskal.test(score ~ ae2)$parameter,
						kruskal.p = kruskal.test(score ~ ae2)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))


pdf("output/2019-03-28_GRS_vs_ae2.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol=3) +
	geom_violin(aes(color = ae2, x = ae2, y = score)) +
	geom_point(data = df2, aes(color = ae2, x = ae2, y = mean)) + 
	geom_errorbar(data = df2, aes(color = ae2, x = ae2, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = ae2, 
														x = ae2, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Sensitisation by age 2", values=c("black", "red")) +
	labs(x = "Sensitisation by age 2", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Interaction with wheeze
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, ae2, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "ae2", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(ae2), !is.na(crwz5), !is.na(score)) %>%
	mutate(ae2 = factor(ae2, labels=c("No Sensitisation", "Sensitisation")),
				 crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(ae2), !is.na(crwz5)) %>%
	group_by(ae2, crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(ae2)) %>%
	group_by(GRS, ae2) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_ae2_by_crwz5.pdf", height=40, width=10)
df %>% ggplot() + 
	facet_grid(factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name) ~ ae2) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Plot GRS vs. smoke2
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, smoke2, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "smoke2"),
			 variable.name = "GRS",
			 value.name = "score") %>%
	filter(!is.na(smoke2), !is.na(score)) %>%
	mutate(smoke2 = factor(smoke2, labels=c("No", 
																		"Yes")))
df2 <- df %>% 
	filter(!is.na(smoke2)) %>%
	group_by(smoke2, GRS) %>%
	summarize(n = n(),
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(smoke2)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ smoke2)$statistic,
						kruskal.df = kruskal.test(score ~ smoke2)$parameter,
						kruskal.p = kruskal.test(score ~ smoke2)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))


pdf("output/2019-03-28_GRS_vs_smoke2.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol=3) +
	geom_violin(aes(color = smoke2, x = smoke2, y = score)) +
	geom_point(data = df2, aes(color = smoke2, x = smoke2, y = mean)) + 
	geom_errorbar(data = df2, aes(color = smoke2, x = smoke2, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = smoke2, 
														x = smoke2, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Smoke exposure at age 2", values=c("black", "red")) +
	labs(x = "Smoke exposure at age 2", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Interaction with wheeze
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, smoke2, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "smoke2", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(smoke2), !is.na(crwz5), !is.na(score)) %>%
	mutate(smoke2 = factor(smoke2, labels=c("No smoke exposure", "Smoke exposure")),
				 crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(smoke2), !is.na(crwz5)) %>%
	group_by(smoke2, crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(smoke2)) %>%
	group_by(GRS, smoke2) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_smoke2_by_crwz5.pdf", height=40, width=10)
df %>% ggplot() + 
	facet_grid(factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name) ~ smoke2) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Plot GRS vs. crcatexp2
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, crcatexp2, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "crcatexp2"),
			 variable.name = "GRS",
			 value.name = "score") %>%
	filter(!is.na(crcatexp2), !is.na(score)) %>%
	mutate(crcatexp2 = factor(crcatexp2, labels=c("No", 
																		"Yes")))
df2 <- df %>% 
	filter(!is.na(crcatexp2)) %>%
	group_by(crcatexp2, GRS) %>%
	summarize(n = n(),
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(crcatexp2)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crcatexp2)$statistic,
						kruskal.df = kruskal.test(score ~ crcatexp2)$parameter,
						kruskal.p = kruskal.test(score ~ crcatexp2)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))


pdf("output/2019-03-28_GRS_vs_crcatexp2.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol=3) +
	geom_violin(aes(color = crcatexp2, x = crcatexp2, y = score)) +
	geom_point(data = df2, aes(color = crcatexp2, x = crcatexp2, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crcatexp2, x = crcatexp2, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crcatexp2, 
														x = crcatexp2, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Cat exposure at age 2", values=c("black", "red")) +
	labs(x = "Cat exposure at age 2", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Interaction with wheeze
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, crcatexp2, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "crcatexp2", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(crcatexp2), !is.na(crwz5), !is.na(score)) %>%
	mutate(crcatexp2 = factor(crcatexp2, labels=c("No cat exposure", "Cat exposure")),
				 crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(crcatexp2), !is.na(crwz5)) %>%
	group_by(crcatexp2, crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(crcatexp2)) %>%
	group_by(GRS, crcatexp2) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_crcatexp2_by_crwz5.pdf", height=40, width=10)
df %>% ggplot() + 
	facet_grid(factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name) ~ crcatexp2) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Plot GRS vs. microbial traj
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, traj, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "traj"),
			 variable.name = "GRS",
			 value.name = "score") %>%
	filter(!is.na(traj), !is.na(score)) %>%
	mutate(traj = factor(traj, labels=c("Traj.A", "Traj.B", "Traj.C1", "Traj.C2")))
df2 <- df %>% 
	filter(!is.na(traj)) %>%
	group_by(traj, GRS) %>%
	summarize(n = n(),
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(traj)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ traj)$statistic,
						kruskal.df = kruskal.test(score ~ traj)$parameter,
						kruskal.p = kruskal.test(score ~ traj)$p.value,
						kruskal.p.trajb = kruskal.test(score ~ as.numeric(traj %in% c("Traj.B")))$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p.trajb, method = "fdr"))


pdf("output/2019-03-28_GRS_vs_traj.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol=3) +
	geom_violin(aes(color = traj, x = traj, y = score)) +
	geom_point(data = df2, aes(color = traj, x = traj, y = mean)) + 
	geom_errorbar(data = df2, aes(color = traj, x = traj, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = traj, 
														x = traj, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("Traj.B p = ", signif(kruskal.p.trajb,2)), 
														x = 2.5, y = -3.5)) +
	scale_color_manual(name = "Microbial trajectory", 
										 values=c("#7FC97F", "#E41A1C", "#BEAED4","#857994")) +
	labs(x = "Microbial trajectory", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Interaction with wheeze
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, traj, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "traj", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(traj), !is.na(crwz5), !is.na(score)) %>%
	mutate(traj = factor(traj, labels=c("Traj.A", "Traj.B", "Traj.C1", "Traj.C2")),
				 crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(traj), !is.na(crwz5)) %>%
	group_by(traj, crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(traj)) %>%
	group_by(GRS, traj) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_traj_by_crwz5.pdf", height=40, width=10)
df %>% ggplot() + 
	facet_grid(factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name) ~ traj) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Plot GRS vs. traj.alt
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, traj, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "traj"),
			 variable.name = "GRS",
			 value.name = "score") %>%
	filter(!is.na(traj), !is.na(score)) %>%
	mutate(traj = factor(traj, levels=c("Traj.A", "Traj.B", "Traj.C1", "Traj.C2"),
											 labels=c("Traj.A", "Traj.B", "Traj.C", "Traj.C")))
df2 <- df %>% 
	filter(!is.na(traj)) %>%
	group_by(traj, GRS) %>%
	summarize(n = n(),
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(traj)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ traj)$statistic,
						kruskal.df = kruskal.test(score ~ traj)$parameter,
						kruskal.p = kruskal.test(score ~ traj)$p.value,
						kruskal.p.trajb = kruskal.test(score ~ as.numeric(traj %in% c("Traj.B")))$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p.trajb, method = "fdr"))


pdf("output/2019-03-28_GRS_vs_traj.alt.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol=3) +
	geom_violin(aes(color = traj, x = traj, y = score)) +
	geom_point(data = df2, aes(color = traj, x = traj, y = mean)) + 
	geom_errorbar(data = df2, aes(color = traj, x = traj, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = traj, 
														x = traj, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("Traj.B p = ", signif(kruskal.p.trajb,2)), 
														x = 2, y = -3.5)) +
	scale_color_manual(name = "Microbial trajectory", 
										 values=c("#7FC97F", "#E41A1C", "#BEAED4")) +
	labs(x = "Microbial trajectory", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Interaction with wheeze
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, traj, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "traj", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(traj), !is.na(crwz5), !is.na(score)) %>%
	mutate(traj = factor(traj, levels=c("Traj.A", "Traj.B", "Traj.C1", "Traj.C2"),
											 labels=c("Traj.A", "Traj.B", "Traj.C", "Traj.C")),
				 crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(traj), !is.na(crwz5)) %>%
	group_by(traj, crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(traj)) %>%
	group_by(GRS, traj) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_traj.alt_by_crwz5.pdf", height=40, width=10)
df %>% ggplot() + 
	facet_grid(factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name) ~ traj) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Plot GRS vs. illness.mpg
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, illness.mpg.healthy.age2.bin, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "illness.mpg.healthy.age2.bin"),
			 variable.name = "GRS",
			 value.name = "score") %>%
	filter(!is.na(illness.mpg.healthy.age2.bin), !is.na(score)) %>%
	mutate(illness.mpg.healthy.age2.bin = factor(illness.mpg.healthy.age2.bin, labels=c("Low", 
																		"High")))
df2 <- df %>% 
	filter(!is.na(illness.mpg.healthy.age2.bin)) %>%
	group_by(illness.mpg.healthy.age2.bin, GRS) %>%
	summarize(n = n(),
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(illness.mpg.healthy.age2.bin)) %>%
	group_by(GRS) %>%
	summarize(kruskal.stat = kruskal.test(score ~ illness.mpg.healthy.age2.bin)$statistic,
						kruskal.df = kruskal.test(score ~ illness.mpg.healthy.age2.bin)$parameter,
						kruskal.p = kruskal.test(score ~ illness.mpg.healthy.age2.bin)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))


pdf("output/2019-03-28_GRS_vs_illness.mpg.healthy.age2.bin.pdf", height=16, width=10)
df %>% ggplot() + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name),
						 ncol=3) +
	geom_violin(aes(color = illness.mpg.healthy.age2.bin, x = illness.mpg.healthy.age2.bin, y = score)) +
	geom_point(data = df2, aes(color = illness.mpg.healthy.age2.bin, x = illness.mpg.healthy.age2.bin, y = mean)) + 
	geom_errorbar(data = df2, aes(color = illness.mpg.healthy.age2.bin, x = illness.mpg.healthy.age2.bin, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = illness.mpg.healthy.age2.bin, 
														x = illness.mpg.healthy.age2.bin, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Sensitisation by age 2", values=c("black", "red")) +
	labs(x = "Sensitisation by age 2", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()

# Interaction with wheeze
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, illness.mpg.healthy.age2.bin, crwz5, one_of(grs_phenotype_std)) %>%
	melt(id.vars = c("SUBJ_ID", "illness.mpg.healthy.age2.bin", "crwz5"),
			 variable.name = "GRS",
			 value.name = "score") %>% 
	filter(!is.na(illness.mpg.healthy.age2.bin), !is.na(crwz5), !is.na(score)) %>%
	mutate(illness.mpg.healthy.age2.bin = factor(illness.mpg.healthy.age2.bin, labels=c("Low", "High")),
				 crwz5 = factor(crwz5, labels=c("No wheeze", "Wheeze")))
df2 <- df %>% 
	filter(!is.na(illness.mpg.healthy.age2.bin), !is.na(crwz5)) %>%
	group_by(illness.mpg.healthy.age2.bin, crwz5, GRS) %>%
	summarize(n = n(), 
						mean = mean(score, na.rm=T),
						sd = sd(score, na.rm=T),
						lci = mean - 1.96*sd/sqrt(n()),
						uci = mean + 1.96*sd/sqrt(n())) %>%
	ungroup()
df3 <- df %>%
	filter(!is.na(illness.mpg.healthy.age2.bin)) %>%
	group_by(GRS, illness.mpg.healthy.age2.bin) %>%
	summarize(kruskal.stat = kruskal.test(score ~ crwz5)$statistic,
						kruskal.df = kruskal.test(score ~ crwz5)$parameter,
						kruskal.p = kruskal.test(score ~ crwz5)$p.value) %>%
	ungroup() %>%
	mutate(kruskal.p.fdr = p.adjust(kruskal.p, method = "fdr"))

pdf("output/2019-03-28_GRS_vs_illness.mpg.healthy.age2.bin_by_crwz5.pdf", height=40, width=10)
df %>% ggplot() + 
	facet_grid(factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name) ~ illness.mpg.healthy.age2.bin) +
	geom_violin(aes(color = crwz5, x = crwz5, y = score)) +
	geom_point(data = df2, aes(color = crwz5, x = crwz5, y = mean)) + 
	geom_errorbar(data = df2, aes(color = crwz5, x = crwz5, 
																ymin = lci, ymax = uci), width = 0.1) +
	geom_text(data = df2, aes(label = n, color = crwz5, 
														x = crwz5, y = 3.5)) +
	geom_text(data = df3, aes(label = paste0("p = ", signif(kruskal.p,2)), 
														x = 1.5, y = -3.5)) +
	scale_color_manual(name = "Age-five wheeze", values=c("Black", "Red")) +
	labs(x = "Age-five wheeze", y = "Score") +
	theme(axis.text.x=element_blank())
dev.off()


# Linear regression function
# ------------------------------------------------------------------------------
regression <- function(df){
	#setting the regression function. 
	reg_fun <- lm(formula=df$len~df$dose) #regression function
	#getting the slope, intercept, R square and adjusted R squared
	slope <- signif(coef(reg_fun)[2],2)  
	intercept <- signif(coef(reg_fun)[1],2) 
	R2 <- signif(as.numeric(summary(reg_fun)[8]),2)
	R2.Adj <- signif(as.numeric(summary(reg_fun)[9]),2)
	p <- signif(as.numeric(summary(reg_fun)$coefficients[2,4]),2)
	c(slope,intercept,R2,R2.Adj,p)
}

# Plot GRS vs. FEV1
# ------------------------------------------------------------------------------
grs_phenotype_std_1 <- sym(grs_phenotype_std[1])

df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, FEV1.pre.y5, one_of(grs_phenotype_std)) %>%
	filter(!is.na(FEV1.pre.y5), !is.na(!!grs_phenotype_std_1)) %>%
	melt(id.vars = c("SUBJ_ID", "FEV1.pre.y5"),
			 variable.name = "GRS",
			 value.name = "score")
df2 <- df %>% 
	group_by(GRS) %>%
	summarize(n = n()) %>%
	ungroup()
n_value <- mean(df2$n, na.rm=T)

regressions_data <- df %>%
	dplyr::rename(len = FEV1.pre.y5,
								dose = score) %>%
	plyr::ddply("GRS",regression)
colnames(regressions_data)<-c ("GRS","slope","intercept","R2","R2.Adj","p.val")

pdf("output/2019-03-28_GRS_vs_fev1.pre.y5.pdf", height=20, width=20)
df %>% ggplot(aes(y = FEV1.pre.y5, x = score)) + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name)) +
	geom_point() + 
	stat_smooth(method = "lm", col = "red") +
	geom_text(data = regressions_data,
						aes(label = paste0("B = ", slope,"; P = ", p.val)),
						x = 0, y = 0.50)
dev.off()

df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, FEV1.pre.y10, one_of(grs_phenotype_std)) %>%
	filter(!is.na(FEV1.pre.y10), !is.na(!!grs_phenotype_std_1)) %>%
	melt(id.vars = c("SUBJ_ID", "FEV1.pre.y10"),
			 variable.name = "GRS",
			 value.name = "score")
df2 <- df %>% 
	group_by(GRS) %>%
	summarize(n = n()) %>%
	ungroup()
n_value <- mean(df2$n, na.rm=T)

regressions_data <- df %>%
	dplyr::rename(len = FEV1.pre.y10,
								dose = score) %>%
	plyr::ddply("GRS",regression)
colnames(regressions_data)<-c ("GRS","slope","intercept","R2","R2.Adj","p.val")

pdf("output/2019-03-28_GRS_vs_fev1.pre.y10.pdf", height=20, width=20)
df %>% ggplot(aes(y = FEV1.pre.y10, x = score)) + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name)) +
	geom_point() + 
	stat_smooth(method = "lm", col = "red") +
	geom_text(data = regressions_data,
						aes(label = paste0("B = ", slope,"; P = ", p.val)),
						x = 0, y = 1.45)
dev.off()

# Plot GRS vs. FVC
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, FVC.pre.y5, one_of(grs_phenotype_std)) %>%
	filter(!is.na(FVC.pre.y5), !is.na(!!grs_phenotype_std_1)) %>%
	melt(id.vars = c("SUBJ_ID", "FVC.pre.y5"),
			 variable.name = "GRS",
			 value.name = "score")
df2 <- df %>% 
	group_by(GRS) %>%
	summarize(n = n()) %>%
	ungroup()
n_value <- mean(df2$n, na.rm=T)

regressions_data <- df %>%
	dplyr::rename(len = FVC.pre.y5,
								dose = score) %>%
	plyr::ddply("GRS",regression)
colnames(regressions_data)<-c ("GRS","slope","intercept","R2","R2.Adj","p.val")

pdf("output/2019-03-28_GRS_vs_FVC.pre.y5.pdf", height=20, width=20)
df %>% ggplot(aes(y = FVC.pre.y5, x = score)) + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name)) +
	geom_point() + 
	stat_smooth(method = "lm", col = "red") +
	geom_text(data = regressions_data,
						aes(label = paste0("B = ", slope,"; P = ", p.val)),
						x = 0, y = 0.50)
dev.off()

df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, FVC.pre.y10, one_of(grs_phenotype_std)) %>%
	filter(!is.na(FVC.pre.y10), !is.na(!!grs_phenotype_std_1)) %>%
	melt(id.vars = c("SUBJ_ID", "FVC.pre.y10"),
			 variable.name = "GRS",
			 value.name = "score")
df2 <- df %>% 
	group_by(GRS) %>%
	summarize(n = n()) %>%
	ungroup()
n_value <- mean(df2$n, na.rm=T)

regressions_data <- df %>%
	dplyr::rename(len = FVC.pre.y10,
								dose = score) %>%
	plyr::ddply("GRS",regression)
colnames(regressions_data)<-c ("GRS","slope","intercept","R2","R2.Adj","p.val")

pdf("output/2019-03-28_GRS_vs_FVC.pre.y10.pdf", height=20, width=20)
df %>% ggplot(aes(y = FVC.pre.y10, x = score)) + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name)) +
	geom_point() + 
	stat_smooth(method = "lm", col = "red") +
	geom_text(data = regressions_data,
						aes(label = paste0("B = ", slope,"; P = ", p.val)),
						x = 0, y = 1.45)
dev.off()

# Plot GRS vs. FEV1.FVC
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, FEV1.FVC.pre.y5, one_of(grs_phenotype_std)) %>%
	filter(!is.na(FEV1.FVC.pre.y5), !is.na(!!grs_phenotype_std_1)) %>%
	melt(id.vars = c("SUBJ_ID", "FEV1.FVC.pre.y5"),
			 variable.name = "GRS",
			 value.name = "score")
df2 <- df %>% 
	group_by(GRS) %>%
	summarize(n = n()) %>%
	ungroup()
n_value <- mean(df2$n, na.rm=T)

regressions_data <- df %>%
	dplyr::rename(len = FEV1.FVC.pre.y5,
								dose = score) %>%
	plyr::ddply("GRS",regression)
colnames(regressions_data)<-c ("GRS","slope","intercept","R2","R2.Adj","p.val")

pdf("output/2019-03-28_GRS_vs_FEV1.FVC.pre.y5.pdf", height=20, width=20)
df %>% ggplot(aes(y = FEV1.FVC.pre.y5, x = score)) + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name)) +
	geom_point() + 
	stat_smooth(method = "lm", col = "red") +
	geom_text(data = regressions_data,
						aes(label = paste0("B = ", slope,"; P = ", p.val)),
						x = 0, y = 0.65)
dev.off()

df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, FEV1.FVC.pre.y10, one_of(grs_phenotype_std)) %>%
	filter(!is.na(FEV1.FVC.pre.y10), !is.na(!!grs_phenotype_std_1)) %>%
	melt(id.vars = c("SUBJ_ID", "FEV1.FVC.pre.y10"),
			 variable.name = "GRS",
			 value.name = "score")
df2 <- df %>% 
	group_by(GRS) %>%
	summarize(n = n()) %>%
	ungroup()
n_value <- mean(df2$n, na.rm=T)

regressions_data <- df %>%
	dplyr::rename(len = FEV1.FVC.pre.y10,
								dose = score) %>%
	plyr::ddply("GRS",regression)
colnames(regressions_data)<-c ("GRS","slope","intercept","R2","R2.Adj","p.val")

pdf("output/2019-03-28_GRS_vs_FEV1.FVC.pre.y10.pdf", height=20, width=20)
df %>% ggplot(aes(y = FEV1.FVC.pre.y10, x = score)) + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name)) +
	geom_point() + 
	stat_smooth(method = "lm", col = "red") +
	geom_text(data = regressions_data,
						aes(label = paste0("B = ", slope,"; P = ", p.val)),
						x = 0, y = 0.675)
dev.off()

# Plot GRS vs. PEF
# ------------------------------------------------------------------------------
df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, PEF.pre.y5, one_of(grs_phenotype_std)) %>%
	filter(!is.na(PEF.pre.y5), !is.na(!!grs_phenotype_std_1)) %>%
	melt(id.vars = c("SUBJ_ID", "PEF.pre.y5"),
			 variable.name = "GRS",
			 value.name = "score")
df2 <- df %>% 
	group_by(GRS) %>%
	summarize(n = n()) %>%
	ungroup()
n_value <- mean(df2$n, na.rm=T)

regressions_data <- df %>%
	dplyr::rename(len = PEF.pre.y5,
								dose = score) %>%
	plyr::ddply("GRS",regression)
colnames(regressions_data)<-c ("GRS","slope","intercept","R2","R2.Adj","p.val")

pdf("output/2019-03-28_GRS_vs_PEF.pre.y5.pdf", height=20, width=20)
df %>% ggplot(aes(y = PEF.pre.y5, x = score)) + 
	facet_wrap(.~factor(GRS, levels=grs_phenotype_std, labels=grs_phenotype_name)) +
	geom_point() + 
	stat_smooth(method = "lm", col = "red") +
	geom_text(data = regressions_data,
						aes(label = paste0("B = ", slope,"; P = ", p.val)),
						x = 0, y = 0.75)
dev.off()


# ------------------------------------------------------------------------------
# GLM associations with traits at specific timepoints
# ------------------------------------------------------------------------------

# Association analyses with GRS
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 grs_phenotype_name = rep(grs_phenotype_name, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs.df <- results.df
results_glm_grs.df$grs_phenotype_name <- factor(results_glm_grs.df$grs_phenotype_name, 
											 levels=unique(results_glm_grs.df$grs_phenotype_name))
results_glm_grs.df$outcome_name <- factor(results_glm_grs.df$outcome_name, 
													 levels=unique(results_glm_grs.df$outcome_name))

write.table(results.df, "output/2019-03-28_glm_grs_results.txt",
						quote=F, sep = "\t", row.names=F)

pdf("output/2019-03-28_glm_GRS.pdf", height=20, width=20)
results_glm_grs.df %>%
	mutate(outcome_name = factor(outcome_name, levels=rev(levels(outcome_name)))) %>%
	ggplot(aes(x = grs_phenotype_name, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", good.signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Early-life outcome in CAS", fill="Log10(OR)", color="p < 0.05?") + 
	geom_hline(yintercept = c(2.5, 7.5, 9.5, 10.5, 12.5, 19.5, 21.5, 22.5)) +
	geom_vline(xintercept = c(2.5, 5.5, 7.5, 8.5, 11.5, 14.5, 17.5, 20.5))
dev.off()

results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 grs_phenotype_name = rep(grs_phenotype_name, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = master_GRS.df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs.df <- results.df
results_glm_cont_grs.df$grs_phenotype_name <- factor(results_glm_cont_grs.df$grs_phenotype_name, 
																				levels=unique(results_glm_cont_grs.df$grs_phenotype_name))
results_glm_cont_grs.df$outcome_name <- factor(results_glm_cont_grs.df$outcome_name, 
																					levels=unique(results_glm_cont_grs.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_results.txt",
						quote=F, sep = "\t", row.names=F)

pdf("output/2019-03-28_glm_cont_GRS.pdf", height=20, width=20)
results_glm_cont_grs.df %>%
	mutate(outcome_name = factor(outcome_name, levels=rev(levels(outcome_name)))) %>%
	ggplot(aes(x = grs_phenotype_name, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, b.exact))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), lim=c(-1,1)) + 
	geom_text(aes(label = paste0("B = ", good.signif(b.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="beta", color="p < 0.05?") +
	geom_vline(xintercept = c(2.5, 5.5, 7.5, 8.5, 11.5, 14.5, 17.5, 20.5))
dev.off()

# Association analyses with GRS - with npEM cluster membership as covariate
# ------------------------------------------------------------------------------
outcome_list.alt <- outcome_list[!outcome_list %in% c("group == 3", "group %in% c(2,3)")]
outcome_name.alt <- outcome_name[!outcome_list %in% c("group == 3", "group %in% c(2,3)")]

results.df <- data.frame(outcome = rep(outcome_list.alt, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list.alt)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list.alt)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list.alt)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list.alt)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list.alt)),
												 outcome_name = rep(outcome_name.alt, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list.alt)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list.alt)),
												 grs_phenotype_name = rep(grs_phenotype_name, times=length(outcome_list.alt)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list.alt)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list.alt)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + I(group == 3) + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_npem3.df <- results.df
results_glm_grs_npem3.df$grs_phenotype_name <- factor(results_glm_grs_npem3.df$grs_phenotype_name, 
																				levels=unique(results_glm_grs_npem3.df$grs_phenotype_name))
results_glm_grs_npem3.df$outcome_name <- factor(results_glm_grs_npem3.df$outcome_name, 
																					levels=unique(results_glm_grs_npem3.df$outcome_name))

write.table(results.df, "output/2019-03-28_glm_grs_results_npem3_covar.txt",
						quote=F, sep = "\t", row.names=F)

pdf("output/2019-03-28_glm_GRS_npem3_covar.pdf", height=20, width=20)
results_glm_grs_npem3.df %>%
	filter(!(outcome_name %in% c("npEM cluster 3", "npEM cluster 2/3", "Traj.A", 
													"Traj.B", "Traj.C1", "Traj.C2", "Traj.C1/2"))) %>%
	mutate(outcome_name = factor(outcome_name, levels=rev(levels(outcome_name)))) %>%
	ggplot(aes(x = grs_phenotype_name, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", good.signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	geom_hline(yintercept = c(2.5, 3.5, 5.5, 12.5, 14.5, 15.5)) +
	geom_vline(xintercept = c(2.5, 5.5, 7.5, 8.5, 11.5, 14.5, 17.5, 20.5))
dev.off()

results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 grs_phenotype_name = rep(grs_phenotype_name, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, "+ I(group == 3) + sex")),
								 data = master_GRS.df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_npem3.df <- results.df
results_glm_cont_grs_npem3.df$grs_phenotype_name <- factor(results_glm_cont_grs_npem3.df$grs_phenotype_name, 
																						 levels=unique(results_glm_cont_grs_npem3.df$grs_phenotype_name))
results_glm_cont_grs_npem3.df$outcome_name <- factor(results_glm_cont_grs_npem3.df$outcome_name, 
																							 levels=unique(results_glm_cont_grs_npem3.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_results_npem3_covar.txt",
						quote=F, sep = "\t", row.names=F)

pdf("output/2019-03-28_glm_cont_GRS_npem3_covar.pdf", height=8, width=20)
results_glm_cont_grs_npem3.df %>%
	mutate(outcome_name = factor(outcome_name, levels=rev(levels(outcome_name)))) %>%
	ggplot(aes(x = grs_phenotype_name, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, b.exact))) + 
 	scale_fill_gradientn(colours = divPalette.alt(100), lim=c(-1,1)) + 
	geom_text(aes(label = paste0("B = ", good.signif(b.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="beta", color="p < 0.05?") +
	geom_vline(xintercept = c(2.5, 5.5, 7.5, 8.5, 11.5, 14.5, 17.5, 20.5))
dev.off()

# Association analyses with npem Cluster - with GRS as covariate
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)

# Add extra rows for no GRS
results.df2 <- data.frame(outcome = outcome_list, 
													GRS = rep("none", length(outcome_list)),
													OR.exact = numeric(length(outcome_list)),
													lci.exact = numeric(length(outcome_list)),
													uci.exact =  numeric(length(outcome_list)),
													p.exact =  numeric(length(outcome_list)),
													outcome_name = outcome_name,
													n_total = numeric(length(outcome_list)),
													n_pheno = numeric(length(outcome_list)),
													score_labs = rep("none", times=length(outcome_list)),
													OR =  character(length(outcome_list)),
													p =  character(length(outcome_list)),
													stringsAsFactors = F)

for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ I(group == 3) + ",
																	 target_GRS, " + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}

for (i in 1:nrow(results.df2)) {
	outcome <- results.df2$outcome[i]
	results.df2$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(sex), na.rm=T))
	results.df2$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ I(group == 3) + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df2$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df2$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df2$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df2$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df2$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df2$p[i] <- signif(p, 2)
}

results_glm_npem3_grs.df <- rbind(results.df, results.df2)

results_glm_npem3_grs.df$score_labs <- factor(results_glm_npem3_grs.df$score_labs, 
																							levels=unique(results_glm_npem3_grs.df$score_labs))
results_glm_npem3_grs.df$outcome_name <- factor(results_glm_npem3_grs.df$outcome_name, 
																								levels=unique(results_glm_npem3_grs.df$outcome_name))

write.table(results.df, "output/2019-03-28_glm_npem3_results_grs_covar.txt",
						quote=F, sep = "\t", row.names=F)

pdf("output/2019-03-28_glm_npem3_GRS_covar.pdf", height=14, width=14)
results_glm_npem3_grs.df %>%
	filter(!(outcome_name %in% c("npEM cluster 3", "npEM cluster 2/3", "Traj.A", 
															 "Traj.B", "Traj.C1", "Traj.C2", "Traj.C1/2", 
															 "Any LRI by age 2", "Allergic sensitisation at age 5", 
															 "Allergic sensitisation at age 2"))) %>%
	mutate(outcome_name = factor(outcome_name, levels=rev(levels(outcome_name)))) %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-3, 3)) + 
	geom_text(aes(label = paste0("OR = ", good.signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	geom_hline(yintercept = c(1.5, 2.5, 7.5, 9.5, 10.5)) + 
	geom_vline(xintercept = 10.5)
dev.off()








# ------------------------------------------------------------------------------
# GLM associations with GRS - stratified by CAS
# ------------------------------------------------------------------------------

# Association analyses with GRS - CAS1
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_CAS1.df <- results.df
results_glm_grs_CAS1.df$score_labs <- factor(results_glm_grs_CAS1.df$score_labs, 
																				levels=unique(results_glm_grs_CAS1.df$score_labs))
results_glm_grs_CAS1.df$outcome_name <- factor(results_glm_grs_CAS1.df$outcome_name, 
																					levels=unique(results_glm_grs_CAS1.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_CAS1.txt",
						quote=F, sep = "\t", row.names=F)


results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_CAS1.df <- results.df
results_glm_cont_grs_CAS1.df$score_labs <- factor(results_glm_cont_grs_CAS1.df$score_labs, 
																						 levels=unique(results_glm_cont_grs_CAS1.df$score_labs))
results_glm_cont_grs_CAS1.df$outcome_name <- factor(results_glm_cont_grs_CAS1.df$outcome_name, 
																							 levels=unique(results_glm_cont_grs_CAS1.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_CAS1_results.txt",
						quote=F, sep = "\t", row.names=F)

# Association analyses with GRS - CAS2
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group == 2)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_CAS2.df <- results.df
results_glm_grs_CAS2.df$score_labs <- factor(results_glm_grs_CAS2.df$score_labs, 
																						 levels=unique(results_glm_grs_CAS2.df$score_labs))
results_glm_grs_CAS2.df$outcome_name <- factor(results_glm_grs_CAS2.df$outcome_name, 
																							 levels=unique(results_glm_grs_CAS2.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_CAS2.txt",
						quote=F, sep = "\t", row.names=F)

results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group == 2)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_CAS2.df <- results.df
results_glm_cont_grs_CAS2.df$score_labs <- factor(results_glm_cont_grs_CAS2.df$score_labs, 
																									levels=unique(results_glm_cont_grs_CAS2.df$score_labs))
results_glm_cont_grs_CAS2.df$outcome_name <- factor(results_glm_cont_grs_CAS2.df$outcome_name, 
																										levels=unique(results_glm_cont_grs_CAS2.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_CAS2_results.txt",
						quote=F, sep = "\t", row.names=F)


# Association analyses with GRS - CAS3
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group == 3)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_CAS3.df <- results.df
results_glm_grs_CAS3.df$score_labs <- factor(results_glm_grs_CAS3.df$score_labs, 
																						 levels=unique(results_glm_grs_CAS3.df$score_labs))
results_glm_grs_CAS3.df$outcome_name <- factor(results_glm_grs_CAS3.df$outcome_name, 
																							 levels=unique(results_glm_grs_CAS3.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_CAS3.txt",
						quote=F, sep = "\t", row.names=F)

results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group == 3)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_CAS3.df <- results.df
results_glm_cont_grs_CAS3.df$score_labs <- factor(results_glm_cont_grs_CAS3.df$score_labs, 
																									levels=unique(results_glm_cont_grs_CAS3.df$score_labs))
results_glm_cont_grs_CAS3.df$outcome_name <- factor(results_glm_cont_grs_CAS3.df$outcome_name, 
																										levels=unique(results_glm_cont_grs_CAS3.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_CAS3_results.txt",
						quote=F, sep = "\t", row.names=F)

# Combine
# ------------------------------------------------------------------------------
results_glm_grs_CAS1.df$group <- "CAS1"
results_glm_cont_grs_CAS1.df$group <- "CAS1"
results_glm_grs_CAS2.df$group <- "CAS2"
results_glm_cont_grs_CAS2.df$group <- "CAS2"
results_glm_grs_CAS3.df$group <- "CAS3"
results_glm_cont_grs_CAS3.df$group <- "CAS3"

results_glm_grs_npem.df <- rbind(results_glm_grs_CAS1.df, results_glm_grs_CAS2.df,
																 results_glm_grs_CAS3.df)
results_glm_cont_grs_npem.df <- rbind(results_glm_cont_grs_CAS1.df, results_glm_cont_grs_CAS2.df,
																			results_glm_cont_grs_CAS3.df)

pdf("output/2019-03-28_glm_GRS_CAS.pdf", height=12, width=48)
results_glm_grs_npem.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~group)
dev.off()

pdf("output/2019-03-28_glm_cont_GRS_CAS.pdf", height=12, width=48)
results_glm_cont_grs_npem.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, b.exact))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("B = ", signif(b.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="B", color="p < 0.05?") +
	facet_grid(.~group)
dev.off()

# Association analyses with GRS - CAS2/3
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group != 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_CAS2.3.df <- results.df
results_glm_grs_CAS2.3.df$score_labs <- factor(results_glm_grs_CAS2.3.df$score_labs, 
																						 levels=unique(results_glm_grs_CAS2.3.df$score_labs))
results_glm_grs_CAS2.3.df$outcome_name <- factor(results_glm_grs_CAS2.3.df$outcome_name, 
																							 levels=unique(results_glm_grs_CAS2.3.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_CAS2.3.txt",
						quote=F, sep = "\t", row.names=F)

results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(group != 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_CAS2.3.df <- results.df
results_glm_cont_grs_CAS2.3.df$score_labs <- factor(results_glm_cont_grs_CAS2.3.df$score_labs, 
																									levels=unique(results_glm_cont_grs_CAS2.3.df$score_labs))
results_glm_cont_grs_CAS2.3.df$outcome_name <- factor(results_glm_cont_grs_CAS2.3.df$outcome_name, 
																										levels=unique(results_glm_cont_grs_CAS2.3.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_CAS2.3_results.txt",
						quote=F, sep = "\t", row.names=F)

results_glm_grs_CAS1.df$group <- "CAS1"
results_glm_cont_grs_CAS1.df$group <- "CAS1"
results_glm_grs_CAS2.3.df$group <- "CAS2/3"
results_glm_cont_grs_CAS2.3.df$group <- "CAS2/3"

results_glm_grs_npem2.df <- rbind(results_glm_grs_CAS1.df, results_glm_grs_CAS2.3.df)
results_glm_cont_grs_npem2.df <- rbind(results_glm_cont_grs_CAS1.df, results_glm_cont_grs_CAS2.3.df)

pdf("output/2019-03-28_glm_GRS_CAS_alt.pdf", height=12, width=32)
results_glm_grs_npem2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~group)
dev.off()


pdf("output/2019-03-28_glm_cont_GRS_CAS_alt.pdf", height=12, width=48)
results_glm_cont_grs_npem2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, b.exact))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("B = ", signif(b.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="B", color="p < 0.05?") +
	facet_grid(.~group)
dev.off()


# ------------------------------------------------------------------------------
# GLM associations with GRS - stratified by allergic sensitisation at age 2
# ------------------------------------------------------------------------------

# Association analyses with GRS - ae2 == 0
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(ae2 == 0)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_ae2_no.df <- results.df
results_glm_grs_ae2_no.df$score_labs <- factor(results_glm_grs_ae2_no.df$score_labs, 
																						 levels=unique(results_glm_grs_ae2_no.df$score_labs))
results_glm_grs_ae2_no.df$outcome_name <- factor(results_glm_grs_ae2_no.df$outcome_name, 
																							 levels=unique(results_glm_grs_ae2_no.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_ae2_no.txt",
						quote=F, sep = "\t", row.names=F)


results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(ae2 == 0)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_ae2_no.df <- results.df
results_glm_cont_grs_ae2_no.df$score_labs <- factor(results_glm_cont_grs_ae2_no.df$score_labs, 
																									levels=unique(results_glm_cont_grs_ae2_no.df$score_labs))
results_glm_cont_grs_ae2_no.df$outcome_name <- factor(results_glm_cont_grs_ae2_no.df$outcome_name, 
																										levels=unique(results_glm_cont_grs_ae2_no.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_ae2_no_results.txt",
						quote=F, sep = "\t", row.names=F)


# Association analyses with GRS - ae2 == 1
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(ae2 == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_ae2_yes.df <- results.df
results_glm_grs_ae2_yes.df$score_labs <- factor(results_glm_grs_ae2_yes.df$score_labs, 
																							 levels=unique(results_glm_grs_ae2_yes.df$score_labs))
results_glm_grs_ae2_yes.df$outcome_name <- factor(results_glm_grs_ae2_yes.df$outcome_name, 
																								 levels=unique(results_glm_grs_ae2_yes.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_ae2_yes.txt",
						quote=F, sep = "\t", row.names=F)

results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(ae2 == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_ae2_yes.df <- results.df
results_glm_cont_grs_ae2_yes.df$score_labs <- factor(results_glm_cont_grs_ae2_yes.df$score_labs, 
																										levels=unique(results_glm_cont_grs_ae2_yes.df$score_labs))
results_glm_cont_grs_ae2_yes.df$outcome_name <- factor(results_glm_cont_grs_ae2_yes.df$outcome_name, 
																											levels=unique(results_glm_cont_grs_ae2_yes.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_ae2_yes_results.txt",
						quote=F, sep = "\t", row.names=F)

# Combine
# ------------------------------------------------------------------------------
results_glm_grs_ae2_no.df$ae2 <- "Not sensitised by age 2"
results_glm_cont_grs_ae2_no.df$ae2 <- "Not sensitised by age 2"
results_glm_grs_ae2_yes.df$ae2 <- "Sensitised by age 2 "
results_glm_cont_grs_ae2_yes.df$ae2 <- "Sensitised by age 2 "

results_glm_grs_ae2.df <- rbind(results_glm_grs_ae2_no.df, results_glm_grs_ae2_yes.df)
results_glm_cont_grs_ae2.df <- rbind(results_glm_cont_grs_ae2_no.df, results_glm_cont_grs_ae2_yes.df)

pdf("output/2019-03-28_glm_GRS_ae2.pdf", height=12, width=32)
results_glm_grs_ae2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~ae2)
dev.off()

pdf("output/2019-03-28_glm_cont_GRS_ae2.pdf", height=12, width=48)
results_glm_cont_grs_ae2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, b.exact))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("B = ", signif(b.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="B", color="p < 0.05?") +
	facet_grid(.~ae2)
dev.off()

# ------------------------------------------------------------------------------
# GLM associations with GRS - stratified by smoke exposure at age 2
# ------------------------------------------------------------------------------

# Association analyses with GRS - smoke2 == 0
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(smoke2 == 0)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_smoke2_no.df <- results.df
results_glm_grs_smoke2_no.df$score_labs <- factor(results_glm_grs_smoke2_no.df$score_labs, 
																							 levels=unique(results_glm_grs_smoke2_no.df$score_labs))
results_glm_grs_smoke2_no.df$outcome_name <- factor(results_glm_grs_smoke2_no.df$outcome_name, 
																								 levels=unique(results_glm_grs_smoke2_no.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_smoke2_no.txt",
						quote=F, sep = "\t", row.names=F)


results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(smoke2 == 0)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_smoke2_no.df <- results.df
results_glm_cont_grs_smoke2_no.df$score_labs <- factor(results_glm_cont_grs_smoke2_no.df$score_labs, 
																										 levels=unique(results_glm_cont_grs_smoke2_no.df$score_labs))
results_glm_cont_grs_smoke2_no.df$outcome_name <- factor(results_glm_cont_grs_smoke2_no.df$outcome_name, 
																											 levels=unique(results_glm_cont_grs_smoke2_no.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_smoke2_no_results.txt",
						quote=F, sep = "\t", row.names=F)

# Association analyses with GRS - smoke2 == 1
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(smoke2 == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_smoke2_yes.df <- results.df
results_glm_grs_smoke2_yes.df$score_labs <- factor(results_glm_grs_smoke2_yes.df$score_labs, 
																								levels=unique(results_glm_grs_smoke2_yes.df$score_labs))
results_glm_grs_smoke2_yes.df$outcome_name <- factor(results_glm_grs_smoke2_yes.df$outcome_name, 
																									levels=unique(results_glm_grs_smoke2_yes.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_smoke2_yes.txt",
						quote=F, sep = "\t", row.names=F)

results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(smoke2 == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_smoke2_yes.df <- results.df
results_glm_cont_grs_smoke2_yes.df$score_labs <- factor(results_glm_cont_grs_smoke2_yes.df$score_labs, 
																											 levels=unique(results_glm_cont_grs_smoke2_yes.df$score_labs))
results_glm_cont_grs_smoke2_yes.df$outcome_name <- factor(results_glm_cont_grs_smoke2_yes.df$outcome_name, 
																												 levels=unique(results_glm_cont_grs_smoke2_yes.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_smoke2_yes_results.txt",
						quote=F, sep = "\t", row.names=F)

# Combine
# ------------------------------------------------------------------------------

results_glm_grs_smoke2_no.df$smoke2 <- "No smoke exposure at age 2"
results_glm_cont_grs_smoke2_no.df$smoke2 <- "No smoke exposure at age 2"
results_glm_grs_smoke2_yes.df$smoke2 <- "Smoke exposure at age 2 "
results_glm_cont_grs_smoke2_yes.df$smoke2 <- "Smoke exposure at age 2 "

results_glm_grs_smoke2.df <- rbind(results_glm_grs_smoke2_no.df, results_glm_grs_smoke2_yes.df)
results_glm_cont_grs_smoke2.df <- rbind(results_glm_cont_grs_smoke2_no.df, results_glm_cont_grs_smoke2_yes.df)

pdf("output/2019-03-28_glm_GRS_smoke2.pdf", height=12, width=32)
results_glm_grs_smoke2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~smoke2)
dev.off()

pdf("output/2019-03-28_glm_cont_GRS_smoke2.pdf", height=12, width=48)
results_glm_cont_grs_smoke2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, b.exact))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("B = ", signif(b.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="B", color="p < 0.05?") +
	facet_grid(.~smoke2)
dev.off()

# ------------------------------------------------------------------------------
# GLM associations with GRS - stratified by cat exposure at age 2
# ------------------------------------------------------------------------------

# Association analyses with GRS - crcatexp2 == 0
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(crcatexp2 == 0)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_crcatexp2_no.df <- results.df
results_glm_grs_crcatexp2_no.df$score_labs <- factor(results_glm_grs_crcatexp2_no.df$score_labs, 
																									levels=unique(results_glm_grs_crcatexp2_no.df$score_labs))
results_glm_grs_crcatexp2_no.df$outcome_name <- factor(results_glm_grs_crcatexp2_no.df$outcome_name, 
																										levels=unique(results_glm_grs_crcatexp2_no.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_crcatexp2_no.txt",
						quote=F, sep = "\t", row.names=F)


results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(crcatexp2 == 0)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_crcatexp2_no.df <- results.df
results_glm_cont_grs_crcatexp2_no.df$score_labs <- factor(results_glm_cont_grs_crcatexp2_no.df$score_labs, 
																											 levels=unique(results_glm_cont_grs_crcatexp2_no.df$score_labs))
results_glm_cont_grs_crcatexp2_no.df$outcome_name <- factor(results_glm_cont_grs_crcatexp2_no.df$outcome_name, 
																												 levels=unique(results_glm_cont_grs_crcatexp2_no.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_crcatexp2_no_results.txt",
						quote=F, sep = "\t", row.names=F)


# Association analyses with GRS - crcatexp2 == 1
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(crcatexp2 == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_crcatexp2_yes.df <- results.df
results_glm_grs_crcatexp2_yes.df$score_labs <- factor(results_glm_grs_crcatexp2_yes.df$score_labs, 
																									 levels=unique(results_glm_grs_crcatexp2_yes.df$score_labs))
results_glm_grs_crcatexp2_yes.df$outcome_name <- factor(results_glm_grs_crcatexp2_yes.df$outcome_name, 
																										 levels=unique(results_glm_grs_crcatexp2_yes.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_crcatexp2_yes.txt",
						quote=F, sep = "\t", row.names=F)


results.df <- data.frame(outcome = rep(outcome_list_cont, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list_cont)),
												 b.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 outcome_name_cont = rep(outcome_name_cont, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list_cont)),
												 score_labs = rep(score_labs, times=length(outcome_list_cont)),
												 b =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list_cont)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(crcatexp2 == 1)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="gaussian")
	results.df$b.exact[i] <- b <- glm.out$coefficients[2]
	results.df$uci.exact[i] <- uci <- confint(glm.out)[2,2]
	results.df$lci.exact[i] <- lci <- confint(glm.out)[2,1]
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$b[i] <- paste0(good.signif(b, 2), " (", 
														good.signif(lci, 2), "-", 
														good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_cont_grs_crcatexp2_yes.df <- results.df
results_glm_cont_grs_crcatexp2_yes.df$score_labs <- factor(results_glm_cont_grs_crcatexp2_yes.df$score_labs, 
																											 levels=unique(results_glm_cont_grs_crcatexp2_yes.df$score_labs))
results_glm_cont_grs_crcatexp2_yes.df$outcome_name <- factor(results_glm_cont_grs_crcatexp2_yes.df$outcome_name, 
																												 levels=unique(results_glm_cont_grs_crcatexp2_yes.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_cont_grs_crcatexp2_yes_results.txt",
						quote=F, sep = "\t", row.names=F)

# Combine
# ------------------------------------------------------------------------------
results_glm_grs_crcatexp2_no.df$crcatexp2 <- "No cat exposure at age 2"
results_glm_cont_grs_crcatexp2_no.df$crcatexp2 <- "No cat exposure at age 2"
results_glm_grs_crcatexp2_yes.df$crcatexp2 <- "Cat exposure at age 2"
results_glm_cont_grs_crcatexp2_yes.df$crcatexp2 <- "Cat exposure at age 2"

results_glm_grs_crcatexp2.df <- rbind(results_glm_grs_crcatexp2_no.df, results_glm_grs_crcatexp2_yes.df)
results_glm_cont_grs_crcatexp2.df <- rbind(results_glm_cont_grs_crcatexp2_no.df, results_glm_cont_grs_crcatexp2_yes.df)

results_glm_grs_crcatexp2.df$crcatexp2 <- 
	factor(results_glm_grs_crcatexp2.df$crcatexp2, 
				 levels=c("No cat exposure at age 2", "Cat exposure at age 2"))
results_glm_cont_grs_crcatexp2.df$crcatexp2 <- 
	factor(results_glm_cont_grs_crcatexp2.df$crcatexp2, 
				 levels=c("No cat exposure at age 2", "Cat exposure at age 2"))


pdf("output/2019-03-28_glm_GRS_crcatexp2.pdf", height=12, width=32)
results_glm_grs_crcatexp2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-2, 2)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~crcatexp2)
dev.off()

pdf("output/2019-03-28_glm_cont_GRS_crcatexp2.pdf", height=12, width=32)
results_glm_cont_grs_crcatexp2.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, b.exact))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("B = ", signif(b.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="B", color="p < 0.05?") +
	facet_grid(.~crcatexp2)
dev.off()


# ------------------------------------------------------------------------------
# GLM associations with GRS - stratified by trajectory
# ------------------------------------------------------------------------------

# Association analyses with GRS - trajectory == Traj.A
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(traj == "Traj.A")
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_traj_A.df <- results.df
results_glm_grs_traj_A.df$score_labs <- factor(results_glm_grs_traj_A.df$score_labs, 
																										 levels=unique(results_glm_grs_traj_A.df$score_labs))
results_glm_grs_traj_A.df$outcome_name <- factor(results_glm_grs_traj_A.df$outcome_name, 
																											 levels=unique(results_glm_grs_traj_A.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_traj_A.txt",
						quote=F, sep = "\t", row.names=F)

# Association analyses with GRS - traj == B
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(traj == "Traj.B")
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_traj_B.df <- results.df
results_glm_grs_traj_B.df$score_labs <- factor(results_glm_grs_traj_B.df$score_labs, 
																											levels=unique(results_glm_grs_traj_B.df$score_labs))
results_glm_grs_traj_B.df$outcome_name <- factor(results_glm_grs_traj_B.df$outcome_name, 
																												levels=unique(results_glm_grs_traj_B.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_traj_B.txt",
						quote=F, sep = "\t", row.names=F)

# Association analyses with GRS - traj == C1
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(traj == "Traj.C1")
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_traj_C1.df <- results.df
results_glm_grs_traj_C1.df$score_labs <- factor(results_glm_grs_traj_C1.df$score_labs, 
																							 levels=unique(results_glm_grs_traj_C1.df$score_labs))
results_glm_grs_traj_C1.df$outcome_name <- factor(results_glm_grs_traj_C1.df$outcome_name, 
																								 levels=unique(results_glm_grs_traj_C1.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_traj_C1.txt",
						quote=F, sep = "\t", row.names=F)

# Association analyses with GRS - traj == C2
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(traj == "Traj.C2")
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_traj_C2.df <- results.df
results_glm_grs_traj_C2.df$score_labs <- factor(results_glm_grs_traj_C2.df$score_labs, 
																								levels=unique(results_glm_grs_traj_C2.df$score_labs))
results_glm_grs_traj_C2.df$outcome_name <- factor(results_glm_grs_traj_C2.df$outcome_name, 
																									levels=unique(results_glm_grs_traj_C2.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_traj_C2.txt",
						quote=F, sep = "\t", row.names=F)

# Association analyses with GRS - traj == all C
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(traj %in% c("Traj.C1", "Traj.C2"))
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_traj_C.df <- results.df
results_glm_grs_traj_C.df$score_labs <- factor(results_glm_grs_traj_C.df$score_labs, 
																								levels=unique(results_glm_grs_traj_C.df$score_labs))
results_glm_grs_traj_C.df$outcome_name <- factor(results_glm_grs_traj_C.df$outcome_name, 
																									levels=unique(results_glm_grs_traj_C.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_traj_C.txt",
						quote=F, sep = "\t", row.names=F)


results_glm_grs_traj_A.df$traj <- "Traj.A"
results_glm_grs_traj_B.df$traj <- "Traj.B"
results_glm_grs_traj_C1.df$traj <- "Traj.C1"
results_glm_grs_traj_C2.df$traj <- "Traj.C2"
results_glm_grs_traj_C.df$traj <- "Traj.C1/2"


results_glm_grs_traj.df <- rbind(results_glm_grs_traj_A.df,
																 results_glm_grs_traj_B.df,
																 results_glm_grs_traj_C1.df,
																 results_glm_grs_traj_C2.df)
results_glm_grs_traj.df$traj <- 
	factor(results_glm_grs_traj.df$traj, 
				 levels=c("Traj.A", "Traj.B", "Traj.C1", "Traj.C2"))


pdf("output/2019-03-28_glm_GRS_traj.pdf", height=10, width=32)
results_glm_grs_traj.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-2, 2)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~traj)
dev.off()

results_glm_grs_traj.alt.df <- rbind(results_glm_grs_traj_A.df,
																 results_glm_grs_traj_B.df,
																 results_glm_grs_traj_C.df)
results_glm_grs_traj.df$traj <- 
	factor(results_glm_grs_traj.df$traj, 
				 levels=c("Traj.A", "Traj.B", "Traj.C1/2"))


pdf("output/2019-03-28_glm_GRS_traj.alt.pdf", height=10, width=32)
results_glm_grs_traj.alt.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-2, 2)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~traj)
dev.off()


# ------------------------------------------------------------------------------
# GLM associations with GRS - stratified by MPG
# ------------------------------------------------------------------------------

# Association analyses with GRS - illness.mpg = Low
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(illness.mpg.healthy.age2.bin == "Low")
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_illness.mpg.low.df <- results.df
results_glm_grs_illness.mpg.low.df$score_labs <- factor(results_glm_grs_illness.mpg.low.df$score_labs, 
																							 levels=unique(results_glm_grs_illness.mpg.low.df$score_labs))
results_glm_grs_illness.mpg.low.df$outcome_name <- factor(results_glm_grs_illness.mpg.low.df$outcome_name, 
																								 levels=unique(results_glm_grs_illness.mpg.low.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_illness.mpg.low.txt",
						quote=F, sep = "\t", row.names=F)

# Association analyses with GRS - illness.mpg = High
# ------------------------------------------------------------------------------
results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
df <- master_GRS.df %>% filter(illness.mpg.healthy.age2.bin == "High")
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " + sex")),
								 data = df,
								 family="binomial")
	results.df$OR.exact[i] <- or <- exp(glm.out$coefficients[2])
	results.df$uci.exact[i] <- uci <- exp(confint(glm.out)[2,2])
	results.df$lci.exact[i] <- lci <- exp(confint(glm.out)[2,1])
	results.df$p.exact[i] <- p <- coef(summary(glm.out))[,4][2]
	results.df$OR[i] <- paste0(good.signif(or, 2), " (", 
														 good.signif(lci, 2), "-", 
														 good.signif(uci, 2), ")")
	results.df$p[i] <- signif(p, 2)
}
results_glm_grs_illness.mpg.high.df <- results.df
results_glm_grs_illness.mpg.high.df$score_labs <- factor(results_glm_grs_illness.mpg.high.df$score_labs, 
																												levels=unique(results_glm_grs_illness.mpg.high.df$score_labs))
results_glm_grs_illness.mpg.high.df$outcome_name <- factor(results_glm_grs_illness.mpg.high.df$outcome_name, 
																													levels=unique(results_glm_grs_illness.mpg.high.df$outcome_name))
write.table(results.df, "output/2019-03-28_glm_grs_results_illness.mpg.high.txt",
						quote=F, sep = "\t", row.names=F)

results_glm_grs_illness.mpg.low.df$illness.mpg.bin <- "Low"
results_glm_grs_illness.mpg.high.df$illness.mpg.bin <- "High"

results_glm_grs_illness.mpg.df <- rbind(results_glm_grs_illness.mpg.low.df,
																				results_glm_grs_illness.mpg.high.df)
results_glm_grs_illness.mpg.df$illness.mpg.bin <- 
	factor(results_glm_grs_illness.mpg.df$illness.mpg.bin, 
				 levels=c("Low", "High"))


pdf("output/2019-03-28_glm_GRS_illness.mpg.pdf", height=10, width=32)
results_glm_grs_illness.mpg.df %>%
	ggplot(aes(x = score_labs, y = outcome_name)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-2, 2)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?") +
	facet_grid(.~illness.mpg.bin)
dev.off()

# ------------------------------------------------------------------------------
# GLM - Test Interaction with group
# ------------------------------------------------------------------------------

results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " * I(group %in% c(2,3)) + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df$OR.score.exact[i] <- or.score <- exp(glm.out$coefficients[2])
	results.df$uci.score.exact[i] <- uci.score <- exp(confint(glm.out)[2,2])
	results.df$lci.score.exact[i] <- lci.score <- exp(confint(glm.out)[2,1])
	results.df$p.score.exact[i] <- p.score <- coef(summary(glm.out))[,4][2]
	results.df$OR.covar.exact[i] <- or.covar <- exp(glm.out$coefficients[3])
	results.df$uci.covar.exact[i] <- uci.covar <- exp(confint(glm.out)[3,2])
	results.df$lci.covar.exact[i] <- lci.covar <- exp(confint(glm.out)[3,1])
	results.df$p.covar.exact[i] <- p.covar <- coef(summary(glm.out))[,4][3]
	results.df$OR.interact.exact[i] <- or.interact <- exp(glm.out$coefficients[5])
	results.df$uci.interact.exact[i] <- uci.interact <- exp(confint(glm.out)[5,2])
	results.df$lci.interact.exact[i] <- lci.interact <- exp(confint(glm.out)[5,1])
	results.df$p.interact.exact[i] <- p.interact <- coef(summary(glm.out))[,4][5]
	results.df$OR.score[i] <- paste0(good.signif(or.score, 2), " (",
																	 good.signif(lci.score, 2), "-",
																	 good.signif(uci.score, 2), ")")
	results.df$p.score[i] <- signif(p.score, 2)
	results.df$OR.covar[i] <- paste0(good.signif(or.covar, 2), " (",
																	 good.signif(lci.covar, 2), "-",
																	 good.signif(uci.covar, 2), ")")
	results.df$p.covar[i] <- signif(p.covar, 2)
	results.df$OR.interact[i] <- paste0(good.signif(or.interact, 2), " (",
																			good.signif(lci.interact, 2), "-",
																			good.signif(uci.interact, 2), ")")
	results.df$p.interact[i] <- signif(p.interact, 2)
}
results_glm_grs_interact.df <- results.df
results_glm_grs_interact.df$score_labs <- factor(results_glm_grs.df$score_labs, 
																								 levels=unique(results_glm_grs.df$score_labs))
results_glm_grs_interact.df$outcome_name <- factor(results_glm_grs.df$outcome_name, 
																									 levels=unique(results_glm_grs.df$outcome_name))

results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.score.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05 &
																		results_glm_grs_interact.df$p.score.exact < 0.05),]

results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_PC1_std"),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_allergic_std"),]

write.csv(results_glm_grs_interact.df[, c((ncol(results_glm_grs_interact.df)-7): 
																						ncol(results_glm_grs_interact.df))],
					"output/2019-03-28_glm_grs_interact_npEM_2.3.csv",
					quote=F, row.names=F)

# ------------------------------------------------------------------------------
# GLM - Test Interaction with crcatexp2
# ------------------------------------------------------------------------------

results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " * crcatexp2 + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df$OR.score.exact[i] <- or.score <- exp(glm.out$coefficients[2])
	results.df$uci.score.exact[i] <- uci.score <- exp(confint(glm.out)[2,2])
	results.df$lci.score.exact[i] <- lci.score <- exp(confint(glm.out)[2,1])
	results.df$p.score.exact[i] <- p.score <- coef(summary(glm.out))[,4][2]
	results.df$OR.covar.exact[i] <- or.covar <- exp(glm.out$coefficients[3])
	results.df$uci.covar.exact[i] <- uci.covar <- exp(confint(glm.out)[3,2])
	results.df$lci.covar.exact[i] <- lci.covar <- exp(confint(glm.out)[3,1])
	results.df$p.covar.exact[i] <- p.covar <- coef(summary(glm.out))[,4][3]
	results.df$OR.interact.exact[i] <- or.interact <- exp(glm.out$coefficients[5])
	results.df$uci.interact.exact[i] <- uci.interact <- exp(confint(glm.out)[5,2])
	results.df$lci.interact.exact[i] <- lci.interact <- exp(confint(glm.out)[5,1])
	results.df$p.interact.exact[i] <- p.interact <- coef(summary(glm.out))[,4][5]
	results.df$OR.score[i] <- paste0(good.signif(or.score, 2), " (",
														 good.signif(lci.score, 2), "-",
														 good.signif(uci.score, 2), ")")
	results.df$p.score[i] <- signif(p.score, 2)
	results.df$OR.covar[i] <- paste0(good.signif(or.covar, 2), " (",
																	 good.signif(lci.covar, 2), "-",
																	 good.signif(uci.covar, 2), ")")
	results.df$p.covar[i] <- signif(p.covar, 2)
	results.df$OR.interact[i] <- paste0(good.signif(or.interact, 2), " (",
																	 good.signif(lci.interact, 2), "-",
																	 good.signif(uci.interact, 2), ")")
	results.df$p.interact[i] <- signif(p.interact, 2)
}
results_glm_grs_interact.df <- results.df
results_glm_grs_interact.df$score_labs <- factor(results_glm_grs.df$score_labs, 
																				levels=unique(results_glm_grs.df$score_labs))
results_glm_grs_interact.df$outcome_name <- factor(results_glm_grs.df$outcome_name, 
																					levels=unique(results_glm_grs.df$outcome_name))

results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.score.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05 &
																		results_glm_grs_interact.df$p.score.exact < 0.05),]

results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_PC1_std"),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_allergic_std"),]

write.csv(results_glm_grs_interact.df[, c((ncol(results_glm_grs_interact.df)-7): 
																						ncol(results_glm_grs_interact.df))],
						"output/2019-03-28_glm_grs_interact_crcatexp.csv",
						quote=F, row.names=F)


# ------------------------------------------------------------------------------
# GLM - Test Interaction with traj
# ------------------------------------------------------------------------------

results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " * I(traj %in% c('Traj.C1', 'Traj.C2')) + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df$OR.score.exact[i] <- or.score <- exp(glm.out$coefficients[2])
	results.df$uci.score.exact[i] <- uci.score <- exp(confint(glm.out)[2,2])
	results.df$lci.score.exact[i] <- lci.score <- exp(confint(glm.out)[2,1])
	results.df$p.score.exact[i] <- p.score <- coef(summary(glm.out))[,4][2]
	results.df$OR.covar.exact[i] <- or.covar <- exp(glm.out$coefficients[3])
	results.df$uci.covar.exact[i] <- uci.covar <- exp(confint(glm.out)[3,2])
	results.df$lci.covar.exact[i] <- lci.covar <- exp(confint(glm.out)[3,1])
	results.df$p.covar.exact[i] <- p.covar <- coef(summary(glm.out))[,4][3]
	results.df$OR.interact.exact[i] <- or.interact <- exp(glm.out$coefficients[5])
	results.df$uci.interact.exact[i] <- uci.interact <- exp(confint(glm.out)[5,2])
	results.df$lci.interact.exact[i] <- lci.interact <- exp(confint(glm.out)[5,1])
	results.df$p.interact.exact[i] <- p.interact <- coef(summary(glm.out))[,4][5]
	results.df$OR.score[i] <- paste0(good.signif(or.score, 2), " (",
																	 good.signif(lci.score, 2), "-",
																	 good.signif(uci.score, 2), ")")
	results.df$p.score[i] <- signif(p.score, 2)
	results.df$OR.covar[i] <- paste0(good.signif(or.covar, 2), " (",
																	 good.signif(lci.covar, 2), "-",
																	 good.signif(uci.covar, 2), ")")
	results.df$p.covar[i] <- signif(p.covar, 2)
	results.df$OR.interact[i] <- paste0(good.signif(or.interact, 2), " (",
																			good.signif(lci.interact, 2), "-",
																			good.signif(uci.interact, 2), ")")
	results.df$p.interact[i] <- signif(p.interact, 2)
}
results_glm_grs_interact.df <- results.df
results_glm_grs_interact.df$score_labs <- factor(results_glm_grs.df$score_labs, 
																								 levels=unique(results_glm_grs.df$score_labs))
results_glm_grs_interact.df$outcome_name <- factor(results_glm_grs.df$outcome_name, 
																									 levels=unique(results_glm_grs.df$outcome_name))

results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.score.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05 &
																		results_glm_grs_interact.df$p.score.exact < 0.05),]

results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_PC1_std"),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_allergic_std"),]

write.csv(results_glm_grs_interact.df[, c((ncol(results_glm_grs_interact.df)-7): 
																						ncol(results_glm_grs_interact.df))],
					"output/2019-03-28_glm_grs_traj.c.mpg.csv",
					quote=F, row.names=F)

# ------------------------------------------------------------------------------
# GLM - Test Interaction with illness.mpg.healthy.age2.bin
# ------------------------------------------------------------------------------

results.df <- data.frame(outcome = rep(outcome_list, each=length(grs_phenotype_std)), 
												 GRS = rep(grs_phenotype_std, times=length(outcome_list)),
												 OR.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.score.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 lci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 uci.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact.exact =  numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_total = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 n_pheno = numeric(length(grs_phenotype_std) * length(outcome_list)),
												 outcome_name = rep(outcome_name, each=length(grs_phenotype_std)),
												 score_labs = rep(score_labs, times=length(outcome_list)),
												 OR.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.score =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.covar =  character(length(grs_phenotype_std) * length(outcome_list)),
												 OR.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 p.interact =  character(length(grs_phenotype_std) * length(outcome_list)),
												 stringsAsFactors = F)
for (i in 1:nrow(results.df)) {
	target_GRS <- results.df$GRS[i]
	outcome <- results.df$outcome[i]
	results.df$n_total[i] <- with(master_GRS.df, sum(
		!is.na(eval(parse(text = outcome))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	results.df$n_pheno[i] <- with(master_GRS.df, sum(
		eval(parse(text = paste0("(", outcome, ") == 1"))) & 
			!is.na(eval(parse(text = target_GRS))) &
			!is.na(sex), na.rm=T))
	glm.out <- glm(as.formula(paste0(outcome, " ~ ",
																	 target_GRS, " * illness.mpg.healthy.age2.bin + sex")),
								 data = master_GRS.df,
								 family="binomial")
	results.df$OR.score.exact[i] <- or.score <- exp(glm.out$coefficients[2])
	results.df$uci.score.exact[i] <- uci.score <- exp(confint(glm.out)[2,2])
	results.df$lci.score.exact[i] <- lci.score <- exp(confint(glm.out)[2,1])
	results.df$p.score.exact[i] <- p.score <- coef(summary(glm.out))[,4][2]
	results.df$OR.covar.exact[i] <- or.covar <- exp(glm.out$coefficients[3])
	results.df$uci.covar.exact[i] <- uci.covar <- exp(confint(glm.out)[3,2])
	results.df$lci.covar.exact[i] <- lci.covar <- exp(confint(glm.out)[3,1])
	results.df$p.covar.exact[i] <- p.covar <- coef(summary(glm.out))[,4][3]
	results.df$OR.interact.exact[i] <- or.interact <- exp(glm.out$coefficients[5])
	results.df$uci.interact.exact[i] <- uci.interact <- exp(confint(glm.out)[5,2])
	results.df$lci.interact.exact[i] <- lci.interact <- exp(confint(glm.out)[5,1])
	results.df$p.interact.exact[i] <- p.interact <- coef(summary(glm.out))[,4][5]
	results.df$OR.score[i] <- paste0(good.signif(or.score, 2), " (",
																	 good.signif(lci.score, 2), "-",
																	 good.signif(uci.score, 2), ")")
	results.df$p.score[i] <- signif(p.score, 2)
	results.df$OR.covar[i] <- paste0(good.signif(or.covar, 2), " (",
																	 good.signif(lci.covar, 2), "-",
																	 good.signif(uci.covar, 2), ")")
	results.df$p.covar[i] <- signif(p.covar, 2)
	results.df$OR.interact[i] <- paste0(good.signif(or.interact, 2), " (",
																			good.signif(lci.interact, 2), "-",
																			good.signif(uci.interact, 2), ")")
	results.df$p.interact[i] <- signif(p.interact, 2)
}
results_glm_grs_interact.df <- results.df
results_glm_grs_interact.df$score_labs <- factor(results_glm_grs.df$score_labs, 
																								 levels=unique(results_glm_grs.df$score_labs))
results_glm_grs_interact.df$outcome_name <- factor(results_glm_grs.df$outcome_name, 
																									 levels=unique(results_glm_grs.df$outcome_name))

results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.score.exact < 0.05),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$p.interact.exact < 0.05 &
																		results_glm_grs_interact.df$p.score.exact < 0.05),]

results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_PC1_std"),]
results_glm_grs_interact.df[which(results_glm_grs_interact.df$GRS %in% "SCORESUM_allergic_std"),]

write.csv(results_glm_grs_interact.df[, c((ncol(results_glm_grs_interact.df)-7): 
																						ncol(results_glm_grs_interact.df))],
					"output/2019-03-28_glm_grs_interact_illness.mpg.csv",
					quote=F, row.names=F)

# ------------------------------------------------------------------------------
# Tertile variables
# ------------------------------------------------------------------------------

grs_phenotype_tert <- c("SCORESUM_allergic_tert", "SCORESUM_asthma_tert",
									"SCORESUM_adult_asthma_tert", "SCORESUM_childhood_asthma_tert",
									"SCORESUM_ar_tert", "SCORESUM_as_tert", "SCORESUM_nar_tert",
									"SCORESUM_combined_tert", "SCORESUM_PC1_tert")
tertile_labs <- c("Allergic disease GRS tertile", "Asthma GRS tertile",
									"Adult asthma GRS tertile", "Childhood asthma GRS tertile",
									"Allergic rhinitis GRS tertile", "Allergic sensitisation GRS tertile",
									"Non-allergic rhinitis GRS tertile", "Combined GRS tertile",
									"PC1 of all GRS tertile")

# ------------------------------------------------------------------------------
# Crwz vs. scores
# ------------------------------------------------------------------------------

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.df <- results.df

# Crwz vs. scores - stratified by CAS group
# ------------------------------------------------------------------------------

# CAS 1

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(group == 1)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.CAS1.df <- results.df

# CAS2

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(group == 2)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.CAS2.df <- results.df

# CAS3

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(group == 3)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.CAS3.df <- results.df

# CAS2/3

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(group %in% c(2,3))

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.CAS2.3.df <- results.df

crwz.score.grs.results.CAS1.df$group <- "CAS1"
crwz.score.grs.results.CAS2.df$group <- "CAS2"
crwz.score.grs.results.CAS3.df$group <- "CAS3"
crwz.score.grs.results.CAS2.3.df$group <- "CAS2/3"

crwz.score.grs.results.CAS.df <- rbind(crwz.score.grs.results.CAS1.df, 
																			 crwz.score.grs.results.CAS2.df,
																			 crwz.score.grs.results.CAS3.df)

pdf("output/2019-03-28_gee_GRS_crwz_npem.pdf", height=10, width=16)
crwz.score.grs.results.CAS.df %>%
	ggplot(aes(x = GRS_name, y = group)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?")
dev.off()

crwz.score.grs.results.CAS.alt.df <- rbind(crwz.score.grs.results.CAS1.df, 
																			 crwz.score.grs.results.CAS2.3.df)

pdf("output/2019-03-28_gee_GRS_crwz_npem_alt.pdf", height=10, width=16)
crwz.score.grs.results.CAS.alt.df %>%
	ggplot(aes(x = GRS_name, y = group)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?")
dev.off()


# Crwz vs. scores - stratified by ae2
# ------------------------------------------------------------------------------

# ae2 == 0

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(ae2 == 0)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.ae2.no.df <- results.df

# ae2 == 1

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(ae2 == 1)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.ae2.yes.df <- results.df

crwz.score.grs.results.ae2.no.df$group <- "No sensitisation by age 2"
crwz.score.grs.results.ae2.yes.df$group <- "Sensitisation by age 2"

crwz.score.grs.results.ae2.df <- rbind(crwz.score.grs.results.ae2.no.df, 
																			 crwz.score.grs.results.ae2.yes.df)

pdf("output/2019-03-28_gee_GRS_crwz_ae2.pdf", height=10, width=16)
crwz.score.grs.results.ae2.df %>%
	ggplot(aes(x = GRS_name, y = group)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?")
dev.off()

# Crwz vs. scores - stratified by smoke2
# ------------------------------------------------------------------------------

# smoke2 == 0

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(smoke2 == 0)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.smoke2.no.df <- results.df

# smoke2 == 1

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(smoke2 == 1)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.smoke2.yes.df <- results.df

crwz.score.grs.results.smoke2.no.df$group <- "No smoke exposure at age 2"
crwz.score.grs.results.smoke2.yes.df$group <- "Smoke exposure at age 2"

crwz.score.grs.results.smoke2.df <- rbind(crwz.score.grs.results.smoke2.no.df, 
																			 crwz.score.grs.results.smoke2.yes.df)

pdf("output/2019-03-28_gee_GRS_crwz_smoke2.pdf", height=10, width=16)
crwz.score.grs.results.smoke2.df %>%
	ggplot(aes(x = GRS_name, y = group)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?")
dev.off()


# Crwz vs. scores - stratified by crcatexp2
# ------------------------------------------------------------------------------

# crcatexp2 == 0

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(crcatexp2 == 0)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.crcatexp2.no.df <- results.df

# crcatexp2 == 1

results.df <- data.frame(outcome = rep("crwz", each=length(grs_phenotype_std)), 
												 GRS = grs_phenotype_std,
												 OR.exact =  numeric(length(grs_phenotype_std)),
												 lci.exact =  numeric(length(grs_phenotype_std)),
												 uci.exact =  numeric(length(grs_phenotype_std)),
												 p.exact =  numeric(length(grs_phenotype_std)),
												 outcome_name = rep("Current wheeze", length(grs_phenotype_std)),
												 GRS_name = score_labs,
												 OR = character(length(grs_phenotype_std)),
												 p = character(length(grs_phenotype_std)),
												 stringsAsFactors = F)

master_GRS.df2 <- master_GRS.df %>% filter(crcatexp2 == 1)

for (i in 1:length(grs_phenotype_std)) {
	var_score <- grs_phenotype_std[i]
	var_score_lab <- score_labs[i]
	df <- with(master_GRS.df2, 
						 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
						 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
						 					 var_score = get(var_score))) %>%
		melt(id.vars = c("SUBJ_ID", "var_score"),
				 variable.name = "timepoint",
				 value.name = "crwz") %>%
		mutate(timepoint = as.numeric(as.character(
			factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
						 labels = c(1, 2, 3, 4, 5, 10)))))
	df <- merge(df, with(master_GRS.df2, data.frame(SUBJ_ID, sex)), by="SUBJ_ID", 
							all.x = T)
	gee.out <- tryCatch(gee(crwz ~ as.numeric(var_score) + sex + timepoint,
													id = SUBJ_ID,
													data = df,
													family = "binomial",
													corstr = "unstructured"),
											error = function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(gee.out)
	}
	results.df$OR.exact[i] <- OR.exact <- exp(s$coefficients[2,1])
	results.df$lci.exact[i] <- lci.exact <- exp(s$coefficients[2,1] - 1.96*(s$coefficients[2,4]))
	results.df$uci.exact[i] <- uci.exact <- exp(s$coefficients[2,1] + 1.96*(s$coefficients[2,4]))
	results.df$p.exact[i] <- p.exact <- 2 * pnorm(abs(s$coefficients[2,5]), lower.tail = FALSE)
	results.df$OR[i] <- paste0(good.signif(OR.exact, 2), " (", 
														 good.signif(lci.exact, 2), "-", 
														 good.signif(uci.exact, 2), ")")
	results.df$p[i] <- signif(p.exact, 2)
}
crwz.score.grs.results.crcatexp2.yes.df <- results.df

crwz.score.grs.results.crcatexp2.no.df$group <- "No cat exposure at age 2"
crwz.score.grs.results.crcatexp2.yes.df$group <- "Cat exposure at age 2"

crwz.score.grs.results.crcatexp2.df <- rbind(crwz.score.grs.results.crcatexp2.no.df, 
																					crwz.score.grs.results.crcatexp2.yes.df)

pdf("output/2019-03-28_gee_GRS_crwz_crcatexp2.pdf", height=10, width=16)
crwz.score.grs.results.crcatexp2.df %>%
	ggplot(aes(x = GRS_name, y = group)) + 
	geom_tile(aes(fill = ifelse(p.exact > 0.05, NA, log10(OR.exact)))) + 
	scale_fill_gradientn(colours = divPalette.alt(100), limits=c(-1, 1)) + 
	geom_text(aes(label = paste0("OR = ", signif(OR.exact,2), 
															 "\n p = ", signif(p.exact, 2)),
								color = (p.exact < 0.05))) +
	scale_color_manual(values=c("grey", "black")) +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
	labs(x = "GRS", y = "Outcome", fill="Log10(OR)", color="p < 0.05?")
dev.off()

# ------------------------------------------------------------------------------
# Mediation analysis
# ------------------------------------------------------------------------------
# install.packages("mediation")
# library(mediation)

# Model 1: X --> (c) --> Y
# Model 2: X --> (a) --> M
# Model 3: X --> (c') --> Y
# 					 M --> (b) ---^

# GRS --> npem --> wheeze/asthma at age 5
# ------------------------------------------------------------------------------

# Crwz5 ~ allergic disease GRS
# ------------------------------------

master_GRS.df$group3 <- as.numeric(master_GRS.df$group == 3)

df <- master_GRS.df %>%
	filter(!is.na(group3), !is.na(SCORESUM_allergic_std), !is.na(sex), !is.na(crwz5)) 

# Model.m = mediator ~ predictor + covariates
model2 <- glm(group3 ~ SCORESUM_allergic_std + sex, 
								data = df, family="binomial")

# Model.y = outcome ~ mediator + predictor + covariates
model3 <- glm(crwz5 ~ group3 + SCORESUM_allergic_std + sex,
						data = df, family="binomial")

# Mediation analysis using bootstrap & simulations to estimate CI
# Number of simulations
ndraws <- 1000

# "treat" = treatment, i.e. exposure or predictor
# "mediator" = mediator

medout_bootstrap <- mediate(model.m = model2, model.y = model3, boot = T, 
														treat = "SCORESUM_allergic_std", mediator = "group3", 
														sims = ndraws)

summ_bootstrap <- summary(medout_bootstrap)

# What you will be interested in:
# Mediation effects and pvalues (ab)
# ACME = average causal mediation effect
summ_bootstrap$d.avg # ACME_boot = 0.043
summ_bootstrap$d.avg.p  # ACME_pval_boot = 0.074# 
# Proportion of mediation
summ_bootstrap$n.avg # = 0.40

# Might be helpful, but i dont use these because when glm is involved, the effect sizes were not easy to interpret
# Direct effects (c')
# ADE = Average direct effect
summ_bootstrap$z.avg # ADE_boot = 0.066
summ_bootstrap$z.avg.p # ADE_pval_boot = 0.068#
# Total effects (c' + a*b)
summ_bootstrap$tau.coef # 0.10
summ_bootstrap$tau.p # 0.008**

# Cras_5 ~ allergic disease GRS
# ------------------------------------

df <- master_GRS.df %>%
	filter(!is.na(group3), !is.na(SCORESUM_allergic_std), !is.na(sex), !is.na(cras_5)) 

# Model.m = mediator ~ predictor + covariates
model2 <- glm(group3 ~ SCORESUM_allergic_std + sex, 
							data = df, family="binomial")

# Model.y = outcome ~ mediator + predictor + covariates
model3 <- glm(cras_5 ~ group3 + SCORESUM_allergic_std + sex,
							data = df, family="binomial")

# Mediation analysis using bootstrap & simulations to estimate CI
# Number of simulations
ndraws <- 1000

# "treat" = treatment, i.e. exposure or predictor
# "mediator" = mediator

medout_bootstrap <- mediate(model.m = model2, model.y = model3, boot = T, 
														treat = "SCORESUM_allergic_std", mediator = "group3", 
														sims = ndraws)

summ_bootstrap <- summary(medout_bootstrap)

# What you will be interested in:
# Mediation effects and pvalues (ab)
# ACME = average causal mediation effect
summ_bootstrap$d.avg # ACME_boot = 0.028
summ_bootstrap$d.avg.p  # ACME_pval_boot = 0.084# 
# Proportion of mediation
summ_bootstrap$n.avg # = 0.35

# Might be helpful, but i dont use these because when glm is involved, the effect sizes were not easy to interpret
# Direct effects (c')
# ADE = Average direct effect
summ_bootstrap$z.avg # ADE_boot = 0.053
summ_bootstrap$z.avg.p # ADE_pval_boot = 0.14
# Total effects (c' + a*b)
summ_bootstrap$tau.coef # 0081
summ_bootstrap$tau.p # 0.008**


# Crwz5 ~ PC1 GRS
# ------------------------------------

master_GRS.df$group3 <- as.numeric(master_GRS.df$group == 3)

df <- master_GRS.df %>%
	filter(!is.na(group3), !is.na(SCORESUM_PC1_std), !is.na(sex), !is.na(crwz5)) 

# Model.m = mediator ~ predictor + covariates
model2 <- glm(group3 ~ SCORESUM_PC1_std + sex, 
							data = df, family="binomial")

# Model.y = outcome ~ mediator + predictor + covariates
model3 <- glm(crwz5 ~ group3 + SCORESUM_PC1_std + sex,
							data = df, family="binomial")

# Mediation analysis using bootstrap & simulations to estimate CI
# Number of simulations
ndraws <- 1000

# "treat" = treatment, i.e. exposure or predictor
# "mediator" = mediator

medout_bootstrap <- mediate(model.m = model2, model.y = model3, boot = T, 
														treat = "SCORESUM_PC1_std", mediator = "group3", 
														sims = ndraws)

summ_bootstrap <- summary(medout_bootstrap)

# What you will be interested in:
# Mediation effects and pvalues (ab)
# ACME = average causal mediation effect
summ_bootstrap$d.avg # ACME_boot = 0.039
summ_bootstrap$d.avg.p  # ACME_pval_boot = 0.086# 
# Proportion of mediation
summ_bootstrap$n.avg # = 0.22

# Might be helpful, but i dont use these because when glm is involved, the effect sizes were not easy to interpret
# Direct effects (c')
# ADE = Average direct effect
summ_bootstrap$z.avg # ADE_boot = 0.14
summ_bootstrap$z.avg.p # ADE_pval_boot = 0.00***
# Total effects (c' + a*b)
summ_bootstrap$tau.coef # 0.18
summ_bootstrap$tau.p # 0.00***

# Cras_5 ~ PC1 GRS
# ------------------------------------

df <- master_GRS.df %>%
	filter(!is.na(group3), !is.na(SCORESUM_PC1_std), !is.na(sex), !is.na(cras_5)) 

# Model.m = mediator ~ predictor + covariates
model2 <- glm(group3 ~ SCORESUM_PC1_std + sex, 
							data = df, family="binomial")

# Model.y = outcome ~ mediator + predictor + covariates
model3 <- glm(cras_5 ~ group3 + SCORESUM_PC1_std + sex,
							data = df, family="binomial")

# Mediation analysis using bootstrap & simulations to estimate CI
# Number of simulations
ndraws <- 1000

# "treat" = treatment, i.e. exposure or predictor
# "mediator" = mediator

medout_bootstrap <- mediate(model.m = model2, model.y = model3, boot = T, 
														treat = "SCORESUM_PC1_std", mediator = "group3", 
														sims = ndraws)

summ_bootstrap <- summary(medout_bootstrap)

# What you will be interested in:
# Mediation effects and pvalues (ab)
# ACME = average causal mediation effect
summ_bootstrap$d.avg # ACME_boot = 0.028
summ_bootstrap$d.avg.p  # ACME_pval_boot = 0.084# 
# Proportion of mediation
summ_bootstrap$n.avg # = 0.35

# Might be helpful, but i dont use these because when glm is involved, the effect sizes were not easy to interpret
# Direct effects (c')
# ADE = Average direct effect
summ_bootstrap$z.avg # ADE_boot = 0.053
summ_bootstrap$z.avg.p # ADE_pval_boot = 0.14
# Total effects (c' + a*b)
summ_bootstrap$tau.coef # 0081
summ_bootstrap$tau.p # 0.008**


# GRS --> allergic sensitisation  --> wheeze/asthma at age 5
# ------------------------------------------------------------------------------

# Crwz5 ~ allergic disease GRS
# ------------------------------------

master_GRS.df$group3 <- as.numeric(master_GRS.df$group == 3)

df <- master_GRS.df %>%
	filter(!is.na(group3), !is.na(SCORESUM_allergic_std), !is.na(sex), !is.na(crwz5)) 

# Model.m = mediator ~ predictor + covariates
model2 <- glm(group3 ~ SCORESUM_allergic_std + sex, 
							data = df, family="binomial")

# Model.y = outcome ~ mediator + predictor + covariates
model3 <- glm(crwz5 ~ group3 + SCORESUM_allergic_std + sex,
							data = df, family="binomial")

# Mediation analysis using bootstrap & simulations to estimate CI
# Number of simulations
ndraws <- 1000

# "treat" = treatment, i.e. exposure or predictor
# "mediator" = mediator

medout_bootstrap <- mediate(model.m = model2, model.y = model3, boot = T, 
														treat = "SCORESUM_allergic_std", mediator = "group3", 
														sims = ndraws)

summ_bootstrap <- summary(medout_bootstrap)

# What you will be interested in:
# Mediation effects and pvalues (ab)
# ACME = average causal mediation effect
summ_bootstrap$d.avg # ACME_boot = 0.043
summ_bootstrap$d.avg.p  # ACME_pval_boot = 0.074# 
# Proportion of mediation
summ_bootstrap$n.avg # = 0.40

# Might be helpful, but i dont use these because when glm is involved, the effect sizes were not easy to interpret
# Direct effects (c')
# ADE = Average direct effect
summ_bootstrap$z.avg # ADE_boot = 0.066
summ_bootstrap$z.avg.p # ADE_pval_boot = 0.068#
# Total effects (c' + a*b)
summ_bootstrap$tau.coef # 0.10
summ_bootstrap$tau.p # 0.008**


# Tertile plots
# ------------------------------------------------------------------------------

# crwz5
# ------------------------------------
df <- with(master_GRS.df, 
					 data.frame(SUBJ_ID, crwz1, crwz2, crwz3, crwz4, 
					 					 crwz5 = as.numeric(as.character(crwz5)), crwz10, 
					 					 SCORESUM_allergic_tert, SCORESUM_asthma_tert,
					 					 SCORESUM_adult_asthma_tert, SCORESUM_childhood_asthma_tert,
					 					 SCORESUM_ar_tert, SCORESUM_as_tert, SCORESUM_nar_tert,
					 					 SCORESUM_combined_tert, SCORESUM_PC1_tert)) %>%
	melt(id.vars = c("SUBJ_ID", "SCORESUM_allergic_tert", "SCORESUM_asthma_tert",
									 "SCORESUM_adult_asthma_tert", "SCORESUM_childhood_asthma_tert",
									 "SCORESUM_ar_tert", "SCORESUM_as_tert", "SCORESUM_nar_tert",
									 "SCORESUM_combined_tert", "SCORESUM_PC1_tert"),
			 variable.name = "timepoint",
			 value.name = "crwz") %>%
	mutate(timepoint = as.numeric(as.character(
		factor(timepoint, levels = c("crwz1", "crwz2", "crwz3", "crwz4", "crwz5", "crwz10"),
					 labels = c(1, 2, 3, 4, 5, 10))))) %>%
	melt(id.vars =c("SUBJ_ID", "timepoint", "crwz"),
			 variable.name = "GRS",
			 value.name = "tertile") %>%
	mutate(tertile = factor(tertile, levels=c("Low", "Med", "High")))
df2 <- data.frame(with(df, table(timepoint, crwz, GRS, tertile))) %>%
	group_by(timepoint, GRS, tertile) %>%
	mutate(sum = sum(Freq),
				 prop = Freq/sum, 
				 se = 1.96*(prop*(1- prop)/sum)^0.5,
				 lower = prop - se, 
				 upper = prop + se) %>%
	filter(crwz == 1)

p.values.df <- data.frame(GRS = factor(c("SCORESUM_allergic_tert", "SCORESUM_asthma_tert",
																	"SCORESUM_adult_asthma_tert", "SCORESUM_childhood_asthma_tert",
																	"SCORESUM_ar_tert", "SCORESUM_as_tert", "SCORESUM_nar_tert",
																	"SCORESUM_combined_tert", "SCORESUM_PC1_tert"), 
																	levels=levels(df2$GRS)),
													p.value = c(NA, NA, NA, NA, "0.025", "0.013", NA, "0.0038", "0.0038"),
													OR = c(NA, NA, NA, NA, "1.2 (1-1.5)", "1.2 (1-1.5)", NA, "1.1 (1-1.1)", 
																 "1.3 (1.1-1.5)"),
													stringsAsFactors = F)

pdf("output/2019-03-28_crwz_all_GRS_lineplot.pdf", width=20, height=5)
df2 %>% 
	ungroup %>%
	ggplot(aes(x = timepoint, y = prop)) + 
	geom_line(aes(group = tertile, color = tertile)) +
	geom_point(aes(group = tertile, color = tertile, shape = tertile)) + 
	geom_ribbon(aes(x = as.numeric(timepoint), 
									ymin = lower, ymax = upper, fill = tertile), alpha = 0.1) +
	scale_color_manual(values = divPalette.alt2(3)) + 
	scale_fill_manual(values = divPalette.alt2(3)) + 
	labs(x = "Timepoint or age", y = "Proportion with current wheeze", 
			 color = "Tertile",  fill = "Tertile", shape = "Tertile") +
	facet_grid(.~factor(GRS, labels=substr(tertile_labs, 1, nchar(tertile_labs)-12))) + 
	geom_text(data = p.values.df, 
						label = ifelse(!is.na(p.values.df$p.value), 
													 paste0("p = ", p.values.df$p.value, 
													 "\n OR = ", p.values.df$OR), 
													 ""),
						x = 3.5, y = 0.5)
dev.off()

# cras
# ------------------------------------
df <- with(master_GRS.df, 
					 data.frame(SUBJ_ID, cras_3, cras_4, cras_5, cras_10,  
					 					 SCORESUM_allergic_tert, SCORESUM_asthma_tert,
					 					 SCORESUM_adult_asthma_tert, SCORESUM_childhood_asthma_tert,
					 					 SCORESUM_ar_tert, SCORESUM_as_tert, SCORESUM_nar_tert,
					 					 SCORESUM_combined_tert, SCORESUM_PC1_tert)) %>%
	melt(id.vars = c("SUBJ_ID", "SCORESUM_allergic_tert", "SCORESUM_asthma_tert",
									 "SCORESUM_adult_asthma_tert", "SCORESUM_childhood_asthma_tert",
									 "SCORESUM_ar_tert", "SCORESUM_as_tert", "SCORESUM_nar_tert",
									 "SCORESUM_combined_tert", "SCORESUM_PC1_tert"),
			 variable.name = "timepoint",
			 value.name = "cras") %>%
	mutate(timepoint = as.numeric(as.character(
		factor(timepoint, levels = c("cras_3", "cras_4", "cras_5", "cras_10"),
					 labels = c(3, 4, 5, 10))))) %>%
	melt(id.vars =c("SUBJ_ID", "timepoint", "cras"),
			 variable.name = "GRS",
			 value.name = "tertile") %>%
	mutate(tertile = factor(tertile, levels=c("Low", "Med", "High")))
df2 <- data.frame(with(df, table(timepoint, cras, GRS, tertile))) %>%
	group_by(timepoint, GRS, tertile) %>%
	mutate(sum = sum(Freq),
				 prop = Freq/sum, 
				 se = 1.96*(prop*(1- prop)/sum)^0.5,
				 lower = prop - se, 
				 upper = prop + se) %>%
	filter(cras == 1)


p.values.df <- data.frame(GRS = factor(c("SCORESUM_allergic_tert", "SCORESUM_asthma_tert",
																				 "SCORESUM_adult_asthma_tert", "SCORESUM_childhood_asthma_tert",
																				 "SCORESUM_ar_tert", "SCORESUM_as_tert", "SCORESUM_nar_tert",
																				 "SCORESUM_combined_tert", "SCORESUM_PC1_tert"), 
																			 levels=levels(df2$GRS)),
													p.value = c("0.0038", NA, NA, NA, "0.028", "0.0074", NA, "0.00038", "0.00031"),
													OR = c("1.6 (1.2-2.1)", NA, NA, NA, "1.4 (1-1.8)", "1.4 (1.1-1.9)", NA, "1.2 (1.1-1.3)", 
																 "1.7 (1.3-2.3)"),
													stringsAsFactors = F)

pdf("output/2019-03-28_cras_all_GRS_lineplot.pdf", width=20, height=5)
df2 %>% 
	ungroup %>%
	ggplot(aes(x = timepoint, y = prop)) + 
	geom_line(aes(group = tertile, color = tertile)) +
	geom_point(aes(group = tertile, color = tertile, shape = tertile)) + 
	geom_ribbon(aes(x = as.numeric(timepoint), 
									ymin = lower, ymax = upper, fill = tertile), alpha = 0.1) +
	scale_color_manual(values = divPalette.alt2(3)) + 
	scale_fill_manual(values = divPalette.alt2(3)) + 
	labs(x = "Timepoint or age", y = "Proportion with current asthma", 
			 color = "Tertile",  fill = "Tertile", shape = "Tertile") +
	facet_grid(.~factor(GRS, labels=substr(tertile_labs, 1, nchar(tertile_labs)-12))) + 
	geom_text(data = p.values.df, 
						label = ifelse(!is.na(p.values.df$p.value), 
													 paste0("p = ", p.values.df$p.value, 
													 			 "\n OR = ", p.values.df$OR), 
													 ""),
						x = 2.5, y = 0.4)
dev.off()

save.image("output/2019-03-28_GRS_associations.Rdata")
