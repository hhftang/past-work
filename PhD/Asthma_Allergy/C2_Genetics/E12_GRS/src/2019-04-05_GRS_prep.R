# =============================================================================
# 2019-04-05 - GRS data preparation
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
											"psych", 
											"NMF")

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

save.image("output/2019-04-05_GRS_prep.Rdata")