# =============================================================================
# 2018-08-07 - GRS analysis
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

# npem cluster data
load("output/2018-08-07_GRS_weight_analysis.Rdata")


allergic_strong_GRS.df <- read.table("data/strong_GRS/allergic.strong.profile",
																			header=T, stringsAsFactors=F)
allergic_strong_GRS.df <- allergic_strong_GRS.df %>%
	rename(SUBJ_ID = FID,
				 SCORESUM_allergic_strong = SCORESUM,
				 CNT_allergic_strong = CNT,
				 CNT2_allergic_strong = CNT2) %>%
	mutate(IID = NULL,
				 PHENO = NULL)

asthma_strong_GRS.df <- read.table("data/strong_GRS/asthma.strong.profile",
														header=T, stringsAsFactors=F)
asthma_strong_GRS.df <- asthma_strong_GRS.df %>%
	rename(SUBJ_ID = FID,
				 SCORESUM_asthma_strong = SCORESUM,
				 CNT_asthma_strong = CNT,
				 CNT2_asthma_strong = CNT2) %>%
	mutate(IID = NULL,
				 PHENO = NULL)

# Merge with CAS Metadata
master_GRS.df <- merge(master_GRS.df, asthma_strong_GRS.df, by="SUBJ_ID", all=T)
master_GRS.df <- merge(master_GRS.df, allergic_strong_GRS.df, by="SUBJ_ID", all=T)

# Standardised scores
master_GRS.df$SCORESUM_asthma_strong_std <- as.numeric(scale(master_GRS.df$SCORESUM_asthma_strong))
master_GRS.df$SCORESUM_allergic_strong_std <- as.numeric(scale(master_GRS.df$SCORESUM_allergic_strong))

# Number of named alleles used in scoring per subject
boxplot(master_GRS.df$CNT2_asthma_strong)
boxplot(master_GRS.df$CNT2_allergic_strong)
# No obvious outliers

# ------------------------------------------------------------------------------
# Correlation of scores
# ------------------------------------------------------------------------------

lm.out <- lm(SCORESUM_allergic_strong_std ~ SCORESUM_asthma_strong_std,
						 data = master_GRS.df)
summary(lm.out) # beta = 0.88, p = <2e-16*** (EXPECTED.)
b <- good.signif(lm.out$coefficients[2], 2)
p <- good.signif(summary(lm.out)[[4]][2,4], 2)

with(master_GRS.df, cor.test(SCORESUM_asthma_strong_std, 
														 SCORESUM_allergic_strong_std, method="pearson"))
# p = 2.2e-16*** (EXPECTED.)

pdf("output/2018-08-07_correlation_of_GRS_strong.pdf", width=8, height=8)
master_GRS.df %>%
	ggplot(aes(x = SCORESUM_asthma_strong_std, y = SCORESUM_allergic_strong_std)) + 
	geom_point(aes(color = group)) + 
	geom_smooth(method = 'lm', formula = y~x) +
	labs(x = "Strong asthma GRS (standardised)", y = "Strong allergic disease GRS (standardised)") +
	scale_color_manual(name = "npEM clusters", values = quaPalette.npem[1:3], na.value = "black")
dev.off()

# Get averaged strong scores
master_GRS.df <- master_GRS.df %>%
	mutate(strong.score.avg = (SCORESUM_asthma_strong_std + SCORESUM_allergic_strong_std)/2)

# ------------------------------------------------------------------------------
# Comparison by and association with npEM cluster
# ------------------------------------------------------------------------------

# Boxplots of scores vs cluster
# ------------------------------------------------------------------------------

summary(aov(strong.score.avg ~ group, data = master_GRS.df))
# p = 0.121

t.test(x=master_GRS.df$strong.score.avg[master_GRS.df$group != 3],
			 y=master_GRS.df$strong.score.avg[master_GRS.df$group == 3])
# p = 0.03141*

pdf("output/2018-08-07_cluster_strong_GRS_boxplot.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(group) & !is.na(strong.score.avg)) %>%
	ggplot(aes(x = group, y = strong.score.avg)) + 
	geom_boxplot(aes(color = group)) + 
	scale_color_manual(values = quaPalette.npem) + 
	labs(x = "npEM cluster", y = "Averaged GRS score from strongly collinear SNPs", color = "npEM cluster") 
dev.off()


# Violins, Means and 95% CIs per cluster
# ------------------------------------------------------------------------------
df1 <- master_GRS.df %>%
	filter(!is.na(group)) %>%
	group_by(group) %>%
	summarize(mean_score_strong = mean(strong.score.avg, na.rm=T),
						se_score_strong = sd(strong.score.avg,na.rm=T)/sqrt(n()))
df2 <- master_GRS.df %>%
	filter(!is.na(group) & !is.na(strong.score.avg))

pdf("output/2018-08-07_cluster_strong_GRS_violin_mean_95CI.pdf", width=8, height=5)
ggplot() + 
	geom_violin(data = df2, 
							aes(x = group,
									y = strong.score.avg, 
									color = group),
							width = 0.5) + 
	labs(x = "npEM cluster", y = "Averaged GRS score from strongly collinear SNPs", color = "npEM cluster") +
	geom_errorbar(data = df1, 
								aes(x = group, 
										ymin=mean_score_strong - 1.96*se_score_strong, 
										max=mean_score_strong + 1.96*se_score_strong,
										color = group), width=0.2) +
	geom_point(data = df1, 
						 aes(x = group, 
						 		y=mean_score_strong,
						 		color = group)) +
	scale_color_manual(values = quaPalette.npem)
dev.off()


# ------------------------------------------------------------------------------
# Comparison by and association with npEM asthma diagnosis
# ------------------------------------------------------------------------------

# Boxplots of scores vs asthma diagnosis
# ------------------------------------------------------------------------------

t.test(x=master_GRS.df$strong.score.avg[master_GRS.df$cras_5 != 1],
			 y=master_GRS.df$strong.score.avg[master_GRS.df$cras_5 == 1])
# p = 0.08177

t.test(x=master_GRS.df$strong.score.avg[master_GRS.df$crwz5 != 1],
			 y=master_GRS.df$strong.score.avg[master_GRS.df$crwz5 == 1])
# p = 0.206


pdf("output/2018-08-07_asthma_diagnosis_strong_GRS_boxplot.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(cras_5) & !is.na(strong.score.avg)) %>%
	ggplot(aes(x = factor(cras_5), y = strong.score.avg)) + 
	geom_boxplot(aes(color = factor(cras_5))) + 
	scale_color_manual(values = c("black", "red")) + 
	labs(x = "Asthma diagnosis at age 5", y = "Strong GRS", color = "Asthma diagnosis at age 5")
dev.off()

pdf("output/2018-08-07_current_wheeze_strong_GRS_boxplot.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(crwz5) & !is.na(strong.score.avg)) %>%
	ggplot(aes(x = factor(crwz5), y = strong.score.avg)) + 
	geom_boxplot(aes(color = factor(crwz5))) + 
	scale_color_manual(values = c("black", "red")) + 
	labs(x = "Current wheeze at age 5", y = "Strong GRS", color = "Current wheeze at age 5")
dev.off()

# Violins, Means and 95% CIs per asthma diagnosis
# ------------------------------------------------------------------------------
df1 <- master_GRS.df %>%
	filter(!is.na(cras_5)) %>%
	mutate(cras_5 = factor(cras_5)) %>%
	group_by(cras_5) %>%
	summarize(mean_score_strong = mean(strong.score.avg, na.rm=T),
						se_score_strong = sd(strong.score.avg,na.rm=T)/sqrt(n()))
df2 <- master_GRS.df %>%
	filter(!is.na(cras_5) & !is.na(strong.score.avg)) %>%
	mutate(cras_5 = factor(cras_5))

pdf("output/2018-08-07_asthma_diagnosis_strong_GRS_violin_mean_95CI.pdf", width=8, height=5)
ggplot() + 
	geom_violin(data = df2, 
							aes(x = cras_5,
									y = strong.score.avg, 
									color = cras_5),
							width = 0.5) + 
	labs(x = "Asthma diagnosis at age 5", y = "Strong GRS", color = "Asthma diagnosis at age 5") +
	geom_errorbar(data = df1, 
								aes(x = cras_5, 
										ymin=mean_score_strong - 1.96*se_score_strong, 
										max=mean_score_strong + 1.96*se_score_strong,
										color = cras_5), width=0.2) +
	geom_point(data = df1, 
						 aes(x = cras_5, 
						 		y=mean_score_strong,
						 		color = cras_5)) +
	scale_color_manual(values = c("black", "red"))
dev.off()

df1 <- master_GRS.df %>%
	filter(!is.na(crwz5)) %>%
	mutate(crwz5 = factor(crwz5)) %>%
	group_by(crwz5) %>%
	summarize(mean_score_strong = mean(strong.score.avg, na.rm=T),
						se_score_strong = sd(strong.score.avg,na.rm=T)/sqrt(n()))
df2 <- master_GRS.df %>%
	filter(!is.na(crwz5) & !is.na(strong.score.avg)) %>%
	mutate(crwz5 = factor(crwz5))

pdf("output/2018-08-07_current_wheeze_strong_GRS_violin_mean_95CI.pdf", width=8, height=5)
ggplot() + 
	geom_violin(data = df2, 
							aes(x = crwz5,
									y = strong.score.avg, 
									color = crwz5),
							width = 0.5) + 
	labs(x = "Wheeze at age 5", y = "Strong GRS", color = "Wheeze at age 5") +
	geom_errorbar(data = df1, 
								aes(x = crwz5, 
										ymin=mean_score_strong - 1.96*se_score_strong, 
										max=mean_score_strong + 1.96*se_score_strong,
										color = crwz5), width=0.2) +
	geom_point(data = df1, 
						 aes(x = crwz5, 
						 		y=mean_score_strong,
						 		color = crwz5)) +
	scale_color_manual(values = c("black", "red"))
dev.off()

# Violins, Means and 95% CIs per asthma and cluster
# ------------------------------------------------------------------------------
df1 <- master_GRS.df %>%
	filter(!is.na(cras_5), !is.na(group)) %>%
	mutate(cras_5 = factor(cras_5)) %>%
	group_by(cras_5, group) %>%
	summarize(mean_score_strong = mean(strong.score.avg, na.rm=T),
						se_score_strong = sd(strong.score.avg,na.rm=T)/sqrt(n()))
df2 <- master_GRS.df %>%
	filter(!is.na(cras_5), !is.na(group), !is.na(strong.score.avg)) %>%
	mutate(cras_5 = factor(cras_5))

pdf("output/2018-08-07_asthma_diagnosis_cluster_strong_GRS_violin_mean_95CI.pdf", width=8, height=5)
ggplot() + 
	geom_violin(data = df2, 
							aes(x = cras_5,
									y = strong.score.avg, 
									color = cras_5),
							width = 0.5) + 
	labs(x = "Asthma diagnosis at age 5", y = "Asthma GRS", color = "Asthma diagnosis at age 5") +
	geom_errorbar(data = df1, 
								aes(x = cras_5, 
										ymin=mean_score_strong - 1.96*se_score_strong, 
										max=mean_score_strong + 1.96*se_score_strong,
										color = cras_5), width=0.2) +
	geom_point(data = df1, 
						 aes(x = cras_5, 
						 		y=mean_score_strong,
						 		color = cras_5)) +
	scale_color_manual(values = c("black", "red")) +
	facet_grid(.~group)
dev.off()

# ------------------------------------------------------------------------------
# Histograms with multiple groups
# ------------------------------------------------------------------------------

# npEM clusters
# ------------------------------------------------------------------------------
bw <- with(master_GRS.df, 2 * IQR(strong.score.avg, na.rm=T) / sum(!is.na(strong.score.avg))^(1/3))

pdf("output/2018-08-07_histogram_strong_GRS_cluster.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(group)) %>%
	ggplot(aes(x = strong.score.avg, fill=group)) + 
	geom_histogram(binwidth = bw, alpha=0.5, position="identity") + 
	scale_fill_manual(values = quaPalette.npem) + 
	labs(x = "Strong GRS", y = "Count", fill = "npEM cluster")
dev.off()

pdf("output/2018-08-07_density_strong_GRS_cluster.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(group)) %>% 
	ggplot(aes(x = strong.score.avg, fill=group)) + 
	geom_density(alpha=0.5) + 
	scale_fill_manual(values = quaPalette.npem) +
	labs(x = "Strong GRS", y = "Density", fill = "npEM cluster")
dev.off()

# Asthma diagnosis
# ------------------------------------------------------------------------------

pdf("output/2018-08-07_histogram_strong_GRS_asthma_diagnosis.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(cras_5)) %>%
	mutate(cras_5 = factor(cras_5)) %>%
	ggplot(aes(x = strong.score.avg, fill=cras_5)) + 
	geom_histogram(binwidth = bw, alpha=0.5, position="identity") + 
	scale_fill_manual(values = c("black", "red")) + 
	labs(x = "Strong GRS", y = "Count", fill = "Asthma diagnosis at age 5")
dev.off()

pdf("output/2018-08-07_density_strong_GRS_asthma_diagnosis.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(cras_5)) %>%
	mutate(cras_5 = factor(cras_5)) %>%
	ggplot(aes(x = strong.score.avg, fill=cras_5)) + 
	geom_density(alpha=0.5) + 
	scale_fill_manual(values = c("black", "red")) + 
	labs(x = "Strong GRS", y = "Density", fill = "Asthma diagnosis at age 5")
dev.off()

# Current wheeze
# ------------------------------------------------------------------------------

pdf("output/2018-08-07_histogram_strong_GRS_current_wheeze.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(crwz5)) %>%
	mutate(crwz5 = factor(crwz5)) %>%
	ggplot(aes(x = strong.score.avg, fill=crwz5)) + 
	geom_histogram(binwidth = bw, alpha=0.5, position="identity") + 
	scale_fill_manual(values = c("black", "red")) + 
	labs(x = "Strong GRS", y = "Count", fill = "Current wheeze at age 5")
dev.off()

pdf("output/2018-08-07_density_strong_GRS_current_wheeze.pdf", width=8, height=5)
master_GRS.df %>% 
	filter(!is.na(crwz5)) %>%
	mutate(crwz5 = factor(crwz5)) %>%
	ggplot(aes(x = strong.score.avg, fill=crwz5)) + 
	geom_density(alpha=0.5) + 
	scale_fill_manual(values = c("black", "red")) + 
	labs(x = "Strong GRS", y = "Density", fill = "Current wheeze at age 5")
dev.off()


# ------------------------------------------------------------------------------
# Association with IgE
# ------------------------------------------------------------------------------

# Total IgE
# ------------------------------------------------------------------------------

df <- with(master_GRS.df, 
					 data.frame(SUBJ_ID, totige6m, totige1, totige2, totige3, totige4, totige5, 
					 					 totige10, strong.score.avg, group) %>%
					 	melt(id.vars = c("SUBJ_ID", "strong.score.avg", "group"),
					 			 variable.name = "timepoint",
					 			 value.name = "totige")) %>%
	mutate(timepoint = as.numeric(as.character(
		factor(timepoint, levels = c("totige6m", "totige1", "totige2", 
																 "totige3", "totige4", "totige5", "totige10"),
					 labels = c(0.5, 1, 2, 3, 4, 5, 10)))))

pdf("output/2018-07-25_correlation_strong_GRS_totige.pdf", width=20, height=8)
df %>%
	filter(!is.na(group)) %>%
	ggplot(aes(x = strong.score.avg, y = totige)) + 
	geom_point(aes(color = group)) + 
	geom_smooth(method = 'lm', formula = y~x) +
	labs(y = "log10(Total IgE)", x = "Strong GRS") + 
	scale_color_manual(name = "npEM Cluster", values = quaPalette.npem) + 
	facet_grid(factor(group != 1, labels=c("Cluster 1", "Cluster 2 or 3"))~timepoint)
dev.off()

gee.out <-  tryCatch(gee(totige ~ strong.score.avg + timepoint,
												 id = SUBJ_ID,
												 data = df,
												 family = "gaussian",
												 corstr = "unstructured"),
										 error = function(err) NA)
if(is.na(gee.out)) {
	s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
} else{
	s  <- summary(gee.out)
}
beta_exact <- s$coefficients[,1]
lci_exact <- s$coefficients[,1] - 1.96*(s$coefficients[,4])
uci_exact <- s$coefficients[,1] + 1.96*(s$coefficients[,4])
p_exact <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
beta  <- good.signif(beta_exact, 2)
lci <- good.signif(lci_exact, 2)
uci <- good.signif(uci_exact, 2)
p <- signif(p_exact, 2)

results.totige.strong.grs.df <- data.frame(beta, lci, uci, p)
# beta = 0.13, p = 4.0e-06***

gee.out <-  tryCatch(gee(totige ~ strong.score.avg + timepoint,
												 id = SUBJ_ID,
												 data = df[which(df$group == 1),],
												 family = "gaussian",
												 corstr = "unstructured"),
										 error = function(err) NA)
if(is.na(gee.out)) {
	s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
} else{
	s  <- summary(gee.out)
}
beta_exact <- s$coefficients[,1]
lci_exact <- s$coefficients[,1] - 1.96*(s$coefficients[,4])
uci_exact <- s$coefficients[,1] + 1.96*(s$coefficients[,4])
p_exact <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
beta  <- good.signif(beta_exact, 2)
lci <- good.signif(lci_exact, 2)
uci <- good.signif(uci_exact, 2)
p <- signif(p_exact, 2)

results.totige.strong.grs.cluster1.df <- data.frame(beta, lci, uci, p)
# beta = -0.036, p = 4.6e-01

gee.out <-  tryCatch(gee(totige ~ strong.score.avg + timepoint,
												 id = SUBJ_ID,
												 data = df[which(df$group != 1),],
												 family = "gaussian",
												 corstr = "unstructured"),
										 error = function(err) NA)
if(is.na(gee.out)) {
	s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
} else{
	s  <- summary(gee.out)
}
beta_exact <- s$coefficients[,1]
lci_exact <- s$coefficients[,1] - 1.96*(s$coefficients[,4])
uci_exact <- s$coefficients[,1] + 1.96*(s$coefficients[,4])
p_exact <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
beta  <- good.signif(beta_exact, 2)
lci <- good.signif(lci_exact, 2)
uci <- good.signif(uci_exact, 2)
p <- signif(p_exact, 2)

results.totige.strong.grs.cluster23.df <- data.frame(beta, lci, uci, p)
# beta = 0.24, p = 2.7e-14***

# Cluster as covariate

gee.out <-  tryCatch(gee(totige ~ strong.score.avg + timepoint + group,
												 id = SUBJ_ID,
												 data = df,
												 family = "gaussian",
												 corstr = "unstructured"),
										 error = function(err) NA)
if(is.na(gee.out)) {
	s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
} else{
	s  <- summary(gee.out)
}
beta_exact <- s$coefficients[,1]
lci_exact <- s$coefficients[,1] - 1.96*(s$coefficients[,4])
uci_exact <- s$coefficients[,1] + 1.96*(s$coefficients[,4])
p_exact <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
beta  <- good.signif(beta_exact, 2)
lci <- good.signif(lci_exact, 2)
uci <- good.signif(uci_exact, 2)
p <- signif(p_exact, 2)

results.totige.strong.grs.cluster.covar.df <- data.frame(beta, lci, uci, p)
# GRS: beta = 0.1, p = 2.8e-04***
# group2: beta = 0.3, p = 1.4e-06***
# group3: beta = 1.2, p = 2.9e-63***

gee.out <-  tryCatch(gee(totige ~ strong.score.avg * timepoint * group,
												 id = SUBJ_ID,
												 data = df,
												 family = "gaussian",
												 corstr = "unstructured"),
										 error = function(err) NA)
if(is.na(gee.out)) {
	s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
} else{
	s  <- summary(gee.out)
}
beta_exact <- s$coefficients[,1]
lci_exact <- s$coefficients[,1] - 1.96*(s$coefficients[,4])
uci_exact <- s$coefficients[,1] + 1.96*(s$coefficients[,4])
p_exact <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
beta  <- good.signif(beta_exact, 2)
lci <- good.signif(lci_exact, 2)
uci <- good.signif(uci_exact, 2)
p <- signif(p_exact, 2)

results.totige.strong.grs.cluster.interact.df <- data.frame(beta, lci, uci, p)
# GRS: beta = -0.01, p = 9.1e-01
# group2: beta = 0.33, p = 2.3e-03**
# group3: beta = 1.4, p = 1.9e-26***
# timepoint x group3: beta = -0.059, p = 4.6e-02*



# ------------------------------------------------------------------------------
# Asthma ~ GRS + cluster
# ------------------------------------------------------------------------------

# GLMs
# ------------------------------------------------------------------------------

# Simple models - GRS only
glm.out <- glm(cras_5 ~ strong.score.avg, data = master_GRS.df)
summary(glm.out) # p = 0.0546

glm.out <- glm(as.numeric(as.character(crwz5))~ strong.score.avg, data = master_GRS.df)
summary(glm.out) # p = 0.19

# GRS + group
glm.out <- glm(cras_5 ~ strong.score.avg + group, data = master_GRS.df)
summary(glm.out) 
# GRS: beta = 0.03388, p = 0.24434
# group3: beta = 0.37607, p = 6.67e-05***

glm.out <- glm(cras_5 ~ SCORESUM_allergic_std + group, data = master_GRS.df)
summary(glm.out)
# GRS: beta = 0.04516, p = 0.106
# group3: beta = 0.36261, p = 0.000138***

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_asthma_std + group, data = master_GRS.df)
summary(glm.out) 
# GRS: beta = 0.02928, p = 0.348
# group3: beta = 0.523, p = 1.91e-06***

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_allergic_std + group, data = master_GRS.df)
summary(glm.out) 
# GRS: beta = 0.05821, p = 0.070
# group3: beta = 0.49803, p = 6.02e-06***

# GRS * group
glm.out <- glm(cras_5 ~ SCORESUM_asthma_std * group, data = master_GRS.df)
summary(glm.out)
# Group3: beta = 0.3035, p = 0.001685**
# GRS x group3: beta = 0.19036, p = 0.046301*

glm.out <- glm(cras_5 ~ SCORESUM_allergic_std * group, data = master_GRS.df)
summary(glm.out) 
# Group3: beta = 0.33602, p = 0.001293**

glm.out <- glm(as.numeric(as.character(crwz5))  ~ SCORESUM_asthma_std * group, data = master_GRS.df)
summary(glm.out) 
# GRS: beta = 0.09532, p = 0.0684
# Group3: beta = 0.53755, p = 3.35e-06***

glm.out <- glm(as.numeric(as.character(crwz5))  ~ SCORESUM_allergic_std * group, data = master_GRS.df)
summary(glm.out) 
# Group3: beta = 0.49049, p = 5.03e-05***

# GLMs - cluster specific
# ------------------------------------------------------------------------------

glm.out <- glm(cras_5 ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 1),])
summary(glm.out) 

glm.out <- glm(cras_5 ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 2),])
summary(glm.out) 

glm.out <- glm(cras_5 ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 3),])
summary(glm.out) # p = 0.063172

glm.out <- glm(cras_5 ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group != 1),])
summary(glm.out) # p = 0.0407*

glm.out <- glm(cras_5 ~ SCORESUM_allergic_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 1),])
summary(glm.out) 

glm.out <- glm(cras_5 ~ SCORESUM_allergic_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 2),])
summary(glm.out) 

glm.out <- glm(cras_5 ~ SCORESUM_allergic_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 3),])
summary(glm.out)

glm.out <- glm(cras_5 ~ SCORESUM_allergic_std, 
							 data = master_GRS.df[which(master_GRS.df$group != 1),])
summary(glm.out) # p = 0.0567

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 1),])
summary(glm.out) # p = 0.0841

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 2),])
summary(glm.out) 

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 3),])
summary(glm.out) 

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_asthma_std, 
							 data = master_GRS.df[which(master_GRS.df$group != 1),])
summary(glm.out) 

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_allergic_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 1),])
summary(glm.out) 

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_allergic_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 2),])
summary(glm.out) # p = 0.0819

glm.out <- glm(as.numeric(as.character(crwz5)) ~ SCORESUM_allergic_std, 
							 data = master_GRS.df[which(master_GRS.df$group == 3),])