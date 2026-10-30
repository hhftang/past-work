# ==============================================================================
# 2018-02-16 - Microbiome global ecosystem analyses
# ==============================================================================

rm(list=ls())

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------

# Bioconductor-specific packages
# ---------------------------------------
desired.biocon.packages <- c("biomformat",
														 "microbiome", 
														 "metagenomeSeq")

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
desired.packages <- c("dplyr", 
											"ggplot2", 
											"GMD", 
											"ecodist",
											"gtools",
											"cluster",
											"gee",
											"reshape2", 
											"psych")

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
if(!("metagMisc" %in% installed.packages()[,"Package"])){
	library(devtools)
	install_github("vmikk/metagMisc")
}
library(metagMisc)

# Ggplot and graphical settings
# ---------------------------------------
theme_set(theme_bw())

seqPalette <- colorRampPalette(c("#fef0d9","#fdd49e","#fdbb84",
																 "#fc8d59","#e34a33","#b30000"))
quaPalette <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
								'#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
# quaPalette <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
# 								'#ff7f00','#ffff33','#a65628','#f781bf',
# 								'#999999')
# quaPalette <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4',
#								'#fed9a6','#ffffcc','#e5d8bd','#fddaec',
#								'#f2f2f2')
divPalette <- colorRampPalette(rev(c("#d73027","#fc8d59","#fee090",
																		 "#e0f3f8","#91bfdb","#4575b4")))

# Functions
# ---------------------------------------
# A better signif
good.signif <- function(x, n) {
	trimws(formatC(signif(x,digits=n), digits=n,format="fg"))
}
# Note that this converts number to character string

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("output/2018-02-20_micro_assoc_analysis.Rdata")

# ------------------------------------------------------------------------------
# Alpha diversity
# ------------------------------------------------------------------------------

diversities.df <- diversities(coast.biom.tr.sum.coast.qc2, index = "all")

# Relative abundance = p

# True diversity = effective number of taxa; number of equally abundant taxa
# needed for average rel abundance to equal observed average rel abundance.

# qD = (sum(p^q))^(1/(1-q)) 
# where sum is sum across all taxa, and q is some non-negative integer

# q = 0:
# 0D = sum(1) = R where R is richness or number of taxa (past detection threshold)

# q = 1: undefined but appraoches
# 1D = exp( - sum(p * ln(p))) (exponential of Shannon diversity)

# q = 2:
# 2D = 1/sum(p^2) = Inverse simpson index

# as q approaches infinity, D approaches largest p in sample

# inverse_simpson = 2D = 1/lambda where lambda = sum(p^2)
# gini_simpson = 1 - lambda = 1 - 1/2D
# Shannon diversity = - sum(p * ln(p)) = ln(1D)
# Fisher alpha
# Coverage (number of species needed to cover 50% ecosystem)

diversities.df$LogNo <- row.names(diversities.df)
log.table.eco.df <- merge(log.table.df, diversities.df, by="LogNo", all=T)

pdf("output/2018-02-16_inverse_simpson_all.pdf", width=8, height=5)
log.table.eco.df %>% 
	filter(!is.na(age.custom), !is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=inverse_simpson)) + 
	geom_boxplot() +
	facet_grid(.~well.sick) +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Diversity (inverse Simpson)") +
	theme(axis.text.x=element_text(angle=90, hjust=1))
dev.off()

pdf("output/2018-02-16_shannon_all.pdf", width=8, height=5)
log.table.eco.df %>% 
	filter(!is.na(age.custom), !is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=shannon)) + 
	geom_boxplot() +
	facet_grid(.~well.sick) +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Diversity (Shannon)") +
	theme(axis.text.x=element_text(angle=90, hjust=1))
dev.off()

# GEE for diversity vs. age
gee.out <- gee(inverse_simpson ~  
							 	age.year + gender + well.sick,
							 id = subj, 
							 data = log.table.eco.df, 
							 family= gaussian,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
					 p = p)

gee.out <- gee(inverse_simpson ~  
							 	age.day + gender + well.sick,
							 id = subj, 
							 data = log.table.eco.df, 
							 family= gaussian,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
					 p = p)

gee.out <- gee(shannon ~  
							 	age.year + gender + well.sick,
							 id = subj, 
							 data = log.table.eco.df, 
							 family= gaussian,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
					 p = p)

gee.out <- gee(shannon ~  
							 	age.day + gender + well.sick,
							 id = subj, 
							 data = log.table.eco.df, 
							 family= gaussian,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
					 p = p)

# Plot diversities ~ age (timepoints), facet by MPG, well vs. unwell
pdf("output/2018-02-16_inverse_simpson_mpg.pdf", width=15, height=6)
log.table.eco.df %>% 
	filter(!is.na(age.custom), !is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=inverse_simpson, color=MPG.name.top)) + 
	geom_boxplot() +
	scale_color_manual(values = levels(log.table.df$MPG.colour.hex.top)) + 
	facet_grid(well.sick~factor(MPG.name.top, levels=rev(levels(MPG.name.top))), scales="free_y", space="free_y") +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Diversity (inverse Simpson)") +
	theme(axis.text.x=element_text(angle=90, hjust=1))
dev.off()

pdf("output/2018-02-16_shannon_mpg.pdf", width=15, height=6)
log.table.eco.df %>% 
	filter(!is.na(age.custom), !is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, y=shannon, color=MPG.name.top)) + 
	geom_boxplot() +
	scale_color_manual(values = levels(log.table.df$MPG.colour.hex.top)) + 
	facet_grid(well.sick~factor(MPG.name.top, levels=rev(levels(MPG.name.top))), scales="free_y", space="free_y") +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Diversity (Shannon)") +
	theme(axis.text.x=element_text(angle=90, hjust=1))
dev.off()

# GEE for diversity vs. MPG + age

MPG.list <- levels(log.table.eco.df$MPG.name.top)
results <- list(length(MPG.list))

for(i in 1:length(MPG.list)){
	gee.out <- gee(inverse_simpson ~  
							 	I(MPG.name.top == MPG.list[i]) * age.day + gender + well.sick,
							 id = subj, 
							 data = log.table.eco.df, 
							 family= gaussian,
							 corstr = "unstructured")
	s  <- summary(gee.out)
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results[[i]] <- data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
						 p = p)
	row.names(results[[i]])[2] <- MPG.list[i]
}
mpg.diversity.list.inverse.simpson <- results
mpg.diversity.inverse.simpson.df <- rbind.fill(mpg.diversity.list.inverse.simpson[[1]][2,],
											mpg.diversity.list.inverse.simpson[[2]][2,],
											mpg.diversity.list.inverse.simpson[[3]][2,],
											mpg.diversity.list.inverse.simpson[[4]][2,],
											mpg.diversity.list.inverse.simpson[[5]][2,],
											mpg.diversity.list.inverse.simpson[[6]][2,],
											mpg.diversity.list.inverse.simpson[[7]][2,],
											mpg.diversity.list.inverse.simpson[[8]][2,],
											mpg.diversity.list.inverse.simpson[[9]][2,])
mpg.diversity.inverse.simpson.df$MPG.names <- MPG.list
mpg.diversity.inverse.simpson.df <- 
	mpg.diversity.inverse.simpson.df[match(rev(MPG.list),
																				 mpg.diversity.inverse.simpson.df$MPG.names),]
write.csv(mpg.diversity.inverse.simpson.df, "output/2018-02-16_mpg_diversity_inverse_simpson.csv", 
					quote=F, row.names=F)

# GEE for diversity vs. MPG + age

MPG.list <- levels(log.table.eco.df$MPG.name.top)
results <- list(length(MPG.list))

for(i in 1:length(MPG.list)){
	gee.out <- gee(shannon ~  
								 	I(MPG.name.top == MPG.list[i]) * age.day + gender + well.sick,
								 id = subj, 
								 data = log.table.eco.df, 
								 family= gaussian,
								 corstr = "unstructured")
	s  <- summary(gee.out)
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results[[i]] <- data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
														 p = p)
	row.names(results[[i]])[2] <- MPG.list[i]
}
mpg.diversity.list.shannon <- results
mpg.diversity.shannon.df <- rbind.fill(mpg.diversity.list.shannon[[1]][2,],
																							 mpg.diversity.list.shannon[[2]][2,],
																							 mpg.diversity.list.shannon[[3]][2,],
																							 mpg.diversity.list.shannon[[4]][2,],
																							 mpg.diversity.list.shannon[[5]][2,],
																							 mpg.diversity.list.shannon[[6]][2,],
																							 mpg.diversity.list.shannon[[7]][2,],
																							 mpg.diversity.list.shannon[[8]][2,],
																							 mpg.diversity.list.shannon[[9]][2,])
mpg.diversity.shannon.df$MPG.names <- MPG.list
mpg.diversity.shannon.df <- 
	mpg.diversity.shannon.df[match(rev(MPG.list),
																				 mpg.diversity.shannon.df$MPG.names),]
write.csv(mpg.diversity.shannon.df, "output/2018-02-16_mpg_diversity_shannon.csv", 
					quote=F, row.names=F)


gee.out <- gee(inverse_simpson ~  
							 	MPG.name.top + age.day + gender + well.sick,
							 id = subj, 
							 data = log.table.eco.df, 
							 family= gaussian,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
					 p = p)

gee.out <- gee(shannon ~  
							 	MPG.name.top + age.day + gender + well.sick,
							 id = subj, 
							 data = log.table.eco.df, 
							 family= gaussian,
							 corstr = "unstructured")
s  <- summary(gee.out)
OR  <-  exp(s$coefficients[,1])
lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
OR <- good.signif(OR,2)
lci <- good.signif(lci,2)
uci <- good.signif(uci,2)
p <- signif(p,2)
data.frame(OR = paste0(OR, " (", lci, "-", uci, ")"), 
					 p = p)

# ------------------------------------------------------------------------------
# Richness
# ------------------------------------------------------------------------------

richness(coast.biom.tr.sum.coast.qc2, detection = 0)[1,]
sum(otu_table(coast.biom.tr.sum.coast.qc2)[,1] > 0)
length(which(otu_table(coast.biom.tr.sum.coast.qc2)[,1] > 0))

# Check with varying detection thresholds
richness(coast.biom.tr.sum.coast.qc2)

# Gives richness (number of taxa) above different thresholds

# ------------------------------------------------------------------------------
# Dominance
# ------------------------------------------------------------------------------
# 'absolute' This is the most simple variant, giving the absolute abundance of 
# the most abundant species (Magurran & McGill 2011). By default, this refers
# to the single most dominant species (rank=1) but it is possible to calculate 
# the absolute dominance with rank n based on the abundances of top-n species 
# by tuning the rank argument.

# 'relative' Relative abundance of the most abundant species. 
# This is with rank=1 by default but can be calculated for other ranks.

# 'DBP' Berger–Parker index, a special case of relative dominance with rank 1; 
# This also equals the inverse of true diversity of the infinite order.

# 'DMN' McNaughton’s dominance. This is the sum of the relative abundance of 
# the two most abundant taxa, or a special case of relative dominance with rank 2

# 'simpson' Simpson's index ($sum(p^2)$) where p are relative abundances has 
# an interpretation as a dominance measure. Also the version 
# ($sum(q * (q-1)) / S(S-1)$) based on absolute abundances q has been proposed 
# by Simpson (1949) but not included here as it is not within [0,1] range, and 
# it is highly correlated with the simpler Simpson dominance. Finally, it is 
# also possible to calculate dominances up to an arbitrary rank by
# setting the rank argument

# 'core_abundance' Relative proportion of the core species that exceed 
# detection level 0.2% in over 50% of the samples

# 'gini' Gini index is calculated with the function inequality.

# dominance(coast.biom.tr.sum.coast.qc2)

