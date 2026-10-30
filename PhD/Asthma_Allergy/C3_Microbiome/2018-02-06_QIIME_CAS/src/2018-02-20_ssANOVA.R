# ==============================================================================
# 2018-02-20 - fitTimeSeries / ssANOVA curves
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
											"gss")

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

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

load("output/2018-02-20_micro_assoc_analysis.Rdata")

# ------------------------------------------------------------------------------
# Convert phyloseq to metagenomeseq object
# ------------------------------------------------------------------------------

# phyloseq_to_metagenomeSeq(cas.biom.method1.css) does not work - spits error
# about cumNormStat(mrobj) having an empty sample

physeq <- cas.biom.tr.cas.qc

# Narrow down to common ESVs / OTUs
physeq.common <- prune_taxa(common.otus, physeq)

# Refresh metadata to include new variables (well.sick)
sample_data(physeq.common)$well.sick <- sample_data(physeq.common)$Category
sample_data(physeq.common)$well.sick <- 
	factor(sample_data(physeq.common)$well.sick,
				 levels=c("control", "infectious"),
				 labels=c("Well", "Sick"))

physeq.common.well.sick <- 
	prune_samples(!is.na(sample_data(physeq.common)$well.sick), physeq.common)
	
# Convert to metagenomseq object
cas.metaseq <- phyloseq_to_metagenomeSeq(physeq.common.well.sick)

# function (physeq, ...) 
# {
# 	# Makes taxas rows
# 	if (!taxa_are_rows(physeq)) {
# 		physeq <- t(physeq)
# 	}
# 	# Get countdata (rounded otu_table as matrix)
# 	countData = round(as(otu_table(physeq), "matrix"), digits = 0)
# 	# Get annotated dataframe of metadata
# 	if (!is.null(sample_data(physeq, FALSE))) {
# 		ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))
# 	}
# 	else {
# 		ADF = NULL
# 	}
# 	# Get annotated dataframe of taxa data
# 	if (!is.null(tax_table(physeq, FALSE))) {
# 		TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
# 																				data.frame(tax_table(physeq)), row.names = taxa_names(physeq)))
# 	}
# 	else {
# 		TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
# 																				row.names = taxa_names(physeq)))
# 	}
# 	if (requireNamespace("metagenomeSeq")) {
# 		mrobj = metagenomeSeq::newMRexperiment(counts = countData, 
# 																					 phenoData = ADF, featureData = TDF, ...)
# 		if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
# 			p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
# 		}
# 		else {
# 			p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
# 		}
# 		mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
# 		return(mrobj)
# 	}
# }

# ------------------------------------------------------------------------------
# Run fitTimeSeries
# ------------------------------------------------------------------------------
cas.metaseq$well.sick
common.otus
length(common.otus)
mpg.otus.df <- data.frame(
	esv.id = c("dd2e15140d20c84712456e4c12ece08b",
						 "cb50e79d177c169ec5b9d5e7d308874e",
						 "f579ddd8d7c6a3d812ca33d72110f0d6",
						 "bc0d904a0469d05c65b270998488ee91",
						 "d253ca966efcd811d057dc2aa6430774",
						 "29eb8e14ce61e9c4549de0f9cf200c19",
						 "4060107ffdf52e44f9d72343b7332609",
						 "357551f644064bb2b49614b1caab2fa5"),
	esv.name = c("Alloiococcus ESV",
							 "Corynebacterium ESV",
							 "Haemophilus ESV1",
							 "Haemophilus ESV2",
							 "Moraxella ESV",
							 "Staphylococcus ESV",
							 "Streptococcus ESV1",
							 "Streptococcus ESV2"),
	MPG.name2 = c("Alloiococcus MPG",
								 "Corynebacterium MPG",
								 "Haemophilus MPG1",
								 "Haemophilus MPG2",
								 "Moraxella MPG",
								 "Staphylococcus MPG",
								 "Streptococcus MPG1",
								 "Mixed MPG"))
	
MPG.colour.scheme.df

mpg.otus.df <- merge(mpg.otus.df, MPG.colour.scheme.df, 
										 by="MPG.name2", all.x=T, all.y=F)
mpg.otus.df$esv.id <- as.character(mpg.otus.df$esv.id)
cas.fts <- list()
for (i in 1:nrow(mpg.otus.df)){
	esv.id <- mpg.otus.df$esv.id[i]
	cas.fts[[i]] <- fitTimeSeries(cas.metaseq,
																	feature = esv.id,
																	class = "well.sick",
																	id = "SUBJ_ID",
																	time="age", 
																	B=100) # B=100 for Shumei's CAS analysis
}

# ------------------------------------------------------------------------------
# Graph ssANOVA curve
# ------------------------------------------------------------------------------
# OTU names

# Get colors for each OTU curve

mpg.otus.df$hex2[which(mpg.otus.df$esv.name == "Streptococcus ESV2")] <- "#ff7f00"

pdf("output/2018-02-05_ssANOVA.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), las=1)

esv.col <- mpg.otus.df$hex2

esv.names <- as.character(mpg.otus.df$esv.name)

for(i in 1:8){
	
	# Get SSAnova fitTimeSeries object for each OTU
	res = cas.fts[[i]]
	# Get color
	col.shaded = esv.col[i]
	# Get fit, se and timepoint data (same length = 1546)
	fit = res$fit$fit
	se = res$fit$se
	timePoints = res$fit$timePoints
	confInt95 = 1.96
	# Retrieve intervals of curve for which fold change of particular OTU is significantly
	# different from 1 (log2(1) = 0).
	sigDiff = data.frame(res$timeIntervals)
	
	# Draw fitted points for OTU
	# Also, Determine plot parameters for first instance (OTU)
	if(i == 1){plot(x = timePoints, y = fit, col = "grey", xaxt = "n", ylim = c(-4,8), ylab = "log2 fold change", xlab = "Age")
		names.h = c("2m", "6m","1yr","1.5yr","2yr","2.5yr","3yr", "3.5yr", "4yr", "4.5yr", "5yr")
		Axis(side=1, at = c(60,0.5*365,365,1.5*365,2*365,2.5*365,3*365,3.5*365,4*365,4.5*365,5*365), labels = names.h)
	}
	if(i > 1)points(x = timePoints, y = fit, col = "grey")
	
	# Plot each significant curve interval 
	for (j in 1:nrow(sigDiff)) {
		begin = sigDiff[j, 1]
		end = sigDiff[j, 2]
		indBegin = which(timePoints == begin)
		indEnd = which(timePoints == end)
		x = timePoints[indBegin:indEnd]
		y = fit[indBegin:indEnd]
		xx = c(x, rev(x))
		yy = c(y, rep(0, length(y)))
		points(x,y, col=col.shaded)
		
	}
	#lines(x = timePoints, y = fit, pch = "", col = "grey")
	
	# Draw 95% CI dashed lines/curves for OTU
	lines(x = timePoints, y = fit + (confInt95 * se), pch = "", 
				lty = 2, col = col.shaded)
	lines(x = timePoints, y = fit - (confInt95 * se), pch = "", 
				lty = 2, col = col.shaded)
	
}
abline(h=0, lwd = 2)

# legend("bottomright", col = esv.col, legend = esv.names, 
#			 lwd = 3, lty = 1, bty = "n")

dev.off()


pdf("output/2018-02-05_ssANOVA_alt.pdf")
par(mar=c(5.1, 4.1, 4.1, 12.1), las=1)

esv.col <- mpg.otus.df$hex2
esv.names <- as.character(mpg.otus.df$esv.name)

for(i in 1:8){
	
	# Get SSAnova fitTimeSeries object for each OTU
	res = cas.fts[[i]]
	# Get color
	col.shaded = esv.col[i]
	# Get fit, se and timepoint data (same length = 1546)
	fit = res$fit$fit
	se = res$fit$se
	timePoints = res$fit$timePoints
	confInt95 = 1.96
	# Retrieve intervals of curve for which fold change of particular OTU is significantly
	# different from 1 (log2(1) = 0).
	sigDiff = data.frame(res$timeIntervals)
	
	# Draw fitted points for OTU
	# Also, Determine plot parameters for first instance (OTU)
	if(i == 1){plot(x = timePoints, y = fit, col = "grey", xaxt = "n", ylim = c(-4,8), ylab = "log2 fold change", xlab = "Age")
		names.h = c("2m", "6m","1yr","1.5yr","2yr","2.5yr","3yr", "3.5yr", "4yr", "4.5yr", "5yr")
		Axis(side=1, at = c(60,0.5*365,365,1.5*365,2*365,2.5*365,3*365,3.5*365,4*365,4.5*365,5*365), labels = names.h)
	}
	if(i > 1)points(x = timePoints, y = fit, col = "grey")
	
	# Plot each significant curve interval 
	for (j in 1:nrow(sigDiff)) {
		begin = sigDiff[j, 1]
		end = sigDiff[j, 2]
		indBegin = which(timePoints == begin)
		indEnd = which(timePoints == end)
		x = timePoints[indBegin:indEnd]
		y = fit[indBegin:indEnd]
		xx = c(x, rev(x))
		yy = c(y, rep(0, length(y)))
		points(x,y, col=col.shaded)
		
	}
	#lines(x = timePoints, y = fit, pch = "", col = "grey")
	
	# Draw 95% CI dashed lines/curves for OTU
	lines(x = timePoints, y = fit + (confInt95 * se), pch = "", 
				lty = 2, col = col.shaded)
	lines(x = timePoints, y = fit - (confInt95 * se), pch = "", 
				lty = 2, col = col.shaded)
	
}
abline(h=0, lwd = 2)

par(xpd=TRUE)

legend(1900,-1, col = esv.col, legend = esv.names, 
			 lwd = 3, lty = 1, bty = "n")

dev.off()

save.image("output/2018-02-05_ssANOVA.Rdata")
