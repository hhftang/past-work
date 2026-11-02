# =============================================================================
# 2019-04-05 - GRS ~ in vitro gene expression associations
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
											"NMF", 
											"MatrixEQTL",
											"ClusterR")

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

# GRS data
load("output/2019-04-05_GRS_prep.Rdata")

# Get GRS columns
# ------------------------------------------------------------------------------
sample_grs.df <- master_GRS.df %>%
	dplyr::select(SUBJ_ID, one_of(grs_phenotype_std)) %>%
	rename(subjid = SUBJ_ID)
sample_grs.t.df <- as.data.frame(base::t(sample_grs.df))
colnames(sample_grs.t.df) <- sample_grs.df$subjid

# Output GRS data as external file 
# (Columns = samples, rows = GRS)
write.table(sample_grs.t.df, "output/grs_all.txt", quote=F, col.names=F)

# For each GRS phenotype in grs_phenotype, create separate output
for (i in 1:length(grs_phenotype)) {
	target <- grs_phenotype[i]
	target_row <- which(rownames(sample_grs.t.df) == paste0("SCORESUM_", target, "_std"))
	target_grs.df <- sample_grs.t.df[c(1,target_row),]
	write.table(target_grs.df, paste0("output/grs_", target, ".txt"), quote=F, col.names=F)
}
# 
# # Get covariates
# # ------------------------------------------------------------------------------
# # Read 10 genotype PCs and gender
# genoPCs <- read.table("data/qinqin_cell_expression/PCs_10_cas_michigan_imp_filtered_135indiv.txt", 
# 											header = T, stringsAsFactors = F)
# gender <- read.table("data/qinqin_cell_expression/Sample_QC_information.txt", 
# 										 header = T, stringsAsFactors = F) %>%
# 	dplyr::select(subjid, sex_fe1ma2) %>%
# 	unique()
# if(identical(genoPCs$subjid, gender$sex_fe1ma2)) {
# 	genoPCs$gender <- gender$sex_fe1ma2    # 2 = male, 1 = female
# }
# cov <- as.data.frame(t(genoPCs))
# write.table(cov, "output/Cov_gender_10_genoPCs.txt", col.names = F, quote = F, sep = "\t")
# 
# # Load covariates for all 135 individuals
# cvrt = SlicedData$new()
# cvrt$fileDelimiter = "\t"
# cvrt$fileOmitCharacters = "NA"
# cvrt$fileSkipRows = 1
# cvrt$fileSkipColumns = 1
# cvrt$LoadFile("output/Cov_gender_10_genoPCs.txt")
# #cvrt$columnNames
# # 
# # Do association analyses
# # ------------------------------------------------------------------------------
# # Loop for GRSs
# for(whichGRS in grs_phenotype) {
# 	# Create a folder for this GRS
# 	cat("    GRS for", whichGRS, "\n")
# 	dirname <- paste0("~/CAS/Analysis/GRS_expression/", gsub(" ", "_" ,whichGRS))
# 	if(!file.exists(dirname)) {dir.create(dirname)}
# 	setwd(dirname)
# 	
# 	# Read GRS
# 	GRS <- read.table(CASGRS$filename[which(CASGRS$GRS == whichGRS)], header = T)
# 	
# 	# Scale
# 	GRS$scaled_GRS <- scale(GRS$SCORESUM)
# 	t_GRS <- as.data.frame(t(GRS[,c("IID","scaled_GRS")]))
# 	write.table(t_GRS, "GRS_scaled.txt", quote = F, sep = "\t", col.names = F)
# 	
# 	# Load in Matrix eQTL
# 	GRS = SlicedData$new()
# 	GRS$fileDelimiter = "\t"
# 	GRS$fileOmitCharacters = "NA"
# 	GRS$fileSkipRows = 1
# 	GRS$fileSkipColumns = 1
# 	GRS$fileSliceSize = 2000   # read file in pieces of 2000 rows
# 	GRS$LoadFile("GRS_scaled.txt")
# 	
# 	# Four conditions
# 	foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA")) %do% {
# 		cat("  * Working on", celltype, treatment, "...\n")
# 		
# 		# Load gene expression data
# 		gene = SlicedData$new()
# 		gene$fileDelimiter = "\t"
# 		gene$fileOmitCharacters = "NA"
# 		gene$fileSkipRows = 1
# 		gene$fileSkipColumns = 1
# 		gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
# 		gene$LoadFile(paste0("~/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/gene_expression_",celltype,"_",treatment,".txt"))
# 		
# 		# Subject names (column names) 
# 		expr_colname <- sapply(gene$columnNames, function(x) {unlist(strsplit(x, split = "_"))[2]})
# 		
# 		# Keep overlapping individuals
# 		# Subset GRS
# 		GRS_curr <- GRS$Clone()
# 		GRS_curr$ColumnSubsample(match(expr_colname, GRS$columnNames))
# 		# Subest covariates
# 		cvrt_curr <- cvrt$Clone()
# 		cvrt_curr$ColumnSubsample(match(expr_colname, cvrt$columnNames))
# 		
# 		if(!all(identical(GRS_curr$columnNames, as.vector(expr_colname)), identical(cvrt_curr$columnNames, as.vector(expr_colname)))) {
# 			cat("  Error - Check column names!\n")
# 			return(NULL)
# 		}
# 		
# 		#----- Run MatrixeQTL on all genes -----
# 		me <- Matrix_eQTL_main(
# 			snps = GRS_curr,
# 			gene = gene,
# 			cvrt = cvrt_curr,
# 			pvOutputThreshold = 1,
# 			output_file_name = NULL,
# 			useModel = modelLINEAR,
# 			errorCovariance = numeric(),
# 			verbose = FALSE,
# 			min.pv.by.genesnp = FALSE,
# 			noFDRsaveMemory = FALSE)
# 		
# 		# Save association test results
# 		write.table(me$all$eqtls, paste0("Geneexpr_GRS_",celltype,"_",treatment,".txt"), quote = F, sep = "\t", row.names = F)
# 		
# 		
# 		#----- Run MatrixeQTL on cis-eGenes -----
# 		# Keep genes that had significant cis eQTLs
# 		gene$RowReorder(which(gene$GetAllRowNames() %in% ciseGene[[paste0(celltype,"_",treatment)]]))
# 		
# 		me <- Matrix_eQTL_main(
# 			snps = GRS_curr,
# 			gene = gene,
# 			cvrt = cvrt_curr,
# 			pvOutputThreshold = 1,
# 			output_file_name = NULL,
# 			useModel = modelLINEAR,
# 			errorCovariance = numeric(),
# 			verbose = FALSE,
# 			min.pv.by.genesnp = FALSE,
# 			noFDRsaveMemory = FALSE)
# 		
# 		# Save association test results
# 		write.table(me$all$eqtls, paste0("CiseGene_Geneexpr_GRS_",celltype,"_",treatment,".txt"), quote = F, sep = "\t", row.names = F)
# 		
# 	} # End of 4 conditions
# 	
# } # End of all GRSs

# ------------------------------------------------------------------------------
# Get Qinqin's results
# ------------------------------------------------------------------------------

sig_gene_fdr_with_grs.df <- read.csv("data/qinqin_cell_expression/sig_fdr_withGRS.csv",
																			 row.names = NULL, stringsAsFactors=F)

# ------------------------------------------------------------------------------
# Load expression data separately
# ------------------------------------------------------------------------------

# Rank normalised
# ------------------------------------------------------------------------------

cell_conditions <- c("m_Ctrl", "m_LPS", "t_Ctrl", "t_PHA")

for (i in 1:length(cell_conditions)){
	x <- data.frame(t(read.table(paste0("data/qinqin_cell_expression/Rank_normal_transformation_gene_expression_",
				 															 cell_conditions[i], ".txt"), 
				 															 header=F)), stringsAsFactors = F)
	colnames(x) <- x[1,]
	colnames(x)[1] <- "SampleID"
	x <- x[-1,]
	assign(paste0("rank_norm_", cell_conditions[i], ".df"), x)
}

# Quantile normalised
# ------------------------------------------------------------------------------

cell_count <- c(130, 141, 142, 144)

for (i in 1:length(cell_conditions)){
	x <- data.frame(t(read.table(paste0("data/qinqin_cell_expression/QuanNorm_alltogether_gene_expression_",
																			cell_conditions[i], "_", cell_count[i], ".txt"), 
															 header=F)), stringsAsFactors = F)
	colnames(x) <- x[1,]
	colnames(x)[1] <- "SampleID"
	x <- x[-1,]
	assign(paste0("quan_norm_", cell_conditions[i], ".df"), x)
}

# Sample info
# ------------------------------------------------------------------------------
sample.info.df <- read.table("data/qinqin_cell_expression/Sample_QC_information.txt", 
										 header = T, stringsAsFactors = F)

# ------------------------------------------------------------------------------
# Battery of tests for all significant genes
# ------------------------------------------------------------------------------

# List of genes
# ------------------------------------------------------------------------------
gene_of_interest_m_LPS <- unique(sig_gene_fdr_with_grs.df$gene[which(sig_gene_fdr_with_grs.df$condition == "LPS M")])
gene_of_interest_t_Ctrl <- unique(sig_gene_fdr_with_grs.df$gene[which(sig_gene_fdr_with_grs.df$condition == "Resting T")])
gene_of_interest_t_PHA <- unique(sig_gene_fdr_with_grs.df$gene[which(sig_gene_fdr_with_grs.df$condition == "PHA T")])

# Phenotypes
# ------------------------------------------------------------------------------
# wheeze, asthma, sensitisation, eczema, hdmige, pntige, cluster 3, lung function
# ari, lri, flri, wlri frequency
# use premaster.df for ae2, master_all.df for rest


phenotypes_bin <- c("crwz5", "cras_5", "ae2", "crez2", "crez5", 
										"hdmige2 >= log10(0.35)", "pntige6m >= -1.25",
										"group == 3")
phenotypes_cont <- c("hdmige2", "pntige6m", "riskMPGs_0_2",
										 "uri_1", "uri_2", "uri_3",
										 "lri_1", "lri_2", "lri_3",
										 "wlri_1", "wlri_2", "wlri_3",
										 "flri_1", "flri_2", "flri_3",
										 "FVC.pre.y5", "FEV1.pre.y5", "FEV1.FVC.pre.y5", "PEF.pre.y5", 
										 "FVC.pre.y10", "FEV1.pre.y10", "FEV1.FVC.pre.y10", "FEF2575.pre.y10")
phenotypes <- c(phenotypes_bin, phenotypes_cont)

# Merge dataframes
# ------------------------------------------------------------------------------
master_all.df <- master_all.df %>%
	merge((premaster.df %>% 
				 	dplyr::select(SUBJ_ID, ae6m, ae1, ae2, 
				 								ae3, ae4, ae5, ae10)), 
				by="SUBJ_ID", all=T)

# Prepare output lists
# ------------------------------------------------------------------------------
results.m_LPS.q.list <- list()
results.m_LPS.r.list <- list()
results.t_Ctrl.q.list <- list()
results.t_Ctrl.r.list <- list()
results.t_PHA.q.list <- list()
results.t_PHA.r.list <- list()

# ------------------------------------------------------------------------------
# Loop
# ------------------------------------------------------------------------------

# m_LPS
# ------------------------------------------------------------------------------
for (i in 1:length(gene_of_interest_m_LPS)) {
	gene_target <- gene_of_interest_m_LPS[i]
	results.q.df <- data.frame(gene = gene_target,
														 phenotype = character(length(phenotypes)), 
													 effect_size = numeric(length(phenotypes)),
													 lci = numeric(length(phenotypes)),
													 uci = numeric(length(phenotypes)),
													 p_value = numeric(length(phenotypes)),
													 stringsAsFactors = F)
	results.r.df <- data.frame(gene = gene_target,
														 phenotype = character(length(phenotypes)), 
														 effect_size = numeric(length(phenotypes)),
														 lci = numeric(length(phenotypes)),
														 uci = numeric(length(phenotypes)),
														 p_value = numeric(length(phenotypes)),
														 stringsAsFactors = F)
	df.q <- quan_norm_t_PHA.df %>%
		dplyr::select(SampleID, gene_target) %>%
		mutate(target_q = as.numeric(get(gene_target)))
	df.r <- rank_norm_t_PHA.df %>%
		dplyr::select(SampleID, gene_target) %>%
		mutate(target_r = as.numeric(get(gene_target)))
	df <- merge(df.q, df.r, by=c("SampleID"))
	df <- sample.info.df %>% 
		dplyr::select(SampleID, subjid) %>%
		rename(SUBJ_ID = subjid) %>%
		merge(df, by="SampleID", all=T) %>%
		filter(!is.na(target_q), !is.na(target_r))
	df <- merge(master_all.df, df, by="SUBJ_ID", all=T)
	df <- merge(df, lungfunction.y5.df, by="SUBJ_ID", all=T)
	df <- merge(df, lungfunction.y10.df, by="SUBJ_ID", all=T)
	
	# Binary traits:
	# ------------------------------------------------------------------------------
	for (j in 1:length(phenotypes_bin)) {
		phenotype_target <- phenotypes_bin[j]
		
		# quantile normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_q + sex"), 
									 family = "binomial",
									 data = df), 
									 error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.q.df$phenotype[j] <- phenotype_target
		results.q.df$effect_size[j] <- s$coefficients[2,1]
		results.q.df$lci[j] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.q.df$uci[j] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.q.df$p_value[j] <-  s$coefficients[2,4]
		
		# Rank-normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_r + sex"), 
														family = "binomial",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.r.df$phenotype[j] <- phenotype_target
		results.r.df$effect_size[j] <- s$coefficients[2,1]
		results.r.df$lci[j] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.r.df$uci[j] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.r.df$p_value[j] <-  s$coefficients[2,4]
	}
	
	# Continuous traits:
	# ------------------------------------------------------------------------------
	for (j in 1:length(phenotypes_cont)) {
		k <- j + length(phenotypes_bin)
		phenotype_target <- phenotypes_cont[j]
		
		# quantile normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_q + sex"), 
														family = "gaussian",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.q.df$phenotype[k] <- phenotype_target
		results.q.df$effect_size[k] <- s$coefficients[2,1]
		results.q.df$lci[k] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.q.df$uci[k] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.q.df$p_value[k] <-  s$coefficients[2,4]
		
		# Rank-normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_r + sex"), 
														family = "gaussian",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.r.df$phenotype[k] <- phenotype_target
		results.r.df$effect_size[k] <- s$coefficients[2,1]
		results.r.df$lci[k] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.r.df$uci[k] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.r.df$p_value[k] <-  s$coefficients[2,4]
	}
	results.q.df$p_fdr <- p.adjust(results.q.df$p_value, method="fdr")
	results.r.df$p_fdr <- p.adjust(results.r.df$p_value, method="fdr")
	results.m_LPS.q.list[[i]] <- results.q.df
	results.m_LPS.r.list[[i]] <- results.r.df
	print(paste0("Finished gene: ", gene_target))
}

results.m_LPS.q.all.df <- do.call("rbind", results.m_LPS.q.list)
results.m_LPS.r.all.df <- do.call("rbind", results.m_LPS.r.list)

# p_fdr: FDR calculated within gene only.
# p_fdr_all: FDR calculated across all comparisons (1519)
results.m_LPS.q.all.df$p_fdr_all <- p.adjust(results.m_LPS.q.all.df$p_value, method="fdr")
results.m_LPS.r.all.df$p_fdr_all <- p.adjust(results.m_LPS.r.all.df$p_value, method="fdr")

results.m_LPS.q.all.df[which(results.m_LPS.q.all.df$p_value < 0.05),]
results.m_LPS.q.all.df[which(results.m_LPS.q.all.df$p_fdr < 0.05),]
results.m_LPS.q.all.df[which(results.m_LPS.q.all.df$p_fdr_all < 0.1),]

results.m_LPS.r.all.df[which(results.m_LPS.r.all.df$p_value < 0.05),]
results.m_LPS.r.all.df[which(results.m_LPS.r.all.df$p_fdr < 0.05),]
results.m_LPS.r.all.df[which(results.m_LPS.r.all.df$p_fdr_all < 0.1),]

# t_PHA
# ------------------------------------------------------------------------------

for (i in 1:length(gene_of_interest_t_PHA)) {
	gene_target <- gene_of_interest_t_PHA[i]
	results.q.df <- data.frame(gene = gene_target,
														 phenotype = character(length(phenotypes)), 
														 effect_size = numeric(length(phenotypes)),
														 lci = numeric(length(phenotypes)),
														 uci = numeric(length(phenotypes)),
														 p_value = numeric(length(phenotypes)),
														 stringsAsFactors = F)
	results.r.df <- data.frame(gene = gene_target,
														 phenotype = character(length(phenotypes)), 
														 effect_size = numeric(length(phenotypes)),
														 lci = numeric(length(phenotypes)),
														 uci = numeric(length(phenotypes)),
														 p_value = numeric(length(phenotypes)),
														 stringsAsFactors = F)
	df.q <- quan_norm_t_PHA.df %>%
		dplyr::select(SampleID, gene_target) %>%
		mutate(target_q = as.numeric(get(gene_target)))
	df.r <- rank_norm_t_PHA.df %>%
		dplyr::select(SampleID, gene_target) %>%
		mutate(target_r = as.numeric(get(gene_target)))
	df <- merge(df.q, df.r, by=c("SampleID"))
	df <- sample.info.df %>% 
		dplyr::select(SampleID, subjid) %>%
		rename(SUBJ_ID = subjid) %>%
		merge(df, by="SampleID", all=T) %>%
		filter(!is.na(target_q), !is.na(target_r))
	df <- merge(master_all.df, df, by="SUBJ_ID", all=T)
	df <- merge(df, lungfunction.y5.df, by="SUBJ_ID", all=T)
	df <- merge(df, lungfunction.y10.df, by="SUBJ_ID", all=T)
	
	# Binary traits:
	# ------------------------------------------------------------------------------
	for (j in 1:length(phenotypes_bin)) {
		phenotype_target <- phenotypes_bin[j]
		
		# quantile normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_q + sex"), 
														family = "binomial",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.q.df$phenotype[j] <- phenotype_target
		results.q.df$effect_size[j] <- s$coefficients[2,1]
		results.q.df$lci[j] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.q.df$uci[j] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.q.df$p_value[j] <-  s$coefficients[2,4]
		
		# Rank-normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_r + sex"), 
														family = "binomial",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.r.df$phenotype[j] <- phenotype_target
		results.r.df$effect_size[j] <- s$coefficients[2,1]
		results.r.df$lci[j] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.r.df$uci[j] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.r.df$p_value[j] <-  s$coefficients[2,4]
	}
	
	# Continuous traits:
	# ------------------------------------------------------------------------------
	for (j in 1:length(phenotypes_cont)) {
		k <- j + length(phenotypes_bin)
		phenotype_target <- phenotypes_cont[j]
		
		# quantile normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_q + sex"), 
														family = "gaussian",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.q.df$phenotype[k] <- phenotype_target
		results.q.df$effect_size[k] <- s$coefficients[2,1]
		results.q.df$lci[k] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.q.df$uci[k] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.q.df$p_value[k] <-  s$coefficients[2,4]
		
		# Rank-normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_r + sex"), 
														family = "gaussian",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.r.df$phenotype[k] <- phenotype_target
		results.r.df$effect_size[k] <- s$coefficients[2,1]
		results.r.df$lci[k] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.r.df$uci[k] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.r.df$p_value[k] <-  s$coefficients[2,4]
	}
	results.q.df$p_fdr <- p.adjust(results.q.df$p_value, method="fdr")
	results.r.df$p_fdr <- p.adjust(results.r.df$p_value, method="fdr")
	results.t_PHA.q.list[[i]] <- results.q.df
	results.t_PHA.r.list[[i]] <- results.r.df
	print(paste0("Finished gene: ", gene_target))
}

results.t_PHA.q.all.df <- do.call("rbind", results.t_PHA.q.list)
results.t_PHA.r.all.df <- do.call("rbind", results.t_PHA.r.list)

# p_fdr: FDR calculated within gene only.
# p_fdr_all: FDR calculated across all comparisons (1519)
results.t_PHA.q.all.df$p_fdr_all <- p.adjust(results.t_PHA.q.all.df$p_value, method="fdr")
results.t_PHA.r.all.df$p_fdr_all <- p.adjust(results.t_PHA.r.all.df$p_value, method="fdr")

results.t_PHA.q.all.df[which(results.t_PHA.q.all.df$p_value < 0.05),]
results.t_PHA.q.all.df[which(results.t_PHA.q.all.df$p_fdr < 0.05),]
results.t_PHA.q.all.df[which(results.t_PHA.q.all.df$p_fdr_all < 0.1),]

results.t_PHA.r.all.df[which(results.t_PHA.r.all.df$p_value < 0.05),]
results.t_PHA.r.all.df[which(results.t_PHA.r.all.df$p_fdr < 0.05),]
results.t_PHA.r.all.df[which(results.t_PHA.r.all.df$p_fdr_all < 0.1),]

# t_Ctrl
# ------------------------------------------------------------------------------

for (i in 1:length(gene_of_interest_t_Ctrl)) {
	gene_target <- gene_of_interest_t_Ctrl[i]
	results.q.df <- data.frame(gene = gene_target,
														 phenotype = character(length(phenotypes)), 
														 effect_size = numeric(length(phenotypes)),
														 lci = numeric(length(phenotypes)),
														 uci = numeric(length(phenotypes)),
														 p_value = numeric(length(phenotypes)),
														 stringsAsFactors = F)
	results.r.df <- data.frame(gene = gene_target,
														 phenotype = character(length(phenotypes)), 
														 effect_size = numeric(length(phenotypes)),
														 lci = numeric(length(phenotypes)),
														 uci = numeric(length(phenotypes)),
														 p_value = numeric(length(phenotypes)),
														 stringsAsFactors = F)
	df.q <- quan_norm_t_Ctrl.df %>%
		dplyr::select(SampleID, gene_target) %>%
		mutate(target_q = as.numeric(get(gene_target)))
	df.r <- rank_norm_t_Ctrl.df %>%
		dplyr::select(SampleID, gene_target) %>%
		mutate(target_r = as.numeric(get(gene_target)))
	df <- merge(df.q, df.r, by=c("SampleID"))
	df <- sample.info.df %>% 
		dplyr::select(SampleID, subjid) %>%
		rename(SUBJ_ID = subjid) %>%
		merge(df, by="SampleID", all=T) %>%
		filter(!is.na(target_q), !is.na(target_r))
	df <- merge(master_all.df, df, by="SUBJ_ID", all=T)
	df <- merge(df, lungfunction.y5.df, by="SUBJ_ID", all=T)
	df <- merge(df, lungfunction.y10.df, by="SUBJ_ID", all=T)
	
	# Binary traits:
	# ------------------------------------------------------------------------------
	for (j in 1:length(phenotypes_bin)) {
		phenotype_target <- phenotypes_bin[j]
		
		# quantile normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_q + sex"), 
														family = "binomial",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.q.df$phenotype[j] <- phenotype_target
		results.q.df$effect_size[j] <- s$coefficients[2,1]
		results.q.df$lci[j] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.q.df$uci[j] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.q.df$p_value[j] <-  s$coefficients[2,4]
		
		# Rank-normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_r + sex"), 
														family = "binomial",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.r.df$phenotype[j] <- phenotype_target
		results.r.df$effect_size[j] <- s$coefficients[2,1]
		results.r.df$lci[j] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.r.df$uci[j] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.r.df$p_value[j] <-  s$coefficients[2,4]
	}
	
	# Continuous traits:
	# ------------------------------------------------------------------------------
	for (j in 1:length(phenotypes_cont)) {
		k <- j + length(phenotypes_bin)
		phenotype_target <- phenotypes_cont[j]
		
		# quantile normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_q + sex"), 
														family = "gaussian",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.q.df$phenotype[k] <- phenotype_target
		results.q.df$effect_size[k] <- s$coefficients[2,1]
		results.q.df$lci[k] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.q.df$uci[k] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.q.df$p_value[k] <-  s$coefficients[2,4]
		
		# Rank-normalised
		glm.out <- tryCatch(glm(paste0("(", phenotype_target, ") ~ target_r + sex"), 
														family = "gaussian",
														data = df), 
												error = function(err) NA)
		if(is.na(glm.out)) {
			s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
		} else{
			s  <- summary(glm.out)
		}
		results.r.df$phenotype[k] <- phenotype_target
		results.r.df$effect_size[k] <- s$coefficients[2,1]
		results.r.df$lci[k] <- s$coefficients[2,1] - 1.96*(s$coefficients[2,2])
		results.r.df$uci[k] <- s$coefficients[2,1] + 1.96*(s$coefficients[2,2])
		results.r.df$p_value[k] <-  s$coefficients[2,4]
	}
	results.q.df$p_fdr <- p.adjust(results.q.df$p_value, method="fdr")
	results.r.df$p_fdr <- p.adjust(results.r.df$p_value, method="fdr")
	results.t_Ctrl.q.list[[i]] <- results.q.df
	results.t_Ctrl.r.list[[i]] <- results.r.df
	print(paste0("Finished gene: ", gene_target))
}

results.t_Ctrl.q.all.df <- do.call("rbind", results.t_Ctrl.q.list)
results.t_Ctrl.r.all.df <- do.call("rbind", results.t_Ctrl.r.list)

# p_fdr: FDR calculated within gene only.
# p_fdr_all: FDR calculated across all comparisons (1519)
results.t_Ctrl.q.all.df$p_fdr_all <- p.adjust(results.t_Ctrl.q.all.df$p_value, method="fdr")
results.t_Ctrl.r.all.df$p_fdr_all <- p.adjust(results.t_Ctrl.r.all.df$p_value, method="fdr")

results.t_Ctrl.q.all.df[which(results.t_Ctrl.q.all.df$p_value < 0.05),]
results.t_Ctrl.q.all.df[which(results.t_Ctrl.q.all.df$p_fdr < 0.05),]
results.t_Ctrl.q.all.df[which(results.t_Ctrl.q.all.df$p_fdr_all < 0.1),]

results.t_Ctrl.r.all.df[which(results.t_Ctrl.r.all.df$p_value < 0.05),]
results.t_Ctrl.r.all.df[which(results.t_Ctrl.r.all.df$p_fdr < 0.05),]
results.t_Ctrl.r.all.df[which(results.t_Ctrl.r.all.df$p_fdr_all < 0.1),]

# ------------------------------------------------------------------------------
# Write results
# ------------------------------------------------------------------------------
write.csv(results.m_LPS.r.all.df[which(results.m_LPS.r.all.df$p_value < 0.05),],
						"output/assoc.phenotype.m_LPS.csv", quote=F, row.names=F)
write.csv(results.t_Ctrl.r.all.df[which(results.t_Ctrl.r.all.df$p_value < 0.05),],
						"output/assoc.phenotype.t_Ctrl.csv", quote=F, row.names=F)
write.csv(results.t_PHA.r.all.df[which(results.t_PHA.r.all.df$p_value < 0.05),],
						"output/assoc.phenotype.t_PHA.csv", quote=F, row.names=F)

# ------------------------------------------------------------------------------
# GRS data no NA
# ------------------------------------------------------------------------------

sample_grs_complete.df <- sample_grs.df %>%
	filter(!is.na(SCORESUM_Meta_PEF_std))
sample_grs_complete.df <- unique(sample_grs_complete.df)
rownames(sample_grs_complete.df) <- sample_grs_complete.df$subjid
sample_grs_complete.df$subjid <- NULL

sample_grs.pca <- prcomp(sample_grs_complete.df)
screeplot(sample_grs.pca, npcs=10)

sample_grs_pcs.df <- as.data.frame(sample_grs.pca$x,
																	 stringsAsFactors = F)

# Plot first PCs
plot(sample_grs_pcs.df[,1], sample_grs_pcs.df[,2])
pairs(sample_grs_pcs.df[,1:6])

# Find optimal number of clusters or components
# -----------------------------------------------------------------------------
Optimal_Clusters_KMeans(sample_grs_complete.df, max_clusters=10,
												criterion="BIC", 
												num_init=10, max_iters=200)
Optimal_Clusters_GMM(sample_grs_complete.df, max_clusters=10,
												criterion="BIC")
# 6 clusters

# kmeans
# -----------------------------------------------------------------------------
sample_grs.kmeans <-kmeans(sample_grs_complete.df, 6)

sample_grs_pcs.df$cluster <- sample_grs.kmeans$cluster
lattice::splom(sample_grs_pcs.df[,1:6],
			col=sample_grs_pcs.df$cluster)

# GMM
# -----------------------------------------------------------------------------
sample_grs.gmm <-GMM(sample_grs_complete.df, 6)
gmm.cluster <- apply(sample_grs.gmm$Log_likelihood, 1, which.min)

sample_grs_pcs.df$g.cluster <- gmm.cluster
lattice::splom(sample_grs_pcs.df[,1:6],
							 col=sample_grs_pcs.df$g.cluster)

with(sample_grs_pcs.df, table(cluster, g.cluster))
