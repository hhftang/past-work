# Compare results

rm(list=ls())

# Load libraries

desired.packages <- c("MASS", "ggplot2", "gplots", "corrplot")
new.packages <- desired.packages[!(desired.packages %in% 
																	 	installed.packages()[,"Package"])]

if(length(new.packages)>0) install.packages(new.packages, 
																						repos='https://cran.ms.unimelb.edu.au/')
for (i in 1:length(desired.packages)){
	library(desired.packages[i], character.only=T)
}

# Set RColorBrewer palettes
# Hex codes extracted using ColorBrewer website: http://colorbrewer2.org/
# Divergent palette (e.g -1 to 0 to 1)
divPalette <- colorRampPalette(rev(c("#d73027","#fc8d59","#fee090",
																		 "#e0f3f8","#91bfdb","#4575b4")))
# Sequential palette (e.g. 0 to 10)
seqPalette <- colorRampPalette(c("#fef0d9","#fdd49e","#fdbb84",
																 "#fc8d59","#e34a33","#b30000"))
# Qualitative palette (e.g. A to E)
quaPalette <- c("#8da0cb","#a6d854","#fc8d62","#e78ac3","#66c2a5","#ffd92f")

# Read results tables
longgwas.df <- read.table("output/longgwas_results.txt", header=T, 
													stringsAsFactors = F)
colnames(longgwas.df) <- c("id", "p.longgwas")

repeatabel.df <- read.table("output/repeatABEL_results.txt", header=T, 
													stringsAsFactors = F)
colnames(repeatabel.df) <- c("id", "p.repeatabel")

pgs.sis.lmm.df <- read.table("output/pgs_sis_lmm_results.txt", header=T, 
														stringsAsFactors = F)
colnames(pgs.sis.lmm.df) <- c("id", "p.sis.lmm")

pgs.sis.gee.df <- read.table("output/pgs_sis_gee_results.txt", header=T, 
														 stringsAsFactors = F)
colnames(pgs.sis.gee.df) <- c("id", "p.sis.gee")

gee.df <- read.table("output/gee_results.txt", header=T, 
														 stringsAsFactors = F)
colnames(gee.df) <- c("id", "p.gee")

results.df <- merge(longgwas.df, repeatabel.df, by="id", all=T)
results.df <- merge(results.df, pgs.sis.lmm.df, by="id", all=T)
results.df <- merge(results.df, pgs.sis.gee.df, by="id", all=T)
results.df <- merge(results.df, gee.df, by="id", all=T)
results.mat <- as.matrix(results.df[,2:6])
results.mat <- -log10(results.mat)


# Heatmap
heatmap.2(results.mat,
					scale = "none", 
					na.color = "grey",
					col = seqPalette,
					trace = "none", 
					margins = c(8,5), 
					density.info = "none",
					notecol = "black"
)

# p-value of correlations
cor.mtest <- function(mat, ...) {
	mat <- as.matrix(mat)
	n <- ncol(mat)
	p.mat<- matrix(NA, n, n)
	diag(p.mat) <- 0
	for (i in 1:(n - 1)) {
		for (j in (i + 1):n) {
			tmp <- cor.test(mat[, i], mat[, j], ...)
			p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
		}
	}
	colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
	p.mat
}

p.mat <- cor.mtest(results.mat)
p.mat
cor.mat <- cor(results.mat)
cor.mat

# Correlation plot
divPalette.alt <- colorRampPalette(
	c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(cor.mat, 
				 method="color", 
				 col=divPalette.alt(200),  
				 type="upper", order="hclust", 
				 addCoef.col = "black", # Add coefficient of correlation
				 tl.col="black", tl.srt=45, #Text label color and rotation
				 # Combine with significance
				 p.mat = p.mat, 
				 sig.level = 0.01, 
				 insig = "blank", 
				 # hide correlation coefficient on the principal diagonal
				 diag=FALSE 
)
