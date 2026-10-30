# ==============================================================================
# 2018-02-06 - Batch effects on raw data
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
											"gtools")

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

load("output/2018-02-06_import_biom_format.Rdata")

# ------------------------------------------------------------------------------
# Check for batch effects
# ------------------------------------------------------------------------------

# Get ordination plot (MDS) of otu data (compositional)
# ---------------------------------------

# Get compositional abundance data
cas.biom.rel <- microbiome::transform(cas.biom, "compositional")
cas.biom.rel # 23441 taxa, 3439 samples

# Simplify data to core only

cas.biom.core <- core(cas.biom.rel, detection = 0.01/100, prevalence = 50/100)
cas.biom.core # 9 taxa, 3439 samples
cas.biom.core.rel <- microbiome::transform(cas.biom.core, "compositional")
cas.biom.core.rel.sub <- subset_samples(cas.biom.core.rel, sample_sums(cas.biom.core.rel) != 0)
cas.biom.core.rel.sub # 9 taxa, 3435 samples
which(sample_sums(cas.biom.core.rel) == 0)

set.seed(1)

# ordinate(method = MDS) is bugged.
ord.plot <- ordinate(cas.biom.core.rel.sub, method = "DCA", distance = "bray")

pdf("output/2018-02-06_dca_runname.pdf", width=16, height=10)
p <- plot_ordination(cas.biom.core.rel.sub, ord.plot, 
										 color = "RunName", 
										 title = "Run ordination")
print(p)
dev.off()

pdf("output/2018-02-06_dca_type.pdf", width=16, height=10)
p <- plot_ordination(cas.biom.core.rel.sub, ord.plot, 
										 color = "type", 
										 title = "Type ordination")
print(p)
dev.off()

pdf("output/2018-02-06_dca_age.pdf", width=16, height=10)
p <- plot_ordination(cas.biom.core.rel.sub, ord.plot, 
										 color = "age", 
										 title = "Age ordination")
print(p)
dev.off()

pdf("output/2018-02-06_dca_year.pdf", width=16, height=10)
p <- plot_ordination(cas.biom.core.rel.sub, ord.plot,
										 color = "year",
										 title = "Age (year) ordination")
print(p)
dev.off()

pdf("output/2018-02-06_dca_date_day.pdf", width=16, height=10)
p <- plot_ordination(cas.biom.core.rel.sub, ord.plot,
										 color = "date.day",
										 title = "Date of sampling ordination")
print(p)
dev.off()

pdf("output/2018-02-06_dca_date_year.pdf", width=16, height=10)
p <- plot_ordination(cas.biom.core.rel.sub, ord.plot,
										 color = "date.year",
										 title = "Date of sampling ordination")
print(p)
dev.off()

pca <- prcomp(t(otu_table(cas.biom.core.rel.sub)))
meta.df <- data.frame(sample_data(cas.biom.core.rel.sub))
pca.df <- data.frame(pca$x)

pdf("output/2018-02-06_screeplot.pdf", width=8, height=5)
screeplot(pca, type="lines")
dev.off()

pdf("output/2018-02-06_pca_runname.pdf", width=16, height=10)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$RunName))
dev.off()

pdf("output/2018-02-06_pca_type.pdf", width=16, height=10)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$type))
dev.off()

pdf("output/2018-02-06_pca_age.pdf", width=16, height=10)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$age))
dev.off()

pdf("output/2018-02-06_pca_year.pdf", width=16, height=10)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = factor(meta.df$year)))
dev.off()

pdf("output/2018-02-06_pca_date.day.pdf", width=16, height=10)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = meta.df$date.day))
dev.off()

pdf("output/2018-02-06_pca_date.year.pdf", width=16, height=10)
pca.df %>% ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color = factor(meta.df$date.year)))
dev.off()

# Test for batch effects
# ---------------------------------------
batch.effects.df <- cbind(meta.df, pca.df)

# Run name
batch.effects.df %>% ggplot(aes(y=PC1, x=RunName)) + geom_boxplot()
fit <- aov(PC1 ~ RunName, data = batch.effects.df)
summary(fit) # significant

batch.effects.df %>% ggplot(aes(y=PC2, x=RunName)) + geom_boxplot()
fit <- aov(PC2 ~ RunName, data = batch.effects.df)
summary(fit) # significant

# Age of child (day)
batch.effects.df %>% ggplot(aes(y=PC1, x=age)) + geom_point()
with(batch.effects.df, cor.test(PC1, age)) # significant but small R

batch.effects.df %>% ggplot(aes(y=PC2, x=age)) + geom_point()
with(batch.effects.df, cor.test(PC2, age)) # not significant

# Age of child (year)
batch.effects.df %>% ggplot(aes(y=PC1, x=factor(year))) + geom_boxplot()
fit <- aov(PC1 ~ year, data = batch.effects.df)
summary(fit) # significant

batch.effects.df %>% ggplot(aes(y=PC2, x=factor(year))) + geom_boxplot()
fit <- aov(PC2 ~ year, data = batch.effects.df)
summary(fit) # not significant

# Day of sampling (past first day of 1996)
batch.effects.df %>% ggplot(aes(y=PC1, x=date.day)) + geom_point()
with(batch.effects.df, cor.test(PC1, date.day)) # not significant

batch.effects.df %>% ggplot(aes(y=PC2, x=date.day)) + geom_point()
with(batch.effects.df, cor.test(PC2, date.day)) # significant

# Year of sampling (past first day of 1996)
batch.effects.df %>% ggplot(aes(y=PC1, x=factor(date.year))) + geom_boxplot()
fit <- aov(PC1 ~ date.year, data = batch.effects.df)
summary(fit) # not significant

batch.effects.df %>% ggplot(aes(y=PC2, x=factor(date.year))) + geom_boxplot()
fit <- aov(PC2 ~ date.year, data = batch.effects.df)
summary(fit) # significant

# Type (test)
batch.effects.df %>% ggplot(aes(y=PC1, x=type)) + geom_boxplot()
fit <- aov(PC1 ~ type, data = batch.effects.df)
summary(fit) # not significant

batch.effects.df %>% ggplot(aes(y=PC2, x=type)) + geom_boxplot()
fit <- aov(PC2 ~ type, data = batch.effects.df)
summary(fit) # not significant

summary(lm(PC1 ~ RunName, data = batch.effects.df)) # p < 1E-07
summary(lm(PC1 ~ type, data = batch.effects.df)) # 
summary(lm(PC1 ~ age, data = batch.effects.df)) #p = 0.001

fit <- lm(PC1 ~ RunName + age, data = batch.effects.df)
summary(fit)

summary(lm(PC2 ~ RunName, data = batch.effects.df)) # p < 1E-07
summary(lm(PC2 ~ type, data = batch.effects.df)) # 
summary(lm(PC2 ~ age, data = batch.effects.df)) #
fit <- lm(PC2 ~ RunName + age, data = batch.effects.df)
summary(fit)

# Conclusions: significant batch effects detected for Run and child age
# independent of each other

save.image("output/2018-02-06_batch_effects.Rdata")
