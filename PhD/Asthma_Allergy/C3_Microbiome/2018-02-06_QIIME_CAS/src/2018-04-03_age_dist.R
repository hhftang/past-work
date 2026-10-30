# ==============================================================================
# 2018-04-03 - Age sample distributions
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

load("output/2018-02-27_core_common_microbiome.Rdata")

# URI vs. LRI
with(log.table.df, table(Category, lri, useNA="ifany"))
with(log.table.df, table(Category, uri, useNA="ifany"))
log.table.df$uri.lri <- ifelse(log.table.df$Category == "control", "Well", 
															 ifelse(log.table.df$lri == 1, "LRI", 
															 			 ifelse(log.table.df$uri == 1, "URI", NA)))
log.table.df$uri.lri <- factor(log.table.df$uri.lri, levels=c("Well", "URI", "LRI"))
with(log.table.df, table(Category, uri.lri, useNA="ifany"))

pdf("output/2018-04-03_age_dist.pdf", width=10, height=5)
log.table.df %>%
	filter(!is.na(uri.lri)) %>%
	ggplot(aes(x=age)) + geom_bar(aes(fill=uri.lri)) + 
	scale_fill_manual(values=c("blue", "red", "brown")) +
	facet_grid(uri.lri~., scales="free_y", space="free_y") +
	scale_x_continuous(name="Age",
										 breaks=c(0, 365/12*2, 365/12*6, 365, 
										 				 365*1.5, 365*2, 365*2.5, 365*3, 365*3.5, 365*4,
										 				 365*4.5, 365*5),
										 labels=c("Birth", "2mth", "6mth","1y", 
										 				 "1.5y", "2y", "2.5y", "3y", "3.5y", "4y", 
										 				 "4.5y", "5y")) + 
	scale_y_continuous(name="Count") + 
	theme(panel.grid.minor = element_blank())
dev.off()
# 
# df <- spec.subj.table.df %>%
# 	filter(type == "CAS", !is.na(Category), age < 365*4.75) %>%
# 	mutate(age.custom = cut(as.numeric(age), 
# 													breaks=c(0, 365/12*4, 365/12*9, 
# 																	 365*1.25, 365*1.75, 365*2.25, 365*2.75, 365*3.25,
# 																	 365*3.75, 365*4.25),
# 													labels=c("2mth", "6mth", "1y", 
# 																	 "1.5y", "2y", "2.5y", "3y", "3.5y", "4y"))) 
# with(df, table(Category, age.custom))
# df2 <- data.frame(with(df, prop.table(table(Category, age.custom),2)))
# 
# pdf("output/2018-02-06_prop_routine.pdf", width=8, height=2)
# df2 %>%
# 	mutate(Type = Category) %>%
# 	ggplot(aes(x=age.custom, y=Freq, group=Type, color=Type)) + geom_point() + 
# 	geom_line() + 
# 	scale_color_manual(values=c("blue", "red")) +
# 	scale_x_discrete(name="Timepoint of collection") + 
# 	scale_y_continuous(name="Proportion of \n samples") + 
# 	theme(panel.grid.minor = element_blank())
# dev.off()
