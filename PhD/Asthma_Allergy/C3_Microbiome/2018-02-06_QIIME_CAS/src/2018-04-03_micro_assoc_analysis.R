# ==============================================================================
# 2018-04-03 - Microbiome distributions and association analyses (GEEs)
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
											"reshape2")

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
quaPalette <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
								'#ff7f00','#ffff33','#a65628','#f781bf',
								'#999999')
# quaPalette <- c('#fbb4ae','#b3cde3','#ccebc5','#decbe4',
#								'#fed9a6','#ffffcc','#e5d8bd','#fddaec',
#								'#f2f2f2')

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

load("output/2018-02-27_core_common_microbiome.Rdata")

# ------------------------------------------------------------------------------
# Define seasons
# ------------------------------------------------------------------------------

# Define by Northern Meteorological convention
# https://www.timeanddate.com/calendar/aboutseasons.html
# According to the meteorological definition, the seasons begin on the first day of the months that include the equinoxes and solstices:
# 	
# Spring runs from March 1 to May 31;
# Summer runs from June 1 to August 31;
# Fall (autumn) runs from September 1 to November 30; and
# Winter runs from December 1 to February 28 (February 29 in a leap year).

log.table.df$month <- as.numeric(format(
	as.POSIXct(as.Date(as.character(log.table.df$Date), format="%Y-%m-%d")), "%m"))
log.table.df$month.name <- 
	factor(log.table.df$month,
				 labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
				 				 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
log.table.df$season <- factor(
	ifelse(log.table.df$month %in% c(9,10,11), "Spring",
				 ifelse(log.table.df$month %in% c(12,1,2), "Summer",
				 			 ifelse(log.table.df$month %in% c(3,4,5), "Autumn",
				 			 			 ifelse(log.table.df$month %in% c(6,7,8), "Winter", NA))))
)
log.table.df$season <- factor(log.table.df$season, 
															levels = c("Winter","Spring", "Summer", "Autumn"))

log.table.df$season2 <- log.table.df$season
levels(log.table.df$season2) <- c("Winter.Spring", "Winter.Spring",
																	"Summer.Fall", "Summer.Fall")

# ------------------------------------------------------------------------------
# MPG distribution by visittype and illness type
# ------------------------------------------------------------------------------

# Calculate MPG distribution for visittype
# ---------------------------------------
# all samples
mpg.prop.df <- with(log.table.df, data.frame(table(MPG.name)))
sum(mpg.prop.df$Freq) # 3120 = correct
mpg.prop.df$Prop <- mpg.prop.df$Freq/sum(mpg.prop.df$Freq)

# For each illness type - **IGNORE MILD** for well/sick
log.table.df$well.sick <- log.table.df$Category
log.table.df$well.sick <- 
	factor(log.table.df$well.sick, 
				 levels=c("control", "infectious"),
				 
				 labels=c("Well", "Sick"))

# URI vs. LRI
log.table.df$uri.lri <- ifelse(log.table.df$Category == "control", "Well", 
															 ifelse(log.table.df$lri == 1, "LRI", 
															 			 ifelse(log.table.df$uri == 1, "URI", NA)))
log.table.df$uri.lri <- factor(log.table.df$uri.lri, levels=c("Well", "URI", "LRI"))
with(log.table.df, table(Category, uri.lri, useNA="ifany"))

mpg.count.per.type.df <- as.data.frame.matrix(
	with(log.table.df, table(MPG.name, uri.lri)))
colnames(mpg.count.per.type.df) <- c("Well.Freq","URI.Freq", "LRI.Freq")
mpg.prop.per.type.df <- as.data.frame.matrix(
	with(log.table.df, prop.table(table(MPG.name, uri.lri),2)))
colnames(mpg.prop.per.type.df) <- c("Well.Prop","URI.Prop", "LRI.Prop")

mpg.distr.df <- data.frame(MPG.name = mpg.prop.df$MPG.name,
													 Well.Freq = mpg.count.per.type.df$Well.Freq,
													 Well.Prop = paste0(good.signif(mpg.prop.per.type.df$Well.Prop*100,2), "%"),
													 URI.Freq = mpg.count.per.type.df$URI.Freq,
													 URI.Prop = paste0(good.signif(mpg.prop.per.type.df$URI.Prop*100,2), "%"), 
													 LRI.Freq = mpg.count.per.type.df$LRI.Freq,
													 LRI.Prop = paste0(good.signif(mpg.prop.per.type.df$LRI.Prop*100,2), "%"), 
													 stringsAsFactors = F)
mpg.distr.df2 <- data.frame(MPG.name = mpg.distr.df$MPG.name,
														Well = paste0(mpg.distr.df$Well.Freq, " (", 
																					mpg.distr.df$Well.Prop, ")"),
														URI = paste0(mpg.distr.df$URI.Freq, " (", 
																				 mpg.distr.df$URI.Prop, ")"),
														LRI = paste0(mpg.distr.df$LRI.Freq, " (", 
																				 mpg.distr.df$LRI.Prop, ")"))

MPG.sick.sig.order <- c("Haemophilus MPG2", "Streptococcus MPG1", "Haemophilus MPG1", 
												"Moraxella MPG", "Mixed MPG", "Corynebacterium MPG", 
												"Alloiococcus MPG", "Staphylococcus MPG", 
												"Moraxellaceae MPG1",  
												"Moraxella rare MPG","Moraxellaceae rare MPG",
												"Corynebacterium rare MPG", "Haemophilus rare MPG",
												"Streptococcus rare MPG")
MPG.alpha.order <- c("Alloiococcus MPG", "Corynebacterium MPG", "Haemophilus MPG1", 
										 "Haemophilus MPG2", "Mixed MPG", "Moraxella MPG", 
										 "Staphylococcus MPG", "Streptococcus MPG1", 
										 "Moraxellaceae MPG1", 
										 "Corynebacterium rare MPG", "Haemophilus rare MPG", 
										 "Moraxella rare MPG", "Moraxellaceae rare MPG", 
										 "Streptococcus rare MPG")

# Run GEE for LRI vs. Well ~ MPG, adjusting for gender, age, season
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- unique(mpg.distr.df2$MPG.name)
MPG.list

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 OR.exact = numeric(length(MPG.list)),
												 lci.exact = numeric(length(MPG.list)),
												 uci.exact = numeric(length(MPG.list)),
												 p.exact = numeric(length(MPG.list)),
												 stringsAsFactors=F)
df <- log.table.df %>%
	filter(uri.lri != "URI")
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(uri.lri=="LRI") ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + season2,
								 id = SUBJ_ID, 
								 data = df, 
								 family= binomial,
								 corstr = "unstructured"),
								 error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	results.df$OR.exact[i] <- OR[2]
	results.df$lci.exact[i] <- lci[2]
	results.df$uci.exact[i] <- uci[2]
	results.df$p.exact[i] <- p[2]
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.distr.df2$OR <- results.df$OR
mpg.distr.df2$p <- results.df$p

mpg.distr.df2 <- mpg.distr.df2[match(MPG.sick.sig.order, mpg.distr.df2$MPG.name),]

write.csv(mpg.distr.df2, "output/2018-04-03_mpg.distribution.results.csv", quote=F,
					row.names=F)

mpg.distr.df$OR.exact <- results.df$OR.exact
mpg.distr.df$lci.exact <- results.df$lci.exact
mpg.distr.df$uci.exact <- results.df$uci.exact
mpg.distr.df$p.exact <- results.df$p.exact

mpg.distr.df <- mpg.distr.df[match(MPG.sick.sig.order, mpg.distr.df$MPG.name),]

write.csv(mpg.distr.df, "output/2018-04-03_mpg.distribution.results.exact.csv", quote=F,
					row.names=F)
write.csv(mpg.distr.df, 
					"/Users/htang/Documents/PhD/Project - CAS/C3_Microbiome/2018-04-03_QIIME_CAS_COAST_meta/data/2018-04-03_cas.lri.well.csv", 
					quote=F,
					row.names=F)

# Run GEE for URI vs. Well ~ MPG, adjusting for gender, age, season
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- unique(mpg.distr.df2$MPG.name)
MPG.list

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 OR.exact = numeric(length(MPG.list)),
												 lci.exact = numeric(length(MPG.list)),
												 uci.exact = numeric(length(MPG.list)),
												 p.exact = numeric(length(MPG.list)),
												 stringsAsFactors=F)
df <- log.table.df %>%
	filter(uri.lri != "LRI")
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(uri.lri=="URI") ~  
														I(MPG.name == MPG.list[i]) + sex + age + season2,
													id = SUBJ_ID, 
													data = df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	results.df$OR.exact[i] <- OR[2]
	results.df$lci.exact[i] <- lci[2]
	results.df$uci.exact[i] <- uci[2]
	results.df$p.exact[i] <- p[2]
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.distr.df2$OR <- results.df$OR
mpg.distr.df2$p <- results.df$p

mpg.distr.df2 <- mpg.distr.df2[match(MPG.sick.sig.order, mpg.distr.df2$MPG.name),]

write.csv(mpg.distr.df2, "output/2018-04-03_mpg.distribution.results.uri.csv", quote=F,
					row.names=F)

mpg.distr.df$OR.exact <- results.df$OR.exact
mpg.distr.df$lci.exact <- results.df$lci.exact
mpg.distr.df$uci.exact <- results.df$uci.exact
mpg.distr.df$p.exact <- results.df$p.exact

mpg.distr.df <- mpg.distr.df[match(MPG.sick.sig.order, mpg.distr.df2$MPG.name),]

write.csv(mpg.distr.df, "output/2018-04-03_mpg.distribution.results.uri.exact.csv", quote=F,
					row.names=F)
write.csv(mpg.distr.df, 
					"/Users/htang/Documents/PhD/Project - CAS/C3_Microbiome/2018-04-03_QIIME_CAS_COAST_meta/data/2018-04-03_cas.uri.well.csv", 
					quote=F,
					row.names=F)

# Run GEE for LRI/URI vs. Well ~ MPG, adjusting for gender, age, season
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- unique(mpg.distr.df2$MPG.name)
MPG.list

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 OR.exact = numeric(length(MPG.list)),
												 lci.exact = numeric(length(MPG.list)),
												 uci.exact = numeric(length(MPG.list)),
												 p.exact = numeric(length(MPG.list)),
												 stringsAsFactors=F)
df <- log.table.df
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(uri.lri!="Well") ~  
														I(MPG.name == MPG.list[i]) + sex + age + season2,
													id = SUBJ_ID, 
													data = df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	results.df$OR.exact[i] <- OR[2]
	results.df$lci.exact[i] <- lci[2]
	results.df$uci.exact[i] <- uci[2]
	results.df$p.exact[i] <- p[2]
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.distr.df2$OR <- results.df$OR
mpg.distr.df2$p <- results.df$p

mpg.distr.df2 <- mpg.distr.df2[match(MPG.sick.sig.order, mpg.distr.df2$MPG.name),]

write.csv(mpg.distr.df2, "output/2018-04-03_mpg.distribution.results.ari.csv", quote=F,
					row.names=F)

mpg.distr.df$OR.exact <- results.df$OR.exact
mpg.distr.df$lci.exact <- results.df$lci.exact
mpg.distr.df$uci.exact <- results.df$uci.exact
mpg.distr.df$p.exact <- results.df$p.exact

mpg.distr.df <- mpg.distr.df[match(MPG.sick.sig.order, mpg.distr.df2$MPG.name),]

write.csv(mpg.distr.df, "output/2018-04-03_mpg.distribution.results.ari.exact.csv", quote=F,
					row.names=F)
write.csv(mpg.distr.df, 
					"/Users/htang/Documents/PhD/Project - CAS/C3_Microbiome/2018-04-03_QIIME_CAS_COAST_meta/data/2018-04-03_cas.ari.well.csv", 
					quote=F,
					row.names=F)

# Plot proportions of MPGs per visittype/illness type
# Breaks based on on sampling scheme (time period around timepoint
# ---------------------------------------

# Plot MPG per visittype
log.table.df$age.custom <- 
	cut(log.table.df$age, 
			breaks=c(0, 365/12*4, 365/12*9, 
							 365*1.25, 365*1.75, 365*2.25, 365*2.75, 365*3.25,
							 365*3.75, 365*4.25),
			labels=c("2mth", "6mth", "1y", 
							 "1.5y", "2y", "2.5y", "3y", "3.5y", "4y"))

# Target MPG order
log.table.df$MPG.name <- 
	factor(log.table.df$MPG.name,
				 levels = rev(c("Haemophilus MPG1", "Haemophilus MPG2", "Moraxella MPG", 
				 							 "Streptococcus MPG1", "Mixed MPG", "Staphylococcus MPG", 
				 							 "Alloiococcus MPG", "Corynebacterium MPG", 
				 							 "Corynebacterium rare MPG", "Haemophilus rare MPG",
				 							 "Moraxellaceae MPG1", "Moraxellaceae rare MPG", "Moraxella rare MPG",
				 							 "Streptococcus rare MPG")))
log.table.df$MPG.name.top <- log.table.df$MPG.name
levels(log.table.df$MPG.name.top) <- 
	rev(c("Haemophilus MPG1", "Haemophilus MPG2", "Moraxella MPG",
				"Streptococcus MPG1", "Mixed MPG", "Staphylococcus MPG",
				"Alloiococcus MPG", "Corynebacterium MPG",
				"Other MPG", "Other MPG", 
				"Other MPG", "Other MPG", "Other MPG",
				"Other MPG"))
log.table.df$MPG.colour.hex <- 
	factor(log.table.df$MPG.colour.hex,
				 levels = rev(c("#ffd700", "#b29600", "#e31a1c", 
				 							 "#fdbf6f", "#808080", "#cab2db", 
				 							 "#b2df8a", "#a6cee3", 
				 							 "#1f78b4", "#e5c100",
				 							 "#ff00ff", "#ffccff", "#fb9a99",
				 							 "#b15928")))
log.table.df$MPG.colour.hex.top <- log.table.df$MPG.colour.hex
levels(log.table.df$MPG.colour.hex.top) <- 
	rev(c("#ffd700", "#b29600", "#e31a1c", 
				"#fdbf6f", "#808080", "#cab2db", 
				"#b2df8a", "#a6cee3", 
				"#d3d3d3", "#d3d3d3", 
				"#d3d3d3", "#d3d3d3","#d3d3d3",
				"#d3d3d3"))

log.table.df <- log.table.df[order(log.table.df$MPG.name, 
																	 log.table.df$age.custom),]

pdf("output/2018-02-20_MPG_distr_visittype_block.pdf", width=8, height=5)
log.table.df %>% 
	filter(!is.na(age.custom), !is.na(well.sick)) %>%
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(log.table.df$MPG.colour.hex.top)) + 
	facet_grid(.~well.sick) +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Proportion") +
	theme(axis.text.x=element_text(angle=90, hjust=1))
dev.off()

# Plot MPG per illness type
log.table.df$illness.type <- factor(
	ifelse(log.table.df$well.sick == "Well", "Well", 
				 ifelse(log.table.df$uri == 1, "URI",
				 			 ifelse(log.table.df$lri == "1", "LRI", NA))),
	levels=c("Well", "URI", "LRI"))
log.table.df$illness.type

pdf("output/2018-02-06_MPG_distr_illness.type_block.pdf", width=8, height=5)
log.table.df %>% 
	filter(!is.na(age.custom), !is.na(illness.type)) %>%
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(log.table.df$MPG.colour.hex.top)) + 
	facet_grid(.~illness.type) +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Proportion") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) 
dev.off()

pdf("output/2018-02-06_MPG_distr_illness.type_block_2y.pdf", width=8, height=5)
log.table.df %>% 
	filter(!is.na(age.custom), 
				 !(age.custom %in% c("2.5y", "3y", "3.5y", "4y")),
				 !is.na(illness.type)) %>%
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(log.table.df$MPG.colour.hex.top)) + 
	facet_grid(.~illness.type) +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Proportion") +
	theme(axis.text.x=element_text(angle=90, hjust=1)) 
dev.off()

# Look for trend between MPGs and time, with well.sick
# ---------------------------------------
MPG.list <- MPG.alpha.order

mpg.prop.well.age.df <- 
	as.data.frame.matrix(with(log.table.df[which(log.table.df$well.sick == "Well"),], 
														prop.table(table(MPG.name, age.custom),2)))
mpg.prop.sick.age.df <- 
	as.data.frame.matrix(with(log.table.df[which(log.table.df$well.sick == "Sick"),], 
														prop.table(table(MPG.name, age.custom),2)))

results.df <- data.frame(mpg = MPG.list, 
												 rho.cor.well = character(length(MPG.list)),
												 p.cor.well = character(length(MPG.list)),
												 OR.gee.well = character(length(MPG.list)),
												 p.gee.well = character(length(MPG.list)),
												 rho.cor.sick = character(length(MPG.list)),
												 p.cor.sick = character(length(MPG.list)),
												 OR.gee.sick = character(length(MPG.list)),
												 p.gee.sick = character(length(MPG.list)),												 
												 stringsAsFactors=F)

for (i in 1:length(MPG.list)){
	# Trend
	df <- data.frame(mpg = t(mpg.prop.well.age.df[
		which(row.names(mpg.prop.well.age.df) == MPG.list[i]),]),
		time = c(2/12, 6/12, 1, 1.5, 2, 2.5, 3, 3.5, 4))
	colnames(df) <- c("mpg", "time")
	cor.results <- with(df, cor.test(mpg, time, method = "spearman"))
	rho.cor <- cor.results$estimate
	p.cor <- cor.results$p.value
	results.df$rho.cor.well[i] <- good.signif(rho.cor, 2)
	results.df$p.cor.well[i] <- good.signif(p.cor, 2) 
	# gee
	df <- log.table.df[which(log.table.df$well.sick == "Well"),]
	# df$age.custom <- 	as.numeric(
	#		as.character(factor(df$age.custom,labels=c(2/12, 4/12, 6/12, 9/12, 1, 1.5, 2))))
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														age + sex + season,
													id = SUBJ_ID_ID, 
													data = df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR.gee.well[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p.gee.well[i] <- paste0(p)[2]
}

for (i in 1:length(MPG.list)){
	# Trend
	df <- data.frame(mpg = t(mpg.prop.sick.age.df[
		which(row.names(mpg.prop.sick.age.df) == MPG.list[i]),]),
		time = c(2/12, 6/12, 1, 1.5, 2, 2.5, 3, 3.5, 4))
	colnames(df) <- c("mpg", "time")
	cor.results <- with(df, cor.test(mpg, time, method = "spearman"))
	rho.cor <- cor.results$estimate
	p.cor <- cor.results$p.value
	results.df$rho.cor.sick[i] <- good.signif(rho.cor, 2)
	results.df$p.cor.sick[i] <- good.signif(p.cor, 2) 
	# gee
	df <- log.table.df[which(log.table.df$well.sick == "Sick"),]
	# df$age.custom <- 	as.numeric(
	#		as.character(factor(df$age.custom,labels=c(2/12, 4/12, 6/12, 9/12, 1, 1.5, 2))))
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														age + sex + season,
													id = SUBJ_ID_ID, 
													data = df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR.gee.sick[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p.gee.sick[i] <- paste0(p)[2]
}

MPG.graph.order <- rev(c("Haemophilus MPG1", "Haemophilus MPG2", "Moraxella MPG", 
												 "Streptococcus MPG1", "Mixed MPG", "Staphylococcus MPG", 
												 "Alloiococcus MPG", "Corynebacterium MPG", 
												 "Moraxellaceae MPG1", "Corynebacterium rare MPG",
												 "Haemophilus rare MPG", "Moraxella rare MPG",
												 "Moraxellaceae rare MPG", "Streptococcus rare MPG"))
results.df <- results.df[match(MPG.graph.order, results.df$mpg),]
write.csv(results.df, 
					"output/2018-02-20_mpg.age.results.csv", quote=F,
					row.names=F)


# Look for trend between MPGs and time, with well.sick (samples restricted to age < 2)
# ---------------------------------------
MPG.list <- MPG.alpha.order

mpg.prop.well.age.df <- 
	as.data.frame.matrix(with(log.table.df[which(log.table.df$well.sick == "Well"),], 
														prop.table(table(MPG.name, age.custom),2)))
mpg.prop.sick.age.df <- 
	as.data.frame.matrix(with(log.table.df[which(log.table.df$well.sick == "Sick"),], 
														prop.table(table(MPG.name, age.custom),2)))

results.df <- data.frame(mpg = MPG.list, 
												 rho.cor.well = character(length(MPG.list)),
												 p.cor.well = character(length(MPG.list)),
												 OR.gee.well = character(length(MPG.list)),
												 p.gee.well = character(length(MPG.list)),
												 rho.cor.sick = character(length(MPG.list)),
												 p.cor.sick = character(length(MPG.list)),
												 OR.gee.sick = character(length(MPG.list)),
												 p.gee.sick = character(length(MPG.list)),												 
												 stringsAsFactors=F)

for (i in 1:length(MPG.list)){
	# Trend
	df <- data.frame(mpg = t(mpg.prop.well.age.df[
		which(row.names(mpg.prop.well.age.df) == MPG.list[i]),]),
		time = c(2/12, 6/12, 1, 1.5, 2, 2.5, 3, 3.5, 4))
	df <- df %>% filter(time <= 2)
	colnames(df) <- c("mpg", "time")
	cor.results <- with(df, cor.test(mpg, time, method = "spearman"))
	rho.cor <- cor.results$estimate
	p.cor <- cor.results$p.value
	results.df$rho.cor.well[i] <- good.signif(rho.cor, 2)
	results.df$p.cor.well[i] <- good.signif(p.cor, 2) 
	# gee
	df <- log.table.df[which(log.table.df$well.sick == "Well"),]
	df$age.custom <- 	as.numeric(
		as.character(factor(df$age.custom,labels=c(2/12, 6/12, 1, 1.5, 2, 2.5, 3, 3.5, 4))))
	df <- df %>% filter(age.custom <= 2)
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														age + sex + season,
													id = SUBJ_ID_ID, 
													data = df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR.gee.well[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p.gee.well[i] <- paste0(p)[2]
}

for (i in 1:length(MPG.list)){
	# Trend
	df <- data.frame(mpg = t(mpg.prop.sick.age.df[
		which(row.names(mpg.prop.sick.age.df) == MPG.list[i]),]),
		time = c(2/12, 6/12, 1, 1.5, 2, 2.5, 3, 3.5, 4))
	df <- df %>% filter(time <= 2)
	colnames(df) <- c("mpg", "time")
	cor.results <- with(df, cor.test(mpg, time, method = "spearman"))
	rho.cor <- cor.results$estimate
	p.cor <- cor.results$p.value
	results.df$rho.cor.sick[i] <- good.signif(rho.cor, 2)
	results.df$p.cor.sick[i] <- good.signif(p.cor, 2) 
	# gee
	df <- log.table.df[which(log.table.df$well.sick == "Sick"),]
	df$age.custom <- 	as.numeric(
		as.character(factor(df$age.custom,labels=c(2/12, 6/12, 1, 1.5, 2, 2.5, 3, 3.5, 4))))
	df <- df %>% filter(age.custom <= 2)
	gee.out <- tryCatch(gee(I(MPG.name == MPG.list[i]) ~  
														age + sex + season,
													id = SUBJ_ID_ID, 
													data = df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR.gee.sick[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p.gee.sick[i] <- paste0(p)[2]
}

MPG.graph.order <- rev(c("Haemophilus MPG1", "Haemophilus MPG2", "Moraxella MPG", 
												 "Streptococcus MPG1", "Mixed MPG", "Staphylococcus MPG", 
												 "Alloiococcus MPG", "Corynebacterium MPG", 
												 "Moraxellaceae MPG1", "Corynebacterium rare MPG",
												 "Haemophilus rare MPG", "Moraxella rare MPG",
												 "Moraxellaceae rare MPG", "Streptococcus rare MPG"))
results.df <- results.df[match(MPG.graph.order, results.df$mpg),]
write.csv(results.df, 
					"output/2018-02-20_mpg.age.results.age2.csv", quote=F,
					row.names=F)


# Run GEE for LRI vs. URI ~ MPG, adjusting for sex, age
# with id = SUBJ_ID
# ---------------------------------------

log.table.unwell.df <- log.table.df %>%
	filter(log.table.df$Category == "infectious") 
levels(log.table.unwell.df$illness.type) <- c("URI", "URI", "LRI")

log.table.unwell.df <- log.table.unwell.df[order(log.table.unwell.df$LogNo),]

mpg.count.per.illness.type.df <- 
	as.data.frame.matrix(
		with(log.table.unwell.df, table(MPG.name, illness.type)))
colnames(mpg.count.per.illness.type.df) <- c("URI.Freq","LRI.Freq")
mpg.prop.per.illness.type.df <- 
	as.data.frame.matrix(
		with(log.table.unwell.df, prop.table(table(MPG.name, illness.type),2)))
colnames(mpg.prop.per.illness.type.df) <- c("URI.Prop","LRI.Prop")

mpg.distr.illness.type.df <- 
	data.frame(MPG.name = row.names(mpg.prop.per.illness.type.df),
						 URI.Freq = mpg.count.per.illness.type.df$URI.Freq,
						 URI.Prop = paste0(good.signif(mpg.prop.per.illness.type.df$URI.Prop*100,2), "%"),
						 LRI.Freq = mpg.count.per.illness.type.df$LRI.Freq,
						 LRI.Prop = paste0(good.signif(mpg.prop.per.illness.type.df$LRI.Prop*100,2), "%"), 
						 stringsAsFactors = F)
mpg.distr.illness.type.df2 <- 
	data.frame(MPG.name = mpg.distr.illness.type.df$MPG.name,
						 URI = paste0(mpg.distr.illness.type.df$URI.Freq, " (", 
						 						 mpg.distr.illness.type.df$URI.Prop, ")"),
						 LRI = paste0(mpg.distr.illness.type.df$LRI.Freq, " (", 
						 						 mpg.distr.illness.type.df$LRI.Prop, ")"))

MPG.list <- row.names(mpg.prop.per.illness.type.df)

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(illness.type=="LRI") ~  
								 	I(MPG.name == MPG.list[i]) + sex + age,
								 id = SUBJ_ID_ID, 
								 data = log.table.unwell.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}


mpg.distr.illness.type.df2$OR <- results.df$OR
mpg.distr.illness.type.df2$p <- results.df$p
mpg.distr.illness.type.df2 <- mpg.distr.illness.type.df2[
	match(MPG.alpha.order, mpg.distr.illness.type.df2$MPG.name),]

write.csv(mpg.distr.illness.type.df2, 
					"output/2018-02-20_mpg.distr.illness.type.results.csv", quote=F,
					row.names=F)

# Check Strep ESV vs. other ESV association with LRI vs. URI
# ---------------------------------------
# Streptococcus ESV2 is Streptococcus..357551f644064bb2b49614b1caab2fa5 
# Other ESVs is Others..rare
# dataframe with LogNo, MPG.name, ESVs
dim(log.esv.df)
log.esv.df.StrepMPG2 <- with(log.esv.df, 
														 data.frame(LogNo, MPG.name,
														 					 Streptococcus..357551f644064bb2b49614b1caab2fa5,
														 					 Others..rare))
colnames(log.esv.df.StrepMPG2) <- c("LogNo", "MPG.name", "StrepESV2", "OtherESVs")
# log.esv.df.StrepMPG2.m <- melt(log.esv.df.StrepMPG2, 
# 															 id.vars = c("LogNo", "MPG.name"),
# 															 variable.name = "ESV",
# 															 value.name = "Rel.Abundance")
log.table.esv.df <- merge(log.table.df, log.esv.df.StrepMPG2, 
													by=c("LogNo", "MPG.name"), all = T)
dim(log.table.esv.df)

log.table.esv.df %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=well.sick, y=StrepESV2)) + geom_boxplot()
log.table.esv.df %>%
	filter(!is.na(well.sick)) %>%
	ggplot(aes(x=well.sick, y=OtherESVs)) + geom_boxplot()

log.table.esv.StrepMPG2.df <- log.table.esv.df %>%
	filter(MPG.name == "Mixed MPG", !is.na(well.sick))

log.table.esv.StrepMPG2.df %>%
	ggplot(aes(x=well.sick, y=StrepESV2)) + geom_boxplot()
log.table.esv.StrepMPG2.df%>%
	ggplot(aes(x=well.sick, y=OtherESVs)) + geom_boxplot()

log.table.esv.StrepMPG2.df %>% filter(well.sick == "Sick") %>%
	ggplot(aes(x=illness.type, y=StrepESV2)) + geom_boxplot()
log.table.esv.StrepMPG2.df %>% filter(well.sick == "Sick") %>%
	ggplot(aes(x=illness.type, y=OtherESVs)) + geom_boxplot()

gee.out <- gee(I(well.sick == "Sick") ~  
							 	StrepESV2 + OtherESVs + sex + age,
							 id = SUBJ_ID_ID, 
							 data = log.table.esv.StrepMPG2.df, 
							 family= binomial,
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
data.frame(variable = row.names(s$coefficients),
					 OR = paste0(OR, " (", lci, "-", uci, ")"),
					 p = paste0(p))

gee.out <- gee(I(well.sick == "Sick") ~  
							 	StrepESV2 * OtherESVs + sex + age,
							 id = SUBJ_ID_ID, 
							 data = log.table.esv.StrepMPG2.df, 
							 family= binomial,
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
data.frame(variable = row.names(s$coefficients),
					 OR = paste0(OR, " (", lci, "-", uci, ")"),
					 p = paste0(p))

# In Mixed MPG, Streptococcus ESV2 is predictive for disease! (p=0.019, OR 7).
# No evidence for interaction

log.table.esv.df %>%
	filter(!is.na(well.sick), well.sick != "Well") %>%
	ggplot(aes(x=illness.type, y=StrepESV2)) + geom_boxplot()
log.table.esv.df %>%
	filter(!is.na(well.sick), well.sick != "Well") %>%
	ggplot(aes(x=illness.type, y=OtherESVs)) + geom_boxplot()

log.table.esv.StrepMPG2.sick.df <- log.table.esv.df %>%
	filter(MPG.name == "Mixed MPG", !is.na(well.sick), 
				 well.sick != "Well")

log.table.esv.StrepMPG2.sick.df %>%
	ggplot(aes(x=illness.type, y=StrepESV2)) + geom_boxplot()
log.table.esv.StrepMPG2.sick.df%>%
	ggplot(aes(x=illness.type, y=OtherESVs)) + geom_boxplot()

gee.out <- gee(I(illness.type == "LRI") ~  
							 	StrepESV2 + OtherESVs + sex + age,
							 id = SUBJ_ID_ID, 
							 data = log.table.esv.StrepMPG2.sick.df, 
							 family= binomial,
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
data.frame(variable = row.names(s$coefficients),
					 OR = paste0(OR, " (", lci, "-", uci, ")"),
					 p = paste0(p))

gee.out <- gee(I(illness.type == "LRI") ~  
							 	StrepESV2 * OtherESVs + sex + age,
							 id = SUBJ_ID_ID, 
							 data = log.table.esv.StrepMPG2.df, 
							 family= binomial,
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
data.frame(variable = row.names(s$coefficients),
					 OR = paste0(OR, " (", lci, "-", uci, ")"),
					 p = paste0(p))

# ------------------------------------------------------------------------------
# Seasonal trends
# ------------------------------------------------------------------------------

# Get proportions by month
# ---------------------------------------

month.mpg.well.df <- 
	as.data.frame.matrix(
		with(log.table.df[which(log.table.df$well.sick == "Well"),], 
				 prop.table(table(month, MPG.name.top),1))
	)
month.mpg.well.df$month.name <- factor(c(1:12),
																			 labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
																			 				 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
month.mpg.well.df.m <-
	melt(month.mpg.well.df, 
			 id.vars = "month.name",
			 variable.name = "MPG.name",
			 value.name = "Prop")
month.mpg.well.df.m$well.sick <- "Well"

month.mpg.sick.df <- 
	as.data.frame.matrix(
		with(log.table.df[which(log.table.df$well.sick == "Sick"),], 
				 prop.table(table(month, MPG.name.top),1))
	)
month.mpg.sick.df$month.name <- 
	factor(c(1:12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
													 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
month.mpg.sick.df.m <-
	melt(month.mpg.sick.df, 
			 id.vars = "month.name",
			 variable.name = "MPG.name",
			 value.name = "Prop")
month.mpg.sick.df.m$well.sick <- "Sick"

month.mpg.well.sick.df <- rbind(month.mpg.well.df.m, month.mpg.sick.df.m)
month.mpg.well.sick.df$well.sick <- 
	factor(month.mpg.well.sick.df$well.sick, 
				 levels=c("Well", "Sick"))

pdf("output/2018-02-20_MPG_distr_month.pdf", width=10, height=5)
month.mpg.well.sick.df %>%
	ggplot(aes(x=month.name, y=Prop, group=MPG.name)) + 
	geom_line(aes(color=MPG.name)) + 
	geom_point(aes(color=MPG.name)) + 
	scale_color_manual(values = levels(log.table.df$MPG.colour.hex.top)) +
	facet_grid(.~well.sick) +
	labs(color = "MPG name",
					 x = "Month of collection",
			 y = "Proportion") +
	theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank()) + 
	ylim(0, 0.7) +
	geom_rect(xmin=0, xmax=2.5, ymin=0.65, ymax=0.7,
						fill="#ff8384") + 	
	geom_rect(xmin=2.5, xmax=5.5, ymin=0.65, ymax=0.7,
						fill="#edb579") + 	
	geom_rect(xmin=5.5, xmax=8.5, ymin=0.65, ymax=0.7,
						fill="#ace4ea") + 
	geom_rect(xmin=8.5, xmax=11.5, ymin=0.65, ymax=0.7,
						fill="#00cd66") + 
	geom_rect(xmin=11.5, xmax=12.5, ymin=0.65, ymax=0.7,
						fill="#ff8384") + 	
	geom_vline(xintercept = c(2.5, 5.5, 8.5, 11.5)) + 
	annotate("text", x=1.5, y=0.675, label="Summer", hjust=0.5) + 
	annotate("text", x=4, y=0.675, label="Fall", hjust=0.5) +
	annotate("text", x=7, y=0.675, label="Winter", hjust=0.5) +
	annotate("text", x=10, y=0.675, label="Spring", hjust=0.5)
dev.off()

a <- month.mpg.well.df.m[,1:3]
colnames(a) <- c("month.name", "MPG.name", "Well.Prop")
b <- month.mpg.sick.df.m[,1:3]
colnames(b) <- c("month.name", "MPG.name", "Sick.Prop")

month.mpg.well.sick.df.or <- 
	merge(a, b, by=c("month.name", "MPG.name"), all=T)
month.mpg.well.sick.df.or$Sick.well.ratio <- 
	month.mpg.well.sick.df.or$Sick.Prop / month.mpg.well.sick.df.or$Well.Prop
month.mpg.well.sick.df.or$Sick.well.ratio[
	is.infinite(month.mpg.well.sick.df.or$Sick.well.ratio)] <- NA
month.mpg.well.sick.df.or$Sick.well.ratio[
	month.mpg.well.sick.df.or$Sick.well.ratio == 0] <- NA
month.mpg.well.sick.df.or$Sick.well.ratio <-
	log10(month.mpg.well.sick.df.or$Sick.well.ratio)

pdf("output/2018-02-20_MPG_or_month.pdf", width=10, height=5)
month.mpg.well.sick.df.or %>%
	ggplot(aes(x=month.name, group=MPG.name)) + 
	geom_line(aes(y=Sick.well.ratio, color=MPG.name)) +
	geom_line(data=month.mpg.well.sick.df.or[
		which(!is.na(month.mpg.well.sick.df.or$Sick.well.ratio)),],
		aes(y=Sick.well.ratio, color=MPG.name), linetype="dotted") +
	geom_point(aes(y=Sick.well.ratio, color=MPG.name)) + 
	scale_color_manual(values = levels(log.table.df$MPG.colour.hex.top)) +
	labs(color = "MPG name",
			 x = "Month of collection",
			 y = "log10(%Sick/%Well)") +
	theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank()) + 
	ylim(-1.5,1.8) + 
	geom_rect(xmin=0, xmax=2.5, ymin=1.5, ymax=1.8,
						fill="#ff8384") + 	
	geom_rect(xmin=2.5, xmax=5.5, ymin=1.5, ymax=1.8,
						fill="#edb579") + 	
	geom_rect(xmin=5.5, xmax=8.5, ymin=1.5, ymax=1.8,
						fill="#ace4ea") + 
	geom_rect(xmin=8.5, xmax=11.5, ymin=1.5, ymax=1.8,
						fill="#00cd66") + 
	geom_rect(xmin=11.5, xmax=12.5, ymin=1.5, ymax=1.8,
						fill="#ff8384") + 	
	geom_vline(xintercept = c(2.5, 5.5, 8.5, 11.5)) + 
	geom_hline(yintercept = 0, linetype="dashed") +
	annotate("text", x=1.5, y=1.65, label="Summer", hjust=0.5) + 
	annotate("text", x=4, y=1.65, label="Fall", hjust=0.5) +
	annotate("text", x=7, y=1.65, label="Winter", hjust=0.5) +
	annotate("text", x=10, y=1.65, label="Spring", hjust=0.5)
dev.off()

# Get proportions by season
# ---------------------------------------
season.mpg.well.df <- 
	as.data.frame.matrix(
		with(log.table.df[which(log.table.df$well.sick == "Well"),], 
				 prop.table(table(season, MPG.name.top),1))
	)
season.mpg.well.df$season.name <- row.names(season.mpg.well.df)
season.mpg.well.df.m <-
	melt(season.mpg.well.df, 
			 id.vars = "season.name",
			 variable.name = "MPG.name",
			 value.name = "Prop")
season.mpg.well.df.m$well.sick <- "Well"

season.mpg.sick.df <- 
	as.data.frame.matrix(
		with(log.table.df[which(log.table.df$well.sick == "Sick"),], 
				 prop.table(table(season, MPG.name.top),1))
	)
season.mpg.sick.df$season.name <- row.names(season.mpg.sick.df)

season.mpg.sick.df.m <-
	melt(season.mpg.sick.df, 
			 id.vars = "season.name",
			 variable.name = "MPG.name",
			 value.name = "Prop")
season.mpg.sick.df.m$well.sick <- "Sick"

season.mpg.well.sick.df <- rbind(season.mpg.well.df.m, season.mpg.sick.df.m)
season.mpg.well.sick.df$well.sick <- 
	factor(season.mpg.well.sick.df$well.sick, 
				 levels=c("Well", "Sick"))
season.mpg.well.sick.df$season.name <- 
	factor(season.mpg.sick.df$season.name, 
				 levels=c("Winter", "Spring", "Summer", "Autumn"),
				 labels=c("Winter", "Spring", "Summer", "Fall"))

pdf("output/2018-02-20_MPG_distr_season.pdf", width=10, height=5)
season.mpg.well.sick.df %>%
	ggplot(aes(x=season.name, y=Prop, group=MPG.name)) + 
	geom_line(aes(color=MPG.name)) + 
	geom_point(aes(color=MPG.name)) + 
	scale_color_manual(values = levels(log.table.df$MPG.colour.hex.top)) +
	facet_grid(.~well.sick) +
	labs(color = "MPG name",
			 x = "Season of collection",
			 y = "Proportion")
dev.off()

a <- season.mpg.well.df.m[,1:3]
colnames(a) <- c("season.name", "MPG.name", "Well.Prop")
b <- season.mpg.sick.df.m[,1:3]
colnames(b) <- c("season.name", "MPG.name", "Sick.Prop")

season.mpg.well.sick.df.or <- 
	merge(a, b, by=c("season.name", "MPG.name"), all=T)

season.mpg.well.sick.df.or$season.name <- 
	factor(season.mpg.well.sick.df.or$season.name, 
				 levels=c("Winter", "Spring", "Summer", "Autumn"),
				 labels=c("Winter", "Spring", "Summer", "Fall"))

season.mpg.well.sick.df.or$Sick.well.ratio <- 
	season.mpg.well.sick.df.or$Sick.Prop / season.mpg.well.sick.df.or$Well.Prop
season.mpg.well.sick.df.or$Sick.well.ratio[
	is.infinite(season.mpg.well.sick.df.or$Sick.well.ratio)] <- NA
season.mpg.well.sick.df.or$Sick.well.ratio[
	season.mpg.well.sick.df.or$Sick.well.ratio == 0] <- NA
season.mpg.well.sick.df.or$Sick.well.ratio <-
	log10(season.mpg.well.sick.df.or$Sick.well.ratio)

pdf("output/2018-02-20_MPG_or_season.pdf", width=10, height=5)
season.mpg.well.sick.df.or %>%
	ggplot(aes(x=season.name, group=MPG.name)) + 
	geom_line(aes(y=Sick.well.ratio, color=MPG.name)) +
	geom_line(data=season.mpg.well.sick.df.or[
		which(!is.na(season.mpg.well.sick.df.or$Sick.well.ratio)),],
		aes(y=Sick.well.ratio, color=MPG.name), linetype="dotted") +
	geom_point(aes(y=Sick.well.ratio, color=MPG.name)) + 
	scale_color_manual(values = levels(log.table.df$MPG.colour.hex.top)) +
	labs(color = "MPG name",
			 x = "Season of collection",
			 y = "log10(%Sick/%Well)") +
	geom_hline(yintercept = 0, linetype="dashed")
dev.off()

# Run GEE for Season (Summer/Fall vs. Winter/Spring) ~ MPG, 
# adjusting for sex, age, with id = SUBJ_ID
# +/- SEGREGATED BY WELL/SICK
# ---------------------------------------

log.table.df$season2 <- log.table.df$season
levels(log.table.df$season2) <- c("Winter.Spring", "Winter.Spring",
																	"Summer.Fall", "Summer.Fall")

mpg.count.per.season.df <- 
	as.data.frame.matrix(
		with(log.table.df, table(MPG.name, season)))
colnames(mpg.count.per.season.df) <- c("Winter.Freq","Spring.Freq","Summer.Freq","Fall.Freq")
mpg.prop.per.season.df <- 
	as.data.frame.matrix(
		with(log.table.df, prop.table(table(MPG.name, season),2)))
colnames(mpg.prop.per.season.df) <- c("Winter.Prop","Spring.Prop","Summer.Prop","Fall.Prop")

mpg.distr.season.df <- 
	data.frame(MPG.name = row.names(mpg.prop.per.season.df),
						 Winter.Freq = mpg.count.per.season.df$Winter.Freq,
						 Winter.Prop = paste0(good.signif(mpg.prop.per.season.df$Winter.Prop*100,2), "%"),
						 Spring.Freq = mpg.count.per.season.df$Spring.Freq,
						 Spring.Prop = paste0(good.signif(mpg.prop.per.season.df$Spring.Prop*100,2), "%"), 
						 Summer.Freq = mpg.count.per.season.df$Summer.Freq,
						 Summer.Prop = paste0(good.signif(mpg.prop.per.season.df$Summer.Prop*100,2), "%"), 	
						 Fall.Freq = mpg.count.per.season.df$Fall.Freq,
						 Fall.Prop = paste0(good.signif(mpg.prop.per.season.df$Fall.Prop*100,2), "%"), 	
						 stringsAsFactors = F)
mpg.distr.season.df2 <- 
	data.frame(MPG.name = mpg.distr.season.df$MPG.name,
						 Winter = paste0(mpg.distr.season.df$Winter.Freq, " (", 
						 								mpg.distr.season.df$Winter.Prop, ")"),
						 Spring = paste0(mpg.distr.season.df$Spring.Freq, " (", 
						 								mpg.distr.season.df$Spring.Prop, ")"),
						 Summer = paste0(mpg.distr.season.df$Summer.Freq, " (", 
						 								mpg.distr.season.df$Summer.Prop, ")"),
						 Fall = paste0(mpg.distr.season.df$Fall.Freq, " (", 
						 							mpg.distr.season.df$Fall.Prop, ")"))

MPG.list <- row.names(mpg.prop.per.season.df)

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(season2=="Winter.Spring") ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + well.sick,
								 id = SUBJ_ID_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.distr.season.df2$OR <- results.df$OR
mpg.distr.season.df2$p <- results.df$p

mpg.distr.season.df2 <- mpg.distr.season.df2[
	match(MPG.alpha.order, mpg.distr.season.df2$MPG.name),]

write.csv(mpg.distr.season.df2, 
					"output/2018-02-20_mpg.distr.season.results.csv", quote=F,
					row.names=F)

# Well/sick x season interaction
# --------------------------------------
# Calculate MPG distribution for visittype
# ---------------------------------------
# samples per season
mpg.count.winter.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Winter"),], 
			 table(MPG.name, well.sick)))
colnames(mpg.count.winter.well.sick.df) <- c("Winter.Well.Freq","Winter.Sick.Freq")

mpg.count.spring.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Spring"),], 
			 table(MPG.name, well.sick)))
colnames(mpg.count.spring.well.sick.df) <- c("Spring.Well.Freq","Spring.Sick.Freq")

mpg.count.summer.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Summer"),], 
			 table(MPG.name, well.sick)))
colnames(mpg.count.summer.well.sick.df) <- c("Summer.Well.Freq","Summer.Sick.Freq")

mpg.count.fall.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Autumn"),], 
			 table(MPG.name, well.sick)))
colnames(mpg.count.fall.well.sick.df) <- c("Fall.Well.Freq","Fall.Sick.Freq")

mpg.prop.winter.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Winter"),], 
			 prop.table(table(MPG.name, well.sick),2)))
colnames(mpg.prop.winter.well.sick.df) <- c("Winter.Well.Prop","Winter.Sick.Prop")

mpg.prop.spring.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Spring"),], 
			 prop.table(table(MPG.name, well.sick),2)))
colnames(mpg.prop.spring.well.sick.df) <- c("Spring.Well.Prop","Spring.Sick.Prop")

mpg.prop.summer.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Summer"),], 
			 prop.table(table(MPG.name, well.sick),2)))
colnames(mpg.prop.summer.well.sick.df) <- c("Summer.Well.Prop","Summer.Sick.Prop")

mpg.prop.fall.well.sick.df <- as.data.frame.matrix(
	with(log.table.df[which(log.table.df$season == "Autumn"),], 
			 prop.table(table(MPG.name, well.sick),2)))
colnames(mpg.prop.fall.well.sick.df) <- c("Fall.Well.Prop","Fall.Sick.Prop")

mpg.season.well.sick.distr.df <- data.frame(
	MPG.name = row.names(mpg.count.winter.well.sick.df),
	Winter.Well.Freq = mpg.count.winter.well.sick.df$Winter.Well.Freq,
	Winter.Well.Prop = paste0(good.signif(mpg.prop.winter.well.sick.df$Winter.Well.Prop*100,2), "%"),
	Winter.Sick.Freq = mpg.count.winter.well.sick.df$Winter.Sick.Freq,
	Winter.Sick.Prop = paste0(good.signif(mpg.prop.winter.well.sick.df$Winter.Sick.Prop*100,2), "%"),
	Spring.Well.Freq = mpg.count.spring.well.sick.df$Spring.Well.Freq,
	Spring.Well.Prop = paste0(good.signif(mpg.prop.spring.well.sick.df$Spring.Well.Prop*100,2), "%"),
	Spring.Sick.Freq = mpg.count.spring.well.sick.df$Spring.Sick.Freq,
	Spring.Sick.Prop = paste0(good.signif(mpg.prop.spring.well.sick.df$Spring.Sick.Prop*100,2), "%"),
	Summer.Well.Freq = mpg.count.summer.well.sick.df$Summer.Well.Freq,
	Summer.Well.Prop = paste0(good.signif(mpg.prop.summer.well.sick.df$Summer.Well.Prop*100,2), "%"),
	Summer.Sick.Freq = mpg.count.summer.well.sick.df$Summer.Sick.Freq,
	Summer.Sick.Prop = paste0(good.signif(mpg.prop.summer.well.sick.df$Summer.Sick.Prop*100,2), "%"),
	Fall.Well.Freq = mpg.count.fall.well.sick.df$Fall.Well.Freq,
	Fall.Well.Prop = paste0(good.signif(mpg.prop.fall.well.sick.df$Fall.Well.Prop*100,2), "%"),
	Fall.Sick.Freq = mpg.count.fall.well.sick.df$Fall.Sick.Freq,
	Fall.Sick.Prop = paste0(good.signif(mpg.prop.fall.well.sick.df$Fall.Sick.Prop*100,2), "%"),
	stringsAsFactors = F)

mpg.season.well.sick.distr.df2 <- 
	data.frame(MPG.name = mpg.season.well.sick.distr.df$MPG.name,
						 Winter.Well = paste0(mpg.season.well.sick.distr.df$Winter.Well.Freq, " (",
						 										 mpg.season.well.sick.distr.df$Winter.Well.Prop, ")"),
						 Winter.Sick = paste0(mpg.season.well.sick.distr.df$Winter.Sick.Freq, " (",
						 										 mpg.season.well.sick.distr.df$Winter.Sick.Prop, ")"),
						 Spring.Well = paste0(mpg.season.well.sick.distr.df$Spring.Well.Freq, " (",
						 										 mpg.season.well.sick.distr.df$Spring.Well.Prop, ")"),
						 Spring.Sick = paste0(mpg.season.well.sick.distr.df$Spring.Sick.Freq, " (",
						 										 mpg.season.well.sick.distr.df$Spring.Sick.Prop, ")"),
						 Summer.Well = paste0(mpg.season.well.sick.distr.df$Summer.Well.Freq, " (",
						 										 mpg.season.well.sick.distr.df$Summer.Well.Prop, ")"),
						 Summer.Sick = paste0(mpg.season.well.sick.distr.df$Summer.Sick.Freq, " (",
						 										 mpg.season.well.sick.distr.df$Summer.Sick.Prop, ")"),
						 Fall.Well = paste0(mpg.season.well.sick.distr.df$Fall.Well.Freq, " (",
						 									 mpg.season.well.sick.distr.df$Fall.Well.Prop, ")"),
						 Fall.Sick = paste0(mpg.season.well.sick.distr.df$Fall.Sick.Freq, " (",
						 									 mpg.season.well.sick.distr.df$Fall.Sick.Prop, ")"))

# Run GEE for Sick vs. Well ~ MPG, adjusting for sex, age and season
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- unique(mpg.season.well.sick.distr.df2$MPG.name)

results.df <- data.frame(OR_MPG = character(length(MPG.list)),
												 p_MPG = character(length(MPG.list)),
												 OR_season = character(length(MPG.list)),
												 p_season = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(well.sick=="Sick") ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + I(season2 == "Winter.Spring"),
								 id = SUBJ_ID_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR_MPG[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p_MPG[i] <- paste0(p)[2]
	results.df$OR_season[i] <- paste0(OR, " (", lci, "-", uci, ")")[5]
	results.df$p_season[i] <- paste0(p)[5]	
}

mpg.season.well.sick.distr.df2$OR_MPG <- results.df$OR_MPG
mpg.season.well.sick.distr.df2$p_MPG <- results.df$p_MPG
mpg.season.well.sick.distr.df2$OR_season <- results.df$OR_season
mpg.season.well.sick.distr.df2$p_season <- results.df$p_season

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(well.sick=="Sick") ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + season,
								 id = SUBJ_ID_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.season.well.sick.distr.df2$OR_MPG_alt <- results.df$OR
mpg.season.well.sick.distr.df2$p_MPG_alt <- results.df$p

write.csv(mpg.season.well.sick.distr.df2, "output/2018-02-20_mpg.distr.season.well.sick.csv", quote=F,
					row.names=F)

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- tryCatch(gee(I(well.sick=="Sick") ~  
														I(MPG.name == MPG.list[i])*I(season2 == "Winter.Spring") + sex + age,
													id = SUBJ_ID_ID, 
													data = log.table.df, 
													family= binomial,
													corstr = "unstructured"),
											error=function(err) NA)
	if(is.na(gee.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	}
	else{
		s  <- summary(gee.out)
	}
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,4]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,4]))
	p <- 2 * pnorm(abs(s$coefficients[,5]), lower.tail = FALSE)
	OR <- good.signif(OR,2)
	lci <- good.signif(lci,2)
	uci <- good.signif(uci,2)
	p <- signif(p,2)
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[6]
	results.df$p[i] <- paste0(p)[6]
}

mpg.season.well.sick.distr.df2$OR_MPG_int <- results.df$OR
mpg.season.well.sick.distr.df2$p_MPG_int <- results.df$p

results.df <- data.frame(OR_sex = character(length(MPG.list)),
												 p_sex = character(length(MPG.list)),
												 OR_age = character(length(MPG.list)),
												 p_age = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(well.sick=="Sick") ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + season,
								 id = SUBJ_ID_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR_sex[i] <- paste0(OR, " (", lci, "-", uci, ")")[3]
	results.df$p_sex[i] <- paste0(p)[3]
	results.df$OR_age[i] <- paste0(OR, " (", lci, "-", uci, ")")[4]
	results.df$p_age[i] <- paste0(p)[4]
}

mpg.season.well.sick.distr.df2$OR_MPG_sex <- results.df$OR_sex
mpg.season.well.sick.distr.df2$p_MPG_sex <- results.df$p_sex
mpg.season.well.sick.distr.df2$OR_MPG_age <- results.df$OR_age
mpg.season.well.sick.distr.df2$p_MPG_age <- results.df$p_age
mpg.season.well.sick.distr.df2 <- mpg.season.well.sick.distr.df2[
	match(MPG.sick.sig.order, mpg.season.well.sick.distr.df2$MPG.name),]
write.csv(mpg.season.well.sick.distr.df2,
					"output/2018-02-20_mpg.season.well.sick.distr.csv", quote=F,
					row.names=F)
# ------------------------------------------------------------------------------
# Associations with pet ownership
# ------------------------------------------------------------------------------

mpg.count.dog.df <- as.data.frame.matrix(
	with(log.table.df, table(MPG.name, dog0)))
colnames(mpg.count.dog.df) <- c("No.Dog.Freq","Dog.Freq")
mpg.prop.dog.df <- as.data.frame.matrix(
	with(log.table.df, prop.table(table(MPG.name, dog0),2)))
colnames(mpg.prop.dog.df) <- c("No.Dog.Prop","Dog.Prop")

mpg.count.cat.df <- as.data.frame.matrix(
	with(log.table.df, table(MPG.name, cat0)))
colnames(mpg.count.cat.df) <- c("No.Cat.Freq","Cat.Freq")
mpg.prop.cat.df <- as.data.frame.matrix(
	with(log.table.df, prop.table(table(MPG.name, cat0),2)))
colnames(mpg.prop.cat.df) <- c("No.Cat.Prop","Cat.Prop")

mpg.pet.df <- 
	data.frame(MPG.name = row.names(mpg.count.dog.df),
						 No.Dog.Freq = mpg.count.dog.df$No.Dog.Freq,
						 No.Dog.Prop = paste0(good.signif(mpg.prop.dog.df$No.Dog.Prop*100,2), "%"),
						 Dog.Freq = mpg.count.dog.df$Dog.Freq,
						 Dog.Prop = paste0(good.signif(mpg.prop.dog.df$Dog.Prop*100,2), "%"),
						 No.Cat.Freq = mpg.count.cat.df$No.Cat.Freq,
						 No.Cat.Prop = paste0(good.signif(mpg.prop.cat.df$No.Cat.Prop*100,2), "%"),
						 Cat.Freq = mpg.count.cat.df$Cat.Freq,
						 Cat.Prop = paste0(good.signif(mpg.prop.cat.df$Cat.Prop*100,2), "%"),
						 stringsAsFactors = F)
mpg.pet.df2 <- 
	data.frame(MPG.name = mpg.pet.df$MPG.name,
						 No.Dog = paste0(mpg.pet.df$No.Dog.Freq, " (", mpg.pet.df$No.Dog.Prop, ")"),
						 Dog = paste0(mpg.pet.df$Dog.Freq, " (",mpg.pet.df$Dog.Prop, ")"),
						 No.Cat = paste0(mpg.pet.df$No.Cat.Freq, " (",mpg.pet.df$No.Cat.Prop, ")"),
						 Cat = paste0(mpg.pet.df$Cat.Freq, " (",mpg.pet.df$Cat.Prop, ")"))

# Run GEE for Sick vs. Well ~ MPG, adjusting for sex, age
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- mpg.pet.df$MPG.name

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(dog0) ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + well.sick,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.pet.df2$OR_dog <- results.df$OR
mpg.pet.df2$p_dog <- results.df$p

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(cat0) ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + well.sick,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.pet.df2$OR_cat <- results.df$OR
mpg.pet.df2$p_cat <- results.df$p

write.csv(mpg.pet.df2, "output/2018-02-20_mpg.distr.pets.csv", quote=F,
					row.names=F)

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(dog0) ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + well.sick + season2,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.pet.df2$OR_dog <- results.df$OR
mpg.pet.df2$p_dog <- results.df$p

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(cat0) ~  
								 	I(MPG.name == MPG.list[i]) + sex + age + well.sick  + season2,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.pet.df2$OR_cat <- results.df$OR
mpg.pet.df2$p_cat <- results.df$p


# ------------------------------------------------------------------------------
# Associations with delivery mode
# ------------------------------------------------------------------------------

mpg.count.birthby.df <- as.data.frame.matrix(
	with(log.table.df, table(MPG.name, birthby)))
colnames(mpg.count.birthby.df) <- c("Caesar.Freq","Vaginal.Freq")
mpg.prop.birthby.df <- as.data.frame.matrix(
	with(log.table.df, prop.table(table(MPG.name, birthby),2)))
colnames(mpg.prop.birthby.df) <- c("Caesar.Prop","Vaginal.Prop")

mpg.birthby.df <- 
	data.frame(MPG.name = row.names(mpg.count.birthby.df),
						 Caesar.Freq = mpg.count.birthby.df$Caesar.Freq,
						 Caesar.Prop = paste0(good.signif(mpg.prop.birthby.df$Caesar.Prop*100,2), "%"),
						 Vaginal.Freq = mpg.count.birthby.df$Vaginal.Freq,
						 Vaginal.Prop = paste0(good.signif(mpg.prop.birthby.df$Vaginal.Prop*100,2), "%"),
						 stringsAsFactors = F)
mpg.birthby.df2 <- 
	data.frame(MPG.name = mpg.birthby.df$MPG.name,
						 Vaginal = paste0(mpg.birthby.df$Vaginal.Freq, " (",mpg.birthby.df$Vaginal.Prop, ")"),
						 Caesar = paste0(mpg.birthby.df$Caesar.Freq, " (", mpg.birthby.df$Caesar.Prop, ")"),
						 stringsAsFactors = F)

# Run GEE for Sick vs. Well ~ MPG, adjusting for sex, age
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- mpg.birthby.df$MPG.name

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(MPG.name == MPG.list[i])  ~  
								 	
								 	I(birthby=="Cesarean section")+ sex + age + well.sick,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.birthby.df2$OR <- results.df$OR
mpg.birthby.df2$p <- results.df$p

write.csv(mpg.birthby.df2, "output/2018-02-20_mpg.distr.birthby.csv", quote=F,
					row.names=F)

# ------------------------------------------------------------------------------
# Associations with sex/sex
# ------------------------------------------------------------------------------

mpg.count.sex.df <- as.data.frame.matrix(
	with(log.table.df, table(MPG.name, sex)))
# colnames(mpg.count.sex.df) <- c("Female","Male")
mpg.prop.sex.df <- as.data.frame.matrix(
	with(log.table.df, prop.table(table(MPG.name, sex),2)))
# colnames(mpg.prop.sex.df) <- c("Female","Male")

mpg.sex.df <- 
	data.frame(MPG.name = row.names(mpg.count.sex.df),
						 Female.Freq = mpg.count.sex.df$Female,
						 Female.Prop = paste0(good.signif(mpg.prop.sex.df$Female*100,2), "%"),
						 Male.Freq = mpg.count.sex.df$Male,
						 Male.Prop = paste0(good.signif(mpg.prop.sex.df$Male*100,2), "%"),
						 stringsAsFactors = F)
mpg.sex.df2 <- 
	data.frame(MPG.name = mpg.sex.df$MPG.name,
						 Female = paste0(mpg.sex.df$Female.Freq, " (",mpg.sex.df$Female.Prop, ")"),
						 Male = paste0(mpg.sex.df$Male.Freq, " (", mpg.sex.df$Male.Prop, ")"),
						 stringsAsFactors = F)

# Run GEE for sex ~ MPG, adjusting for age, sick/well
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- mpg.sex.df$MPG.name

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(sex=="Male") ~
								 	I(MPG.name == MPG.list[i]) + age + well.sick,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.sex.df2$OR <- results.df$OR
mpg.sex.df2$p <- results.df$p

write.csv(mpg.sex.df2, "output/2018-02-20_mpg.distr.sex.csv", quote=F,
					row.names=F)

# Plot proportions of MPGs per sex
# ---------------------------------------

pdf("output/2018-02-20_MPG_distr_sex_block.pdf", width=8, height=5)
log.table.df %>% 
	filter(!is.na(age.custom)) %>%
	ggplot(aes(x=age.custom, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_bar(position = "fill") + 
	scale_fill_manual(values = levels(log.table.df$MPG.colour.hex.top)) + 
	facet_grid(well.sick~sex) +
	labs(fill = "MPG name",
			 x = "Timepoint of collection",
			 y = "Proportion") +
	theme(axis.text.x=element_text(angle=90, hjust=1))
dev.off()


# ------------------------------------------------------------------------------
# Associations with delivery mode
# ------------------------------------------------------------------------------

mpg.count.birthby.df <- as.data.frame.matrix(
	with(log.table.df, table(MPG.name, birthby)))
colnames(mpg.count.birthby.df) <- c("Caesar.Freq","Vaginal.Freq")
mpg.prop.birthby.df <- as.data.frame.matrix(
	with(log.table.df, prop.table(table(MPG.name, birthby),2)))
colnames(mpg.prop.birthby.df) <- c("Caesar.Prop","Vaginal.Prop")

mpg.birthby.df <- 
	data.frame(MPG.name = row.names(mpg.count.birthby.df),
						 Caesar.Freq = mpg.count.birthby.df$Caesar.Freq,
						 Caesar.Prop = paste0(good.signif(mpg.prop.birthby.df$Caesar.Prop*100,2), "%"),
						 Vaginal.Freq = mpg.count.birthby.df$Vaginal.Freq,
						 Vaginal.Prop = paste0(good.signif(mpg.prop.birthby.df$Vaginal.Prop*100,2), "%"),
						 stringsAsFactors = F)
mpg.birthby.df2 <- 
	data.frame(MPG.name = mpg.birthby.df$MPG.name,
						 Vaginal = paste0(mpg.birthby.df$Vaginal.Freq, " (",mpg.birthby.df$Vaginal.Prop, ")"),
						 Caesar = paste0(mpg.birthby.df$Caesar.Freq, " (", mpg.birthby.df$Caesar.Prop, ")"),
						 stringsAsFactors = F)

# Run GEE for delivery ~ MPG, adjusting for sex, age, sick/well
# with id = SUBJ_ID
# ---------------------------------------

MPG.list <- mpg.birthby.df$MPG.name

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(I(MPG.name == MPG.list[i])  ~  
								 	
								 	I(birthby=="Cesarean section")+ sex + age + well.sick,
								 id = SUBJ_ID, 
								 data = log.table.df, 
								 family= binomial,
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.birthby.df2$OR <- results.df$OR
mpg.birthby.df2$p <- results.df$p

write.csv(mpg.birthby.df2, "output/2018-02-20_mpg.distr.birthby.csv", quote=F,
					row.names=F)

# ------------------------------------------------------------------------------
# Associations with number of siblings
# ------------------------------------------------------------------------------

mpg.livsib.df <- log.table.df %>%
	group_by(MPG.name) %>%
	summarize(mean = mean(livsib, na.rm=T),
						se = sqrt(var(livsib, na.rm=T)/length(na.omit(livsib))))

mpg.livsib.df2 <- 
	data.frame(MPG.name = mpg.livsib.df$MPG.name,
						 Siblings = paste0(good.signif(mpg.livsib.df$mean,2), " (", 
						 									good.signif(mpg.livsib.df$se,2), ")"),
						 stringsAsFactors = F)

# Run GEE for sibling_count ~ MPG, adjusting for sex, age
# with id = SUBJ_ID
# note: family = gaussian
# ---------------------------------------

MPG.list <- mpg.livsib.df2$MPG.name

results.df <- data.frame(OR = character(length(MPG.list)),
												 p = character(length(MPG.list)),
												 stringsAsFactors=F)
for(i in 1:length(MPG.list)){
	gee.out <- gee(livsib ~
								 	I(MPG.name == MPG.list[i]) + age + well.sick,
								 id = SUBJ_ID, 
								 data = log.table.df, 
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
	results.df$OR[i] <- paste0(OR, " (", lci, "-", uci, ")")[2]
	results.df$p[i] <- paste0(p)[2]
}

mpg.livsib.df2$OR <- results.df$OR
mpg.livsib.df2$p <- results.df$p

write.csv(mpg.livsib.df2, "output/2018-02-20_mpg.distr.sibling.csv", quote=F,
					row.names=F)

save.image("output/2018-04-03_micro_assoc_analysis.Rdata")
