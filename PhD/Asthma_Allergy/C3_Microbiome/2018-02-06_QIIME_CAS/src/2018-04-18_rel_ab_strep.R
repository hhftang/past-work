rm(list=ls())
load("output/2018-04-18_hclust_MPG_trajectories.Rdata")

otu.df <- otu_table(cas.biom.method1.css.common)

df <- data.frame(
	LogNo = colnames(otu.df),
	haemoESV2.rel.ab = c(otu.df[
		rownames(otu.df) == "bc0d904a0469d05c65b270998488ee91",]),
	strepESV1.rel.ab = c(otu.df[
		rownames(otu.df) == "4060107ffdf52e44f9d72343b7332609",]),
	haemoESV1.rel.ab = c(otu.df[
		rownames(otu.df) == "f579ddd8d7c6a3d812ca33d72110f0d6",]),	
	moraxESV.rel.ab = c(otu.df[
		rownames(otu.df) == "d253ca966efcd811d057dc2aa6430774",]),
	strepESV2.rel.ab = c(otu.df[
		rownames(otu.df) == "357551f644064bb2b49614b1caab2fa5",]),
	coryESV.rel.ab = c(otu.df[
		rownames(otu.df) == "cb50e79d177c169ec5b9d5e7d308874e",]),
	alloESV.rel.ab = c(otu.df[
		rownames(otu.df) == "dd2e15140d20c84712456e4c12ece08b",]),
	staphESV.rel.ab = c(otu.df[
		rownames(otu.df) == "29eb8e14ce61e9c4549de0f9cf200c19",])
)

df2 <- merge(log.table.df, df, by="LogNo", all=T)

df2 %>%
	ggplot(aes(x=strepESV1.rel.ab, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_histogram(position = "stack", binwidth=0.1) + 
	scale_fill_manual(values =  levels(log.table.df$MPG.colour.hex.top))
df2 %>%
	ggplot(aes(x=strepESV2.rel.ab, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_histogram(position = "stack", binwidth=0.1) + 
	scale_fill_manual(values =  levels(log.table.df$MPG.colour.hex.top))
df2 %>%
	ggplot(aes(x=staphESV.rel.ab, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_histogram(position = "stack", binwidth=0.1) + 
	scale_fill_manual(values =  levels(log.table.df$MPG.colour.hex.top))

df2 %>%
	filter(age.custom == "2mth") %>%
	ggplot(aes(x=strepESV1.rel.ab, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_histogram(position = "stack", binwidth=0.1) + 
	scale_fill_manual(values =  levels(log.table.df$MPG.colour.hex.top))
df2 %>%
	filter(age.custom == "2mth") %>%
	ggplot(aes(x=strepESV2.rel.ab, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_histogram(position = "stack", binwidth=0.1) + 
	scale_fill_manual(values =  levels(log.table.df$MPG.colour.hex.top))
df2 %>%
	filter(age.custom == "2mth") %>%
	ggplot(aes(x=staphESV.rel.ab, group=MPG.name.top, fill=MPG.name.top)) + 
	geom_histogram(position = "stack", binwidth=0.1) + 
	scale_fill_manual(values =  levels(log.table.df$MPG.colour.hex.top))

# Cut-off at 20% makes sense, but I would rather use streptococcus MPG at 2mth 
# as predictor instead - compare anyway

# Get logs at 2mth timepoint that precedes earliest infection/sick sample
df3 <- df2 %>%
	filter(age.custom == "2mth")
# All subj with sick samples in 2mth timepoint
sick.subj <- unique(df3[which(df3$Category == "infectious"),]$SUBJ_ID)

bad.subj.log <- as.numeric()
for (i in 1:length(sick.subj)) {
	age.well <- df3[which(df3$SUBJ_ID == sick.subj[i] & df3$Category == "control"),]$age
	age.sick <- df3[which(df3$SUBJ_ID == sick.subj[i] & df3$Category == "infectious"),]$age
	if (identical(age.well, integer(0))) {
		bad.subj.log[i] <- NA
	} else if (age.well > min(age.sick,na.rm=T)) {
		bad.subj.log[i] <- df3[which(df3$SUBJ_ID == sick.subj[i] & df3$age == age.well),]$LogNo
	} else {
		bad.subj.log[i] <- NA
	}
}

df3.well <- df3 %>% filter(Category == "control")
length(df3.well$SUBJ_ID)
length(unique(df3.well$SUBJ_ID))

df4.well <- df3.well %>%
	filter(!(LogNo %in% bad.subj.log))
# removes 14 "bad" samples that occurred after a preceding illness

# Get average rel ab of all samples fitting criteria (if needed - doubt it)
length(df4.well$SUBJ_ID)
length(unique(df4.well$SUBJ_ID))
# All unique, so no need

# GLM asthma ~ 2-4-6mth MPGs, levels, 
# for strep, staph, and others
# with season, sex as covariates

ESV.list <- c("haemoESV2.rel.ab", "strepESV1.rel.ab", "haemoESV1.rel.ab", "moraxESV.rel.ab", 
							"strepESV2.rel.ab", "coryESV.rel.ab", "alloESV.rel.ab", "staphESV.rel.ab")

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("cras_5 == 1 ~ ",
													ESV.list[i], " + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
asthma.yr5.esv.results.df <- results.df

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("crwz5 == 1 ~ ",
													ESV.list[i], " + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
crwz.yr5.esv.results.df <- results.df

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("cras_10 == 1 ~ ",
													ESV.list[i], " + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
asthma.yr10.esv.results.df <- results.df

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("crwz10 == 1 ~ ",
													ESV.list[i], " + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
crwz.yr5.esv.results.df <- results.df

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("cras_5 == 1 ~ (",
													ESV.list[i], " > 0.2) + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
asthma.yr5.esv.0.2.results.df <- results.df

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("crwz5 == 1 ~ (",
													ESV.list[i], " > 0.2) + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
crwz.yr5.esv.0.2.results.df <- results.df

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("cras_10 == 1 ~ (",
													ESV.list[i], " > 0.2) + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
asthma.yr10.esv.0.2.results.df <- results.df

results.df <- data.frame(ESV.name = ESV.list,
												 OR = numeric(length(ESV.list)),
												 lci = numeric(length(ESV.list)),
												 uci = numeric(length(ESV.list)),
												 p = numeric(length(ESV.list)),
												 stringsAsFactors=F)

for(i in 1:length(ESV.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("crwz10 == 1 ~ (",
													ESV.list[i], " > 0.2) + sex + season")),
				data = df3.well, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
crwz.yr10.esv.0.2.results.df <- results.df

subject.table.df <- premaster.df
colnames(subject.table.df)
colnames(subj.table.routine.df3)
subj.mpg.2mth.df <- subj.table.routine.df3[,1:14]
subj.mpg.2mth.df$SUBJ_ID <- rownames(subj.table.routine.df3)
colnames(subj.mpg.2mth.df) <- sub("2mth.", "", colnames(subj.mpg.2mth.df))
	
subj.mpg.2mth.df <- merge(subject.table.df, subj.mpg.2mth.df, by="SUBJ_ID", all=T)
subj.mpg.2mth.df <- subj.mpg.2mth.df %>% filter(!is.na(`Haemophilus MPG2`))

subj.table.routine.df3$SUBJ_ID <- rownames(subj.table.routine.df3)
subj.mpg.all.df <-  merge(subject.table.df, subj.table.routine.df3, by="SUBJ_ID", all=T)
	
season.for.mpg.2mth.df <- with(df3.well, data.frame(SUBJ_ID, well.sick, season))
subj.mpg.2mth.df <- merge(subj.mpg.2mth.df, season.for.mpg.2mth.df, 
													by="SUBJ_ID", all.x=F, all.y=T)


MPG.top.list <- c("Haemophilus MPG2", "Streptococcus MPG1", 
									"Haemophilus MPG1", "Moraxella MPG", 
									"Mixed MPG", "Corynebacterium MPG", 
									"Alloiococcus MPG", "Staphylococcus MPG")

results.df <- data.frame(ESV.name = MPG.top.list,
												 OR = numeric(length(MPG.top.list)),
												 lci = numeric(length(MPG.top.list)),
												 uci = numeric(length(MPG.top.list)),
												 p = numeric(length(MPG.top.list)),
												 stringsAsFactors=F)

for(i in 1:length(MPG.top.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("cras_5 == 1 ~ `",
													MPG.top.list[i], "` + sex + season")),
				data = subj.mpg.2mth.df, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
asthma.yr5.mpg.results.df <- results.df

results.df <- data.frame(ESV.name = MPG.top.list,
												 OR = numeric(length(MPG.top.list)),
												 lci = numeric(length(MPG.top.list)),
												 uci = numeric(length(MPG.top.list)),
												 p = numeric(length(MPG.top.list)),
												 stringsAsFactors=F)

for(i in 1:length(MPG.top.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("crwz5 == 1 ~ `",
													MPG.top.list[i], "` + sex + season")),
				data = subj.mpg.2mth.df, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
crwz.yr5.mpg.results.df <- results.df

results.df <- data.frame(ESV.name = MPG.top.list,
												 OR = numeric(length(MPG.top.list)),
												 lci = numeric(length(MPG.top.list)),
												 uci = numeric(length(MPG.top.list)),
												 p = numeric(length(MPG.top.list)),
												 stringsAsFactors=F)

for(i in 1:length(MPG.top.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("cras_10 == 1 ~ `",
													MPG.top.list[i], "` + sex + season")),
				data = subj.mpg.2mth.df, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
asthma.yr10.mpg.results.df <- results.df

results.df <- data.frame(ESV.name = MPG.top.list,
												 OR = numeric(length(MPG.top.list)),
												 lci = numeric(length(MPG.top.list)),
												 uci = numeric(length(MPG.top.list)),
												 p = numeric(length(MPG.top.list)),
												 stringsAsFactors=F)

for(i in 1:length(MPG.top.list)){
	glm.out <- tryCatch(
		glm(as.formula(paste0("crwz10 == 1 ~ `",
													MPG.top.list[i], "` + sex + season")),
				data = subj.mpg.2mth.df, 
				family= binomial),
		error=function(err) NA)
	if(is.na(glm.out)) {
		s <- list(coefficients = matrix(data=NA, nrow=6, ncol=5))
	} else{
		s  <- summary(glm.out)
	}	
	OR  <-  exp(s$coefficients[,1])
	lci <- exp(s$coefficients[,1] - 1.96*(s$coefficients[,2]))
	uci <- exp(s$coefficients[,1] + 1.96*(s$coefficients[,2]))
	p <- s$coefficients[,4]
	results.df$OR[i] <- signif(OR[2],2)
	results.df$lci[i] <- signif(lci[2],2)
	results.df$uci[i] <- signif(uci[2],2)
	results.df$p[i] <- signif(p[2],2)
}
crwz.yr10.mpg.results.df <- results.df

