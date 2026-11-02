library(ggplot2)
df <- read.table("data/rs7216389.ped", header=F, stringsAsFactors =F)

df$genotype <- with(df, paste0(V7, V8))
df$genotype <- factor(df$genotype, levels=c("CC", "TC", "TT"))

df <- df %>%
	rename(SUBJ_ID = V1,
				 rs7216389_genotype = genotype) %>%
	mutate(rs7216389_T = as.numeric(as.character(factor(rs7216389_genotype, 
																				 labels=c(0,1,2))))) %>%
	select(SUBJ_ID, rs7216389_genotype, rs7216389_T)

master_GRS.df <- merge(master_GRS.df, df, by="SUBJ_ID", all.x=T)

master_GRS.df %>%
	filter(!is.na(rs7216389_genotype)) %>%
	ggplot(aes(x = rs7216389_genotype, y = hdmige2)) + 
	geom_violin()

with(master_GRS.df, 
		 prop.table(table(rs7216389_genotype, crwz5),1))
with(master_GRS.df, 
		 prop.table(table(rs7216389_genotype, cras_5),1))

with(master_GRS.df, 
		 prop.table(table(rs7216389_genotype, crwz10),1))
with(master_GRS.df, 
		 prop.table(table(rs7216389_genotype, cras_10),1))

glm.out <- glm(crwz5 ~ rs7216389_T + sex,
		family=binomial,
		data = master_GRS.df)
summary(glm.out)

master_GRS.df$rs7216389_genotype <- NULL
master_GRS.df$rs7216389_genotype.x <- NULL
master_GRS.df$rs7216389_genotype.y <- NULL
master_GRS.df$rs7216389_T <- NULL
master_GRS.df$rs7216389_T.x <- NULL
master_GRS.df$rs7216389_T.y <- NULL
