# Get pca for each timepoint, healthy samples
# ---------------------------------------
library(tidyr)

MPG.list

pca <- list()
pca.df.list <- list()

subj.table.healthy.df2 <- 
	with(subj.table.healthy.df,
			 as.data.frame.matrix(xtabs(prop ~ subj + age.MPG)))

for (i in 1:length(MPG.list)) {
	mpg <- MPG.list[i]
	df <- log.table.df %>% 
		filter((visittype == "Well" | visittype == "Well/Sick"),
					 uri.lri == "Well") %>%
		group_by(age.custom, subj, MPG.name) %>%
		summarize(n = n())
	with(df, as.data.frame.matrix(xtabs(n ~ subj + age.custom)))
		spread(age.custom, n)
	
	df[is.na(df)] <- 0
	mat <- as.matrix(df[,-1])
	row.names(mat) <- df$subj
	
	pca[[i]] <- prcomp(mat)
	pca.df <- data.frame(pca[[i]]$x)
	pca.df$subj <- row.names(pca.df)
	pca.df.list[[i]] <- pca.df
}