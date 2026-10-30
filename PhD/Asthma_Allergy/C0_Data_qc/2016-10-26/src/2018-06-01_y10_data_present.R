wimmunol_10.df
polyic_y4_10.df

library(GMD)
library(dplyr)

df0 <- merge(data.frame(SUBJ_ID = premaster.df$SUBJ_ID), immunol_10.df, by="SUBJ_ID", all=T)
row.names(df0) <- df0$SUBJ_ID
df0 <- df0[order(as.character(df0$SUBJ_ID)),]
df <- data.frame(!is.na(df0))
heatmap.3(data.matrix(df),
					diss = FALSE, 
					Rowv = FALSE,
					Colv = TRUE, 
					dendrogram = "none", 
					trace="none", 
					density.info="none",
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					keysize = 1, mapsize = 9, 
					main = ""
)

colnames.na.pattern <- function(df, var) {
	colnames(df)[which(apply(df, 2, function(x) all.equal(x, var,
																												check.attributes=F) == T))]
}

A <- colnames.na.pattern(df, df$ae10)
B <- colnames.na.pattern(df, df$catspt10)
C <- colnames.na.pattern(df, df$DOB)
D <- colnames.na.pattern(df, df$ht10)

results1.df <- with(df, data.frame(SUBJ_ID = row.names(df), 
																	 A = as.numeric(ae10),
																	 B = as.numeric(catspt10),
																	 C = as.numeric(DOB),
																	 D = as.numeric(ht10)))

df2 <- merge(data.frame(SUBJ_ID = premaster.df$SUBJ_ID), 
						 polyic_y4_10.df[,c(1:3, 20:27)], by="SUBJ_ID", all=T)
row.names(df2) <- df2$SUBJ_ID
df2 <- df2 %>% rename("y10_icip" = "X.y10_icip")
df2 <- data.frame(!is.na(df2))
heatmap.3(data.matrix(df2),
					diss = FALSE, 
					Rowv = FALSE,
					Colv = TRUE, 
					dendrogram = "none", 
					trace="none", 
					density.info="none",
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					keysize = 1, mapsize = 9, 
					main = ""
)

E <- colnames.na.pattern(df2, df2$y10_ic10)
G <- colnames.na.pattern(df2, df2$DOB)

results2.df <- with(df2, data.frame(SUBJ_ID = row.names(df2),
																		E = as.numeric(y10_ic10),
																		G = as.numeric(DOB)))

lungfunction_y10.df <- read.csv("data/10yr LFT and SPT.csv", header=T, stringsAsFactors = F)
lungfunction_y10.df$SUBJ_ID_part <- formatC(lungfunction_y10.df$idno, width=4, format="d", flag="0")

df3 <- merge(data.frame(SUBJ_ID = premaster.df$SUBJ_ID,
												SUBJ_ID_part = substring(premaster.df$SUBJ_ID, 0, 4)), 
												lungfunction_y10.df, by="SUBJ_ID_part", all=T)
row.names(df3) <- df3$SUBJ_ID
df3 <- data.frame(!is.na(df3))

heatmap.3(data.matrix(df3),
					diss = FALSE, 
					Rowv = FALSE,
					Colv = TRUE, 
					dendrogram = "none", 
					trace="none", 
					density.info="none",
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					keysize = 1, mapsize = 9, 
					main = ""
)

H <- colnames.na.pattern(df3, df3$DOB)
I <- colnames.na.pattern(df3, df3$FEV1pred_A)
J <- colnames.na.pattern(df3, df3$Heightcms)
K <- colnames.na.pattern(df3, df3$FEV12)
L <- colnames.na.pattern(df3, df3$FEV1B2)

results3.df <- with(df3, data.frame(SUBJ_ID = row.names(df3), 
										 H = as.numeric(DOB),
										 I = as.numeric(FEV1pred_A),
										 J = as.numeric(Heightcms),
										 K = as.numeric(FEV12),
										 L = as.numeric(FEV1B2),
										 Age = as.numeric(Age),
										 Gender = as.numeric(gender),
										 Asthma = as.numeric(Asthma)))

results.df <- merge(results1.df, results2.df, by="SUBJ_ID", all=T)
results.df <- merge(results.df, results3.df , by="SUBJ_ID", all=T)
row.names(results.df) <- results.df$SUBJ_ID
heatmap.3(results.df[,-1],
					diss = FALSE, 
					Rowv = FALSE,
					Colv = FALSE, 
					dendrogram = "none", 
					trace="none", 
					density.info="none",
					color.FUN = colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF","lightgrey"))), 
					keysize = 1, mapsize = 9, 
					main = ""
)

write.csv(results.df, "output/2018-06-01_y10_data_present.csv", quote=F)
