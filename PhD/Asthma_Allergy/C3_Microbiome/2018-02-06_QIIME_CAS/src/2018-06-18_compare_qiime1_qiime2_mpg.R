library(dplyr)
library(ggplot2)
theme_set(theme_bw())

write.csv(with(log.table.df, data.frame(LogNo = LogNo_old, MPG.name)), 
					"output/cas_qiime2_groups.csv", quote=F, row.names=F)

qiime1_mpg.df <- 
	read.csv("output/cas_qiime1_groups.csv", header=T, stringsAsFactors = F)

qiime2_mpg.df <- with(log.table.df, data.frame(LogNo = LogNo_old, MPG.name))

qiime1_mpg.df$MPG.name <- factor(qiime1_mpg.df$MPG.name)
qiime2_mpg.df$MPG.name <- factor(qiime2_mpg.df$MPG.name,
																 levels=sort(levels(qiime2_mpg.df$MPG.name)))

qiime1_2_mpg.df <- merge(qiime1_mpg.df, qiime2_mpg.df, by="LogNo", all=T)
df <- as.data.frame(with(qiime1_2_mpg.df, table(MPG.name.x, MPG.name.y,
																								useNA="ifany")))

df <- df %>% rename(qiime1.mpg.name = MPG.name.x, 
							qiime2.mpg.name = MPG.name.y)

pdf("output/2018-06-18_qiime1_2_mpg_count.pdf", width=10, height=8)
ggplot(df, aes(x = qiime1.mpg.name, y= qiime2.mpg.name)) + 
	geom_tile(aes(fill = Freq)) + 
	theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
	scale_fill_gradientn(colours = seqPalette(100),
											 name = "Frequency") + 
	labs(x = "QIIME1 MPG", y = "QIIME2 MPG")
dev.off()
	
df2 <- as.data.frame(with(qiime1_2_mpg.df, 
												 prop.table(table(MPG.name.x, MPG.name.y, 
												 								 useNA="ifany"), 1)))

df2 <- df2 %>% rename(qiime1.mpg.name = MPG.name.x, 
										qiime2.mpg.name = MPG.name.y,
										prop = Freq)

pdf("output/2018-06-18_qiime1_2_mpg_prop.pdf",width=10, height=8)
ggplot(df2, aes(x = qiime1.mpg.name, y= qiime2.mpg.name)) + 
	geom_tile(aes(fill = prop)) + 
	theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
	scale_fill_gradientn(colours = seqPalette(100),
											 name = "Prop. QIIME1 MPG") + 
	labs(x = "QIIME1 MPG", y = "QIIME2 MPG")
dev.off()

save.image("output/2018-06-18_compare_qiime1_qiime2_mpg.Rdata")
