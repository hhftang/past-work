# Check read count of repeated samples

meta.df <- sample_data(coast.biom)
meta.df$LogNo[which(meta.df$LogNo == 1369)] <- 230
meta.df$LogNo[which(meta.df$LogNo == 3116)] <- 19

repeated.logs <- meta.df$LogNo[which(duplicated(meta.df$LogNo))]


# Get total number of reads per sample
read.sum.df <- data.frame(SampleID = names(sample_sums(coast.biom)),
													read.sum = sample_sums(coast.biom), 
													row.names = NULL)
# Merge with metadata
read.sum.merge.df <- merge(read.sum.df, 
													 data.frame(sample_data(coast.biom), 
													 					 SampleID = names(sample_sums(coast.biom))), 
													 by="SampleID")


df <- read.sum.merge.df %>%
	filter(read.sum.merge.df$LogNo %in% repeated.logs,
				 type %in% "COAST") %>%
	mutate(LogNo = factor(as.character(LogNo)))
df2 <- df %>%
	group_by(LogNo) %>%
	summarize(max_read = max(read.sum),
						min_read = min(read.sum))
df$LogNo <-factor(df$LogNo,
									levels = as.character(df2[order(df2$max_read),]$LogNo))


pdf("output/2018-07-12_compare_repeated_sample_read_count.pdf", width=25, height=8)
df %>%
	ggplot(aes(x=LogNo, y=log10(read.sum))) + 
	geom_boxplot() + 
	geom_hline(yintercept=log10(4000), linetype="dashed", col="red") +
	xlab("Log Number") + 
	ylab("Total number of reads (log 10)") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

length(which(df2$max_read < 4000))
length(df2$max_read)

length(which(df2$min_read < 4000 & df2$max_read > 4000))

# I think QC before summation, as if they failed once and worked another time 
# then it might be an amplification problem, in which case the low read data may 
# be junk. Check with louise if she redid PCR for these repeat samples; 
# if she did then it would make sense to only use data from the "good" amplification.
# But if they are just reruns of the same amplicons at higher conc then it 
# would make sense to combine before QC.

meta.df <- data.frame(sample_data(coast.biom))
meta.subset.df <- meta.df[which((meta.df$LogNo %in% repeated.logs) &
																	(meta.df$type == "COAST")),]
df <- with(meta.subset.df, data.frame(sampleid = rownames(meta.subset.df),
																			specimen,
																			RunName,
																			LogNo,
																			subj))
write.csv(df, "output/repeated_samples.csv", quote=F, row.names=F)
