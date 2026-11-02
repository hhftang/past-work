rm(list=ls())
desired.packages <- c("gee", "reshape2")
new.packages <- desired.packages[!(desired.packages %in% 
																	 	installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages, 
																						repos='https://cran.ms.unimelb.edu.au/')
for (i in 1:length(desired.packages)){
	library(desired.packages[i], character.only=T)
}

# Use test data from longGWAS.R
snps <- as.matrix(read.table('data/snps.tab')) # 100 x 500 matrix of SNP by individual
phenos <- as.matrix(read.table('data/phenos.tab')) # 500 x 2 matrix of Individual by time point

snps.t <- t(snps)
snps.df <- data.frame(snps.t)
snps.df$subjid <- seq(1,nrow(snps.df))

# Get pc's for genetic relatedness
start.time.1 <- Sys.time()
pc <- prcomp(snps, rank=7)
plot(pc$rotation[,1], pc$rotation[,2])
pc.df <- data.frame(pc$rotation)
pc.df$subjid <- seq(1,nrow(pc.df))
end.time.1 <- Sys.time()

# Get phenos
phenos.df <- data.frame(phenos)
colnames(phenos.df) <- c("T1", "T2")
phenos.df$subjid <- rownames(phenos.df)
phenos.m.df <- melt(phenos.df, id.vars="subjid", 
										variable.name="time", 
										value.name="response")

start.time.2 <- Sys.time()
all.df <- merge(pc.df, phenos.m.df, by="subjid", all=T)
all.df <- merge(all.df, snps.df, by="subjid", all=T)

# create results matrix
results.mat <- matrix(nrow=nrow(snps), ncol=1) 

for (i in 1:nrow(snps)) {
	j <- i + 10 
	test.gee <- gee(response ~ get(colnames(all.df)[j]) + time + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, 
			id=subjid, 
			data=all.df, 
			corstr="exchangeable")
	z <- coef(summary(test.gee))[2,5]
	p <- 2 * pnorm(abs(z), lower.tail = FALSE)
	results.mat[i,1] <- p
}
end.time.2 <- Sys.time()

pdf("output/gee_man_plot.pdf")
plot(-log10(results.mat))
dev.off()

pdf("output/gee_qq_plot.pdf")
qq(results.mat)
dev.off()

chisq <- qchisq(1-results.mat,1)
lambda <- median(chisq)/qchisq(0.5,1)
write.table(lambda, "output/gee_lambda.txt")

total.time <- (end.time.1 - start.time.1) + 
	(end.time.2 - start.time.2)
write.table(total.time, "output/gee_compute_time.txt")

results.df <- data.frame(id = seq(1:nrow(snps)),
												 p = results.mat)
write.table(results.df, "output/gee_results.txt",
						row.names=F, quote=F)

