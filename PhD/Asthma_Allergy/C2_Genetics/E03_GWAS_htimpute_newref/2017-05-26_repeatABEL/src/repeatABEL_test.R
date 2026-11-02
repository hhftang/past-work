rm(list=ls())
desired.packages <- c("GenABEL", "RepeatABEL")
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

# Create tped format of snps data
snpid <- seq(1:nrow(snps))
snpid2 <- paste0("rs",snpid)
snps.tped <- data.frame(chr = rep(1, nrow(snps)),
												locus = snpid2,
												gen_dist = rep(0, nrow(snps)),
												pos = seq(1:nrow(snps)))

snps.allele.df <- as.data.frame(matrix(data = 0, nrow=nrow(snps), ncol=ncol(snps)*2))
for (i in 1:ncol(snps)){
	snps.allele.df[which(snps[,i] == 0),2*i-1] <- "A"
	snps.allele.df[which(snps[,i] == 0),2*i] <- "A"
	snps.allele.df[which(snps[,i] == 0.5),2*i-1] <- "G"
	snps.allele.df[which(snps[,i] == 0.5),2*i] <- "A"
	snps.allele.df[which(snps[,i] == 1),2*i-1] <- "G"
	snps.allele.df[which(snps[,i] == 1),2*i] <- "G"
}

snps.tped2 <- cbind(snps.tped, snps.allele.df)
write.table(snps.tped2, "data/snps.tped", row.names = F, col.names=F,
						quote=F)

# Create tfam file
snps.tfam <- data.frame(FID = seq(1:ncol(snps)),
									 IID = seq(1:ncol(snps)),
									 pat = rep(0, ncol(snps)),
									 mat = rep(0, ncol(snps)),
									 sex = rep(1, ncol(snps)),
									 pheno = rep(-9, ncol(snps)))
write.table(snps.tfam, "data/snps.tfam", row.names = F, col.names=F,
						quote=F)

# Create sex.phe file
sex.phe.df <- data.frame(id = seq(1:ncol(snps)),
													 sex = rep(1, ncol(snps)))
write.table(sex.phe.df, "data/sex.phe", row.names = F, col.names=T,
						quote=F)

start.time.1 <- Sys.time()
# Export tped to GenABEL raw format
convert.snp.tped(tpedfile="data/snps.tped", 
								 tfamfile="data/snps.tfam", 
								 outfile="data/snps.raw", strand="u", bcast=10000)

# Reload GenABEL R object
gen.data <- load.gwaa.data(phenofile = "data/sex.phe", 
													 genofile = "data/snps.raw")

end.time.1 <- Sys.time()

# Create pheno.phe file
phenos.df <- data.frame(id = seq(1:ncol(snps)),
												sex = rep(1, ncol(snps)),
												phenos)
colnames(phenos.df) <- c("id", "sex", 1,2)
phenos.m.df <- melt(phenos.df, 
										id.vars = c("id", "sex"),
										variable.name = "age",
										value.name = "y")

start.time.2 <- Sys.time()
test.gwas <- rGLS(y ~ age, genabel.data = gen.data, 
									phenotype.data = phenos.m.df)

test.gwas.total.df <- cbind(test.gwas@results, test.gwas@annotation) 
test.gwas.total.df$SNP <- rownames(test.gwas.total.df)
colnames(test.gwas.total.df)[which(colnames(test.gwas.total.df) == "Chromosome")] <- "CHR"
colnames(test.gwas.total.df)[which(colnames(test.gwas.total.df) == "Position")] <- "BP"
colnames(test.gwas.total.df)[which(colnames(test.gwas.total.df) == "P1df")] <- "P"
test.gwas.total.df$CHR <- as.numeric(as.character(test.gwas.total.df$CHR))
test.gwas.total.df <- test.gwas.total.df[order(test.gwas.total.df$P, na.last=T),]
test.gwas.total.df <- test.gwas.total.df[order(test.gwas.total.df$BP),]
end.time.2 <- Sys.time()

pdf("output/repeatABEL_man_plot.pdf")
plot(-log10(test.gwas.total.df$P))
dev.off()

pdf("output/repeatABEL_qq_plot.pdf")
qq(test.gwas.total.df$P)
dev.off()

chisq <- qchisq(1-test.gwas.total.df$P,1)
lambda <- median(chisq)/qchisq(0.5,1)
write.table(lambda, "output/repeatABEL_lambda.txt")

total.time <- (end.time.1 - start.time.1) + 
	(end.time.2 - start.time.2)
write.table(total.time, "output/repeatABEL_compute_time.txt")

results.df <- data.frame(id = seq(1:nrow(test.gwas.total.df)),
												 p = test.gwas.total.df$P)
write.table(results.df, "output/repeatABEL_results.txt",
						row.names=F, quote=F)

