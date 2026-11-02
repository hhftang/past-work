	rm(list=ls())
	source('src/longGWAS.R')
	
	desired.packages <- c("qqman")
	new.packages <- desired.packages[!(desired.packages %in% 
																		 	installed.packages()[,"Package"])]
	
	if(length(new.packages)>0) install.packages(new.packages, 
																							repos='https://cran.ms.unimelb.edu.au/')
	for (i in 1:length(desired.packages)){
		library(desired.packages[i], character.only=T)
	}
	
	
	
	# First load the SNP and phenotype files
	snps = as.matrix(read.table('data/snps.tab')) # 100 x 500 matrix of SNP by individual
	phenos = as.matrix(read.table('data/phenos.tab')) # 500 x 2 matrix of Individual by time point
	
	n = nrow(phenos)
	m = ncol(phenos)
	
	start.time.1 <- Sys.time()
	
	# Load pre-computed K matrix
	K = as.matrix(read.table('data/kinship.tab'))
	
	# GWAS
	phenos.gwas  = longGWAS.gwas(phenos, n, m, snps = snps, K = K)
	
	end.time.1 <- Sys.time()
	
	pdf("output/longgwas_man_plot.pdf")
	plot(-log10(phenos.gwas$ps))
	dev.off()
	
	pdf("output/longgwas_qq_plot.pdf")
	qq(phenos.gwas$ps)
	dev.off()
	
	chisq <- qchisq(1-phenos.gwas$ps,1)
	lambda <- median(chisq)/qchisq(0.5,1)
	write.table(lambda, "output/longgwas_lambda.txt")
	
	total.time <- end.time.1 - start.time.1
	write.table(total.time, "output/longgwas_compute_time.txt")
	
	results.df <- data.frame(id = seq(1:length(phenos.gwas$ps)),
													 	p = phenos.gwas$ps)
	write.table(results.df, "output/longgwas_results.txt",
							row.names=F, quote=F)
	
	
	
	
