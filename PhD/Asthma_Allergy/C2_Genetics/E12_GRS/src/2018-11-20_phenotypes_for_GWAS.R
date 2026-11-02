# ------------------------------------------------------------------------------
# Get GRS phenotypes for plink GWAS
# ------------------------------------------------------------------------------

load("output/2018-11-19_GRS_associations.Rdata")

# Create phenotype files
pheno.df <- read.table("data/ari_1_num.phe", header=F, stringsAsFactors = F)

# Get GRS

score_vars <- c("SCORESUM_allergic_std", "SCORESUM_asthma_std",
								"SCORESUM_adult_asthma_std", "SCORESUM_childhood_asthma_std",
								"SCORESUM_ar_std", "SCORESUM_as_std", "SCORESUM_nar_std",
								"SCORESUM_combined_std", "SCORESUM_PC1_std")
for (i in 1:length(score_vars)) {
	score_var <- score_vars[i]
	pheno.df$V3 <- master_GRS.df[match(pheno.df$V1, master_GRS.df$SUBJ_ID), 
															 match(score_var, colnames(master_GRS.df))]
	write.table(pheno.df,  paste0("output/pheno/", score_var, ".phe"), 
							quote=F, sep = "\t", row.names=F, col.names=F)
}

