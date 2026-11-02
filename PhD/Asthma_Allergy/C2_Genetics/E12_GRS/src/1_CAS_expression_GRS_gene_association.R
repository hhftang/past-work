#--------------------------------------------------
# 2018-11-07
# Aim: to identify genes associated with GRS
#  Model: Gene ~ GRS + gender + 10 genotype PCs
# 
# Data are under: ~/CAS/Analysis/GRS_expression/
#
# Gene expression: INT
# GRS: scaled (215) and subset
# Genotype PCs: first 10 using geno matrix of 135
# overlapping individuals (12.7m SNPs).
#--------------------------------------------------

# GRSs
# (1) asthma GRS:
#     /projects/Microbiome/CAS_GRS/asthma_GRS/CAS_asthma_GRS.profile
# (2) allergic disease GRS
#     /projects/Microbiome/CAS_GRS/allergic_disease_GRS/CAS_allergic_disease_GRS.profile
# (3) GRS for childhood onset asthma:
#     /projects/Microbiome/CAS_GRS/childhood_asthma/childhood_asthma.profile
# (4) GRS for adult onset asthma:
#     /projects/Microbiome/CAS_GRS/adult_asthma/adult_asthma.profile
# (5) GRS for allergic rhinitis:
#     /projects/Microbiome/CAS_GRS/AR/AR.profile
# (6) GRS for allergic sensitisation:
#     /projects/Microbiome/CAS_GRS/AS/AS.profile
# (7) GRS for non-allergic rhinitis:
#     /projects/Microbiome/CAS_GRS/NAR/NAR.profile

library(foreach)
library(MatrixEQTL)
source("~/CAS/Analysis/Scripts/CAS_My_Functions_load_data.R")

# Working directory
setwd("~/CAS/Analysis/GRS_expression/")

# Load cis eQTLs
ciseQTL <- load_cis_eQTLs()
# Significant cis eGene list
ciseGene <- list(m_Ctrl = sort(unique(ciseQTL[[1]]$gene)),
                 m_LPS = sort(unique(ciseQTL[[2]]$gene)),
                 t_Ctrl = sort(unique(ciseQTL[[3]]$gene)),
                 t_PHA = sort(unique(ciseQTL[[4]]$gene)) )

# GRS table
CASGRS <- data.frame(GRS = c("Asthma","Allergic disease","Childhood asthma","Adult asthma","Allergic rhinitis","Allergic sensitisation","Non-allergic sensitisation"),
                     filename = paste0("/projects/Microbiome/CAS_GRS/", c("asthma_GRS/CAS_asthma_GRS.profile", "allergic_disease_GRS/CAS_allergic_disease_GRS.profile", "childhood_asthma/childhood_asthma.profile", "adult_asthma/adult_asthma.profile", "AR/AR.profile", "AS/AS.profile", "NAR/NAR.profile")))
all(file.exists(CASGRS$filename))
write.table(CASGRS, "Available_CAS_GRS.txt", quote = F, sep = "\t", row.names = F)

# Read 10 genotype PCs and gender
genoPCs <- read.table("~/CAS/Genotype_Data/clean_data/PCs_10_cas_michigan_imp_filtered_135indiv.txt", header = T)
gender <- read.table("~/CAS/Genotype_Data/clean_data/cas_michigan_imp_filtered_135.fam", header = F)
if(identical(genoPCs$subjid, gender$V1)) {
  genoPCs$gender <- gender$V5    # 1 = male, 2 = female
}
cov <- as.data.frame(t(genoPCs))
write.table(cov, "Cov_gender_10_genoPCs.txt", col.names = F, quote = F, sep = "\t")

# Load covariates for all 135 individuals
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile("Cov_gender_10_genoPCs.txt")
#cvrt$columnNames

# Loop for GRSs
for(whichGRS in CASGRS$GRS) {
  # Create a folder for this GRS
  cat("    GRS for", whichGRS, "\n")
  dirname <- paste0("~/CAS/Analysis/GRS_expression/", gsub(" ", "_" ,whichGRS))
  if(!file.exists(dirname)) {dir.create(dirname)}
  setwd(dirname)
  
  # Read GRS
  GRS <- read.table(CASGRS$filename[which(CASGRS$GRS == whichGRS)], header = T)
  
  # Scale
  GRS$scaled_GRS <- scale(GRS$SCORESUM)
  t_GRS <- as.data.frame(t(GRS[,c("IID","scaled_GRS")]))
  write.table(t_GRS, "GRS_scaled.txt", quote = F, sep = "\t", col.names = F)
  
  # Load in Matrix eQTL
  GRS = SlicedData$new()
  GRS$fileDelimiter = "\t"
  GRS$fileOmitCharacters = "NA"
  GRS$fileSkipRows = 1
  GRS$fileSkipColumns = 1
  GRS$fileSliceSize = 2000   # read file in pieces of 2000 rows
  GRS$LoadFile("GRS_scaled.txt")
  
  # Four conditions
  foreach(celltype = c("m","m","t","t"), treatment = c("Ctrl","LPS","Ctrl","PHA")) %do% {
    cat("  * Working on", celltype, treatment, "...\n")
    
    # Load gene expression data
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
    gene$LoadFile(paste0("~/CAS/Analysis/eQTL_mapping/",celltype,"_",treatment,"_gene/gene_expression_",celltype,"_",treatment,".txt"))
    
    # Subject names (column names) 
    expr_colname <- sapply(gene$columnNames, function(x) {unlist(strsplit(x, split = "_"))[2]})
    
    # Keep overlapping individuals
    # Subset GRS
    GRS_curr <- GRS$Clone()
    GRS_curr$ColumnSubsample(match(expr_colname, GRS$columnNames))
    # Subest covariates
    cvrt_curr <- cvrt$Clone()
    cvrt_curr$ColumnSubsample(match(expr_colname, cvrt$columnNames))
    
    if(!all(identical(GRS_curr$columnNames, as.vector(expr_colname)), identical(cvrt_curr$columnNames, as.vector(expr_colname)))) {
      cat("  Error - Check column names!\n")
      return(NULL)
    }
    
    #----- Run MatrixeQTL on all genes -----
    me <- Matrix_eQTL_main(
      snps = GRS_curr,
      gene = gene,
      cvrt = cvrt_curr,
      pvOutputThreshold = 1,
      output_file_name = NULL,
      useModel = modelLINEAR,
      errorCovariance = numeric(),
      verbose = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)
    
    # Save association test results
    write.table(me$all$eqtls, paste0("Geneexpr_GRS_",celltype,"_",treatment,".txt"), quote = F, sep = "\t", row.names = F)
    
    
    #----- Run MatrixeQTL on cis-eGenes -----
    # Keep genes that had significant cis eQTLs
    gene$RowReorder(which(gene$GetAllRowNames() %in% ciseGene[[paste0(celltype,"_",treatment)]]))
    
    me <- Matrix_eQTL_main(
      snps = GRS_curr,
      gene = gene,
      cvrt = cvrt_curr,
      pvOutputThreshold = 1,
      output_file_name = NULL,
      useModel = modelLINEAR,
      errorCovariance = numeric(),
      verbose = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE)
    
    # Save association test results
    write.table(me$all$eqtls, paste0("CiseGene_Geneexpr_GRS_",celltype,"_",treatment,".txt"), quote = F, sep = "\t", row.names = F)
    
  } # End of 4 conditions
  
} # End of all GRSs








