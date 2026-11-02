# RepeatABEL Vignette
# ====

library(RepeatABEL)
#GenABEL object including IDs and marker genotypes
data(gen.data)
#Phenotype data with repeated observations
data(Phen.Data)
str(gen.data)
str(Phen.Data)

GWAS1 <- rGLS(y ~ age + sex, genabel.data = gen.data, 
							 phenotype.data = Phen.Data)

# To use with plink data: see http://forum.genabel.org/viewtopic.php?t=740
# and http://www.genabel.org/GenABEL/convert.snp.ped.html

summary(GWAS1)
plot(GWAS1)
