
# This file contains some very basic usage examples for longGWAS.R.
# Nick Furlotte (nfurlott@cs.ucla.edu)

source('src/longGWAS.R')

# First load the SNP and phenotype files
snps = as.matrix(read.table('data/snps.tab')) # 100 x 500 matrix of SNP by individual
phenos = as.matrix(read.table('data/phenos.tab')) # 500 x 2 matrix of Individual by time point

n = nrow(phenos)
m = ncol(phenos)

# -------------------    THE KINSHIP MATRIX 
# Load pre-computed K matrix
K = as.matrix(read.table('data/kinship.tab'))

# To compute your own load a full SNP matrix and then call myFastLMM.kinship(snps)
# You can try it with the given snp matrix, but since it is just a small number of 
# snps the kinship matrix will not be accurate. NOTE: This calculation may take some time if the nuber
# of snps is large.

# -------------------    GWAS SCAN
# The simpilest GWAS scan can be accomplished with the following, in which the function will
# calculate the kinship from the input SNPs (this case snps.full). SNPs should be encoded with 0, 0.5 and 1.
#
# phenos.gwas  = longGWAS.gwas(phenos, n, m, snps = snps.full)

# Or we can pass in the pre-computed kinship, which is what we want in this case.
# Note this snp set has been 'cleaned' (NAs have been replaced with minor allele frequency)
# so it will not work for the kinship calculation.

phenos.gwas  = longGWAS.gwas(phenos, n, m, snps = snps, K = K)

# This takes a couple of minutes to run because it has to estimate variance components.
# The function returns a few items.  You can see what by executing names(phenos.gwas)
# The element phenos.gwas$vg is the variance component for the genetic effect, phenos.gwas$ve 
# is the variance component for the environmental effect and phenos.gwas$vep is the component for error.
# The element phenos.gwas$ps gives the p-values for each SNP.
# If you are only interested in the variance components, for example, you only want to estimate heritability,
# then you can add the parameter 'fitOnly=True' to the function and it will only return these values.  
# An estimate of heritability is also returned in phenos.gwas$pseudoH.

# An important aspect of the longitudinal model is the calculation of the D matrix.
# We include a very basic method with longGWAS.estimateDMatrix(...).  
# However, you can create your own function with similar arugments to ours and pass it to .gwas 
# with the parameter D.estimator.
# If you do not specify the D parameter or a function to compute it, the .gwas function will
# just estimate it from the data using the default estimator.

# -------------------    Predicting Random Effects
# The longitudinal model assumes that a measurement at a given time point can be explained
# by the mean effect, background genetic effect and residual error.  
# Typically, we don't look directly at the effects due to genetic background and error 
# as we are only interested in performing a hypothesis test on the main effects.
# We show in the paper that having predictions of the random effect values (realized values)
# might be interesting.  Here is how to get those.

# I am passing in the varComps parameter because we just computed them.  If it is 
# not specified then this function will run the estimation for you.

phenos.pred = longGWAS.predictRandomEffects(phenos, n, m, K, varComps = phenos.gwas)

# -------------------    Calculating Analytical Power 
# In the paper we introduce a way to compute the power for a given SNP with a given effect size.
# Typically, this is done through simulations which can take a lot of time, so the ability to 
# calculate power directly is really nice.  This is useful for example to evaluate the potential
# power that can be achieved by certain experimental designs.

# First get a sample D matrix
E = as.matrix(read.table('data/D.1.tab'))
D = kronecker(diag(n), E)

# get some SNP and set an effect
snp = snps[50,]
snpEffect = 0.8
# Set the variance components 
v = list(vg=0.1, ve=0.8, vep=0.1) 
p = longGWAS.analyticalPower(n, m, snp, snpEffect,  K = K, D = D, varComps = v)







