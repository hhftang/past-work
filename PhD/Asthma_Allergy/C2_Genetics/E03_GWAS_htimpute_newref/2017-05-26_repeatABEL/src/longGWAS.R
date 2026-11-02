
# longGWAS - Longitudinal Genome-wide Association
# This software was written by Nicholas A. Furlotte (nfurlott@cs.ucla.edu)
# and implements the methods found in the Genetic Epi paper found here: (not published just yet)
# Naturally, there is no warranty and the author cannot be held responsible for 
# any damage or anything strange.

# You are invited to use this code for whatever you wish and to alter it as you see fit and are 
# asked to simply give credit to the author when using it.

# The basic idea behind the paper is that by utlizing phenotype measurements from multiple time points for each indivivdual in a genome-wide association study you can achieve increased statistical power when performing association.
# This .R file provides 3 main tools
#
# 1) Longitudinal genome-wide association (gwas) -  longGWAS.gwas(...)
# 2) Prediction of Random effects - longGWAS.predictRandomEffects(...)
# 3) Calculation of analytical power - longGWAS.analyticalPower(...)
#
# Each function takes very similar arguments.  At a bare minimum you will need a phenotype matrix where the individuals in the study are put on the rows and the time point values for their phenotypes are on the columns. 
# You will also need a SNP matrix, which may be used to calculate a genetic relatedness matrix or kinship.
# However, since the kinship calculation may take some time, you are encouraged to use the function myFastLMM.kinship(...) to calculate it and then save it to a file.

# Below is a description of the variables used within this module.  Each function has its own set of pointers on how to use it.  I am working on better documentation now and it will be available soon.  Please see R/example.R for a quick intro to each of the important functions.

# LIST OF VARIABLE NAMEs and DESCRIPTIONS - these are common throughout the code, so I just summarize them here
#  Y - an n x m phenotype matrix with individuals on the rows and timepoints on the columns
#  n - number of individuals
#  m - number of time points

#  snps - p x n SNP matrix with snps on the rows and individuals on the columns.  Heterozygous minor allele = 0.5, homozygous minor allele=1 and homozygous major allele=0

#  K - n x n kinship matrix - if NULL it will be calculated with myFastLMM.kinship

#  D - the nm x nm environmental covariance matrix - create the matrix as if there were no missing values in Y
#      If NULL then the function will calculate it.

#  D.estimator - if D is not provided then D.estimator is a function that will estimate it.
#                If NULL the default is used.

#  X0 - the nm x q covariate matrix, where q is the number of covariates.  
#        Details on how to build this matrix are found in longGWAS.buildCovariance
#        If NULL then a global mean is fit and that is all.

# snpsCLEAN - indicates if the SNP matrix has already been run through the cleaning function

# varComps  - a list of variance components that have been pre-computed.
#             must contain: vg - genetic variance, ve - environmental variance and vep - residual error variance.  In some cases such as in longGWAS.gwas it must also contain delta and w.
#             All of these values are returned when running longGWAS.gwas by default.

# fitOnly - flag indicating if you ONLY want to find the optimal variance components.
#           If this flag is set to true then snps can be null and no GWAS scan is performed
# simpleSearch - The default search technique for identifying optimal variance components is fairly simple.
#                The more elegant one has not been implemented here yet.

# NOTE: Please remove individuals with no observations (ie. all values are missing for that individual).  The missing value removal
# process does not currently handle that case and it has the potential to cause some strange behavior.

DEBUG = FALSE

# Perform a genome scan and return p-values for each SNP.  If fitOnly=TRUE, it just returns the optimal variance components and associated parameters.
# NOTE: For speed, it makes a lot of sense to pre-compute the Kinship matrix and save it in a file.
#  If there are many SNPs it might also make sense to pre-compute the variance components and then run all SNPs in a batch (that functionality is not provided here).
longGWAS.gwas <- function(Y, n, m, 
                          snps = NULL,
                          K = NULL,
                          D = NULL,
                          D.estimator = NULL,
                          X0 = NULL,
                          STEP = 0.1,
                          simpleSearch = TRUE,
                          snpsCLEAN = FALSE,
                          varComps = NULL,
                          fitOnly = FALSE){
    stopifnot(n == nrow(Y))
    stopifnot(m == ncol(Y))
    N = n * m

    if(!simpleSearch){ 
      cat("Only simpleSearch is available at this time\n")
      return()
    }

    if(!fitOnly){
      if(is.null(dim(snps)) && length(snps) == n){ snps = matrix(snps,nrow=1)  }
      stopifnot(n == ncol(snps))

      if(snpsCLEAN){ mafs = apply(snps,1,mean); snps.clean = list(snps=snps,mafs=mafs) }
      else{ snps.clean = myFastLMM.cleanSNPs(snps) }
      nsnps = nrow(snps.clean$snps)
    }
    
    if(is.null(K)){
       if(is.null(snps)){
          cat("ERROR: Cannot build kinship without SNP matrix\n")
          return(NA)
       }
       cat('Calculating the ', n,'x',n,' kinship matrix with myFastLMM.kinship\n')
       K = myFastLMM.kinship(snps)
    }
    if(is.null(X0)){
       X = matrix(1,N,1)
    }
    else{ X = X0 }
    q = ncol(X)

    stopifnot(nrow(X) == N)

    if(is.null(D)){
      if(m > 1){
        D = longGWAS.estimateDMatrix(Y,n,m,K,X,D.estimator)
      }
      if( m == 1){ D = diag(N) }
    }
    
    # Double check that the D estimating function is doing the right thing
    stopifnot(nrow(D) == N)
    stopifnot(ncol(D) == N)

    # Fit the null model to find variance components

    # vec(Y') - results in a vector such that the first m elements are the time points for indiv 1.
    #           The next m elements are for indiv 2 and so on...
    Y.c = c(t(Y))
    # ZKZ'
    Kz = kronecker(K, matrix(1,m,m))

    Y.missing = which(is.na(Y.c))
    if(length(Y.missing) > 0){
       cat("Removing a total of ", length(Y.missing)," entries with missing values from the Y matrix\n")
       Y.c = Y.c[-Y.missing]
       X = X[-Y.missing,]
       D  = D[-Y.missing,-Y.missing]
       Kz = Kz[-Y.missing, -Y.missing]
       N = N - length(Y.missing)
    }
    
    # END of all of the Clean up Stuff - Now we fit the model

    if(is.null(varComps)){
       LLs = c()
       opts = list()
       cat("Finding optimal variance components with step size of ", STEP, "\n")
       WS = seq(0.01, 0.99, by=STEP)
       for(w in WS){ 
         cat("\t Step w = ",w,"\n")
         KzS = w * Kz + (1.0 - w) * D
         Y.c.t= myFastLMM.transform(Y.c,X,N,q,KzS)
         Y.c.t.o = myFastLMM.optimizeDeltaREML(Y.c.t$Ys, Y.c.t$Xs, N, q, Y.c.t$H, Y.c.t$S,KzS)
         LLs = append(LLs,Y.c.t.o$LL)
         opts[[length(opts) + 1]] <-  Y.c.t.o

         if(DEBUG){
            cat("LL = ",Y.c.t.o$LL,"\n")
            cat("vg = ",Y.c.t.o$vg*w,"\n")
            cat("ve = ",Y.c.t.o$vg*(1.0 - w),"\n")
            cat("vp = ",Y.c.t.o$ve,"\n")
            cat("delta = ",Y.c.t.o$delta,"\n")
         }
       }

       if(! simpleSearch ){
         longGWAS.error('This search technique has not been implemented yet')
         return(NA)
       }
       else{
         x = which.max(LLs)
         if(is.na(x)){
            cat("Maximum Likelihood Procedure was unsuccessful.  Please verify data.\n")
            return(NA)
         }

         maxLL = LLs[x]
         maxW = WS[x]
         VG = opts[[x]]$vg * maxW
         VE = opts[[x]]$vg * (1.0 - maxW)
         VEp = opts[[x]]$ve
         pseudoH = VG/(VG+VE+VEp)
         delta = opts[[x]]$delta

       }

       cat("Identified maximum likelihood parameters with maxLL = ",maxLL,": \n")
       cat("\t w = ",maxW,"\n")
       cat("\t genetic variance = ",VG,"\n")
       cat("\t environmental variance = ",VE,"\n")
       cat("\t error variance = ",VEp,"\n")
       cat("Pseudo-heritability is given by ", pseudoH, "\n")

       if(fitOnly){
          return(list(vg=VG,ve=VE,vep=VEp,delta=delta,w=maxW,pseudoH=pseudoH)) 
       }
    }
    else{  # varComps is NULL
         if(sum(is.na(match(c("vg","ve","vep","delta","w"),names(varComps))))){
            cat("When providing the varComps argument please include these elements:\n")
            cat("\t vg, ve, vep, delta, w\n")
            return(NA)
         }
         maxW = varComps$w
         VG = varComps$vg
         VE = varComps$ve
         VEp = varComps$vep
         pseudoH = VG/(VG+VE+VEp)
         delta = varComps$delta
    }

    # This is annoying because we have to perform the same transformation twice when we have performed the search.
    # It is either this or we save everything and use a lot of memory.
    KzS = maxW * Kz + (1.0 - maxW) * D
    Y.c.t= myFastLMM.transform(Y.c,X,N,q,KzS)

    # BEGIN the GWAS Scan
    cat("Starting GWAS Scan....\n")

    ps = vector(length=nsnps)
    ts = vector(length=nsnps)

    Ht = t(Y.c.t$H)

    for(i in 1:nsnps){
       if(DEBUG){ cat("SNP ",i,"\n") } 
       snp = kronecker(matrix(snps.clean$snps[i,],ncol=1), matrix(1,nrow=m,ncol=1))
       if(length(Y.missing) > 0) { snp = snp[-Y.missing,] }
       snp.t = Ht %*% snp
       # Check for linear dependencies in the covariates
       if( var(snp.t) == 0){ ps[i] <- NA; ts[i] <- NA; next }
       Xs = cbind(Y.c.t$Xs, snp.t)
       if(det(crossprod(Xs,Xs)) <= 1e-7){ ps[i] <- NA; ts[i] <- NA; next}

       R = myFastLMM.BetadBeta(Y.c.t$Ys,Xs,N,q+1,Y.c.t$H,Y.c.t$S,delta) 
       beta = R$beta
       sigma = myFastLMM.sigma(Y.c.t$Ys,Xs,N,q+1,Y.c.t$H,Y.c.t$S,delta,beta=beta)[1,1]
       sigma = (sigma*N)/(N-(q+1))
       stat = beta[q+1] / sqrt(R$stderrSQ[q+1,q+1] * sigma)
       pss = 2*pt(abs(stat), N-(q+1),lower.tail=F)

       ts[i] <- stat
       ps[i] <- pss
    }

    return(list(ts=ts,ps=ps,vg=VG,ve=VE,vep=VEp,delta=delta,w=maxW,pseudoH=pseudoH)) 
}

# This function obtains predictions of the random effects.
# You can pass K.inv from the result of solve(K) to save time if you have it.

# NOTE: please be sure to remove individuals with no valid measurements (all NAs).  These will be a problem
# here for computing the genetic contributions.
longGWAS.predictRandomEffects <- function(Y, n, m, K, D=NULL, X0=NULL, K.inv=NULL, varComps=NULL, D.estimator=NULL){

    stopifnot(n == nrow(Y))
    stopifnot(m == ncol(Y))
    N = n * m

    if(is.null(X0)){
       X = matrix(1,N,1)
    }
    else{ X = X0 }
    q = ncol(X)
    stopifnot(nrow(X) == N)

    if(is.null(D)){
      D = longGWAS.estimateDMatrix(Y,n,m,K,X,D.estimator)
    }

    stopifnot(nrow(D) == N)
    stopifnot(ncol(D) == N)

    Y.c = c(t(Y))
    Z = kronecker(diag(n),matrix(1,m,1))

    Y.missing = which(is.na(Y.c))
    if(length(Y.missing) > 0){
       Y.c = Y.c[-Y.missing]
       X = X[-Y.missing,]
       Z = Z[-Y.missing,]
       D  = D[-Y.missing,-Y.missing]
       N = N - length(Y.missing)
    }
    
    if(is.null(varComps)){
       cat("Obtaining variance component estimates...\n")
       varComps = longGWAS.gwas(Y,n,m,D=D,X0=X0,K=K,fitOnly=TRUE)
    }
    
    vg = varComps$vg
    ve = varComps$ve
    vep = varComps$vep
       
    if(is.null(K.inv)){ K.inv = solve(K) }
    Xt = t(X)
    Zt = t(Z)
    W = diag(nrow(D))
    Wt = t(W) 

    cat("Inverting D\n")
    D.inv = solve(D)
   
    alpha1 = vep/vg
    alpha2= vep/ve

    # Obtain the predictions

    # Left matrix
    cat("Obtaining Left Hand Matrix\n")
    betaRow =  cbind(Xt%*%X, Xt%*%Z, Xt%*%W)
    gRow    =  cbind(Zt%*%X, Zt%*%Z+K.inv*alpha1, Zt%*%W)
    eRow   =  cbind(Wt%*%X, Wt%*%Z, Wt%*%W+D.inv*alpha2)

    cat("Solving Left hand matrix\n")
    M = solve(rbind(betaRow,gRow,eRow))
 
    # Right Matrix
    Q = rbind(Xt%*%Y.c, Zt%*%Y.c, Wt%*%Y.c)

    # Solve
    cat("Solving...\n")
    S = as.double(M %*% Q) 

    betaHat = S[1:q]
    gHat = S[(q+1):(n+q)] # Assuming that you have one value for each individual
    eHat = S[(n+q+1):length(S)]

    gHat = matrix(gHat,ncol=1)
    eHatX = matrix(NA,nrow=N,ncol=1)
    if(length(Y.missing)){ eHatX[(1:(n*m))[-Y.missing]] <- eHat }
    else{ eHatX = matrix(eHat,ncol=1) }
    eHat = matrix(eHatX,nrow=n,ncol=m,byrow=T)

    return(list(betaHat=betaHat,gHat=gHat,eHat=eHat,S=S,varComps=varComps))
}

# The following function follows almost the same structure as longGWAS.gwas except for a few differences.
# It assumes that Y is an n x 2 matrix where the first vector of n 
# is the measured phenotype under condition 1 and the second is under condition two for n individuals.
# When this is the case, we can fit a GxE model as a special case of the longitudinal model.
# This model has an additional fixed effect for the gene-by-environment (GxE) interaction
# as well as an additional random effect that accounts for the difference in variance between the two conditions.
longGWAS.gwas_GxE <- function(Y, n,  
                          snps = NULL,
                          K = NULL,
                          X0 = NULL,
                          STEP = 0.1,
                          simpleSearch = TRUE,
                          snpsCLEAN = FALSE,
                          varComps = NULL,
                          fitOnly = FALSE){

    stopifnot(n == nrow(Y))
    # Only support two conditions now
    stopifnot(2 == ncol(Y))
    m = 2
    N = n * m

    if(!simpleSearch){ 
      cat("Only simpleSearch is available at this time\n")
      return()
    }

    if(!fitOnly){
      if(is.null(dim(snps)) && length(snps) == n){ snps = matrix(snps,nrow=1)  }
      stopifnot(n == ncol(snps))

      if(snpsCLEAN){ mafs = apply(snps,1,mean); snps.clean = list(snps=snps,mafs=mafs) }
      else{ snps.clean = myFastLMM.cleanSNPs(snps) }
      nsnps = nrow(snps.clean$snps)
    }
    
    if(is.null(K)){
       if(is.null(snps)){
          cat("ERROR: Cannot build kinship without SNP matrix\n")
          return(NA)
       }
       cat('Calculating the ', n,'x',n,' kinship matrix with myFastLMM.kinship\n')
       K = myFastLMM.kinship(snps)
    }
    if(is.null(X0)){
       # default X is a matrix with a global mean and treatement effect
       X = cbind(matrix(1,N,1),c(rep(0,n),rep(1,n)))
    }
    else{ X = X0 }
    q = ncol(X)

    stopifnot(nrow(X) == N)

    # This structure means that between conditions have zero covariance in this component
    # and within condition have covariance proportional to the kinship coefficient.
    D = kronecker(K, matrix(c(1,0,0,1), 2,2))

    # Fit the null model to find variance components

    # vec(Y') - results in a vector such that the first m elements are the time points for indiv 1.
    #           The next m elements are for indiv 2 and so on...
    Y.c = c(t(Y))
    # ZKZ'
    Kz = kronecker(K, matrix(1,m,m))

    Y.missing = which(is.na(Y.c))
    if(length(Y.missing) > 0){
       cat("Removing a total of ", length(Y.missing)," entries with missing values from the Y matrix\n")
       Y.c = Y.c[-Y.missing]
       X = X[-Y.missing,]
       D  = D[-Y.missing,-Y.missing]
       Kz = Kz[-Y.missing, -Y.missing]
       N = N - length(Y.missing)
    }
    
    # END of all of the Clean up Stuff - Now we fit the model

    if(is.null(varComps)){
       LLs = c()
       opts = list()
       cat("Finding optimal variance components with step size of ", STEP, "\n")
       WS = seq(0.01, 0.99, by=STEP)
       for(w in WS){ 
         cat("\t Step w = ",w,"\n")
         KzS = w * Kz + (1.0 - w) * D
         Y.c.t= myFastLMM.transform(Y.c,X,N,q,KzS)
         Y.c.t.o = myFastLMM.optimizeDeltaREML(Y.c.t$Ys, Y.c.t$Xs, N, q, Y.c.t$H, Y.c.t$S,KzS)
         LLs = append(LLs,Y.c.t.o$LL)
         opts[[length(opts) + 1]] <-  Y.c.t.o

         if(DEBUG){
            cat("LL = ",Y.c.t.o$LL,"\n")
            cat("vg = ",Y.c.t.o$vg*w,"\n")
            cat("ve = ",Y.c.t.o$vg*(1.0 - w),"\n")
            cat("vp = ",Y.c.t.o$ve,"\n")
            cat("delta = ",Y.c.t.o$delta,"\n")
         }
       }

       if(! simpleSearch ){
         longGWAS.error('This search technique has not been implemented yet')
         return(NA)
       }
       else{
         x = which.max(LLs)
         if(is.na(x)){
            cat("Maximum Likelihood Procedure was unsuccessful.  Please verify data.\n")
            return(NA)
         }

         maxLL = LLs[x]
         maxW = WS[x]
         VG = opts[[x]]$vg * maxW
         VE = opts[[x]]$vg * (1.0 - maxW)
         VEp = opts[[x]]$ve
         pseudoH = VG/(VG+VE+VEp)
         delta = opts[[x]]$delta

       }

       cat("Identified maximum likelihood parameters with maxLL = ",maxLL,": \n")
       cat("\t w = ",maxW,"\n")
       cat("\t genetic variance = ",VG,"\n")
       cat("\t environmental variance = ",VE,"\n")
       cat("\t error variance = ",VEp,"\n")
       cat("Pseudo-heritability is given by ", pseudoH, "\n")

       if(fitOnly){
          return(list(vg=VG,ve=VE,vep=VEp,delta=delta,w=maxW,pseudoH=pseudoH)) 
       }
    }
    else{  # varComps is NULL
         if(sum(is.na(match(c("vg","ve","vep","delta","w"),names(varComps))))){
            cat("When providing the varComps argument please include these elements:\n")
            cat("\t vg, ve, vep, delta, w\n")
            return(NA)
         }
         maxW = varComps$w
         VG = varComps$vg
         VE = varComps$ve
         VEp = varComps$vep
         pseudoH = VG/(VG+VE+VEp)
         delta = varComps$delta
    }

    # This is annoying because we have to perform the same transformation twice when we have performed the search.
    # It is either this or we save everything and use a lot of memory.
    KzS = maxW * Kz + (1.0 - maxW) * D
    Y.c.t= myFastLMM.transform(Y.c,X,N,q,KzS)
    if(sum(Y.c.t$S < 0)){
      Y.c.t$S[which(Y.c.t$S < 0)] <- 1e-6
    }

    # BEGIN the GWAS Scan
    cat("Starting GWAS Scan....\n")

    ps = vector(length=nsnps)
    ts = vector(length=nsnps)
    ps.gxe = vector(length=nsnps)
    ts.gxe= vector(length=nsnps)

    Ht = t(Y.c.t$H)
    treatmentVector = matrix(c(rep(1,n),rep(2,n)),ncol=1)
    q = q + 1

    for(i in 1:nsnps){
       cat('at snp ',i,"\n")
       snp = kronecker(matrix(snps.clean$snps[i,],ncol=1), matrix(1,nrow=m,ncol=1))
       snp.gxe = treatmentVector * snp
       
       if(length(Y.missing) > 0) { snp = snp[-Y.missing,]; snp.gxe = snp.gxe[-Y.missing,] }
       snp.t = Ht %*% snp
       snp.gxe.t = Ht %*% snp.gxe
       if( var(snp.t) == 0 || var(snp.gxe.t) == 0){ ps[i] <- NA; ts[i] <- NA; ps.gxe[i] <- NA; ts.gxe[i] <- NA; next }
       Xs = cbind(Y.c.t$Xs, snp.gxe.t, snp.t)
       XsXs = crossprod(Xs,Xs)
       XsXs.e = eigen(XsXs)
       if(det(XsXs) <= 1e-7){ ps[i] <- NA; ts[i] <- NA; ps.gxe[i] <- NA; ts.gxe[i] <- NA; next }
       if(prod(XsXs.e$values) <= 1e-7){ ps[i] <- NA; ts[i] <- NA; ps.gxe[i] <- NA; ts.gxe[i] <- NA; next }
       if(min(XsXs.e$values) <= 1e-10){ ps[i] <- NA; ts[i] <- NA; ps.gxe[i] <- NA; ts.gxe[i] <- NA; next }

       R = myFastLMM.BetadBeta(Y.c.t$Ys,Xs,N,q+1,Y.c.t$H,Y.c.t$S,delta)
       #tryCatch(R = myFastLMM.BetadBeta(Y.c.t$Ys,Xs,N,q+1,Y.c.t$H,Y.c.t$S,delta), finally={ps[i] <- NA; ts[i] <- NA; ps.gxe[i] <- NA; ts.gxe[i] <- NA; next})
       beta = R$beta
       sigma = myFastLMM.sigma(Y.c.t$Ys,Xs,N,q+1,Y.c.t$H,Y.c.t$S,delta,beta=beta)[1,1]
       sigma = (sigma*N)/(N-(q+1))

       stat = beta[q+1] / sqrt(R$stderrSQ[q+1,q+1] * sigma)
       pss = 2*pt(abs(stat), N-(q+1),lower.tail=F)
       ts[i] <- stat
       ps[i] <- pss

       stat = beta[q] / sqrt(R$stderrSQ[q,q] * sigma)
       pss = 2*pt(abs(stat), N-(q+1),lower.tail=F)
       ts.gxe[i] <- stat
       ps.gxe[i] <- pss
    }

    return(list(ts=ts,ps=ps,ps.gxe=ps.gxe,ts.gxe=ts.gxe,vg=VG,ve=VE,vep=VEp,delta=delta,w=maxW,pseudoH=pseudoH)) 
}

longGWAS.estimateDMatrix <- function(Y, n, m, K, X, D.estimator=NULL){

   if(! is.null(D.estimator) ){  return(D.estimator(Y,n,m,K,X)) }
   
   Ehat = cor(Y,use='pairwise.complete')
   D = kronecker(diag(n), Ehat)

   return(D)
}

# This function fits each measurement separately and then adjusts by the predicted genetic background.
# The correlation between measurements is then calculated in the transformed space.
# This sometimes works better if the genetic effect is large.
longGWAS.estimateDMatrix_conditionOnGenetics <- function(Y, n, m, K, X){

      Yt = matrix(NA,nrow=n, ncol=m)
      K.e = eigen(K)
      # This function only fits a global mean.
      # It should really take into consideration other covariates if they exist.
      # For the sake of simplicity I am not going to do that here, but using this basic function 
      # as a framework a custom function is easy to write.
     
      X = matrix(1,n,1)
      q = ncol(X)
      for( i in 1:m){
        x = which(is.na(Y[,i])) 
        if(length(x) > 0){ y = Y[-x,i]; Xs = X[-x,]; Ks = K[-x,-x]; Ks.e=NULL}
        else{ y = Y[,i]; Xs = X; Ks = K; Ks.e = K.e; }
        y = matrix(y, ncol=1)
      
        Y.t = myFastLMM.transform(y,Xs,n-length(x),q,Ks,Ks.e)
        Y.t.o = myFastLMM.optimizeDeltaREML(Y.t$Ys, Y.t$Xs, n-length(x), q, Y.t$H, Y.t$S,Ks)
        yt <- t(Y.t$K.eigen$vectors %*% diag(1/sqrt(Y.t$K.eigen$values*Y.t.o$vg + Y.t.o$ve))) %*% (y - Xs %*% Y.t.o$beta)
        if(length(x) > 0){ Yt[(1:n)[-x],i] <- yt  } 
        else{ Yt[,i] <- yt }
      }
      Ehat = cov(Yt, use='pairwise.complete')
      D = kronecker(diag(n), Ehat)

      return(D)
}

longGWAS.estimateDMatrix_conditionOnGenetics_iterative <- function(Y, n, m, K, X){

    # This function fits the full model on the measurements beginning with the observed covariance between
    # measurements as the D matrix.  The predicted genetic contributions are "regressed" from the phenotypes
    # and then the procedure is repeated.  
    # This is an alternative approach to the conditionOnGentics function and sometimes may work better.

   # COMING SOON...

}

longGWAS.buildCovariance <- function(n,m,X){
   
   # For the most basic case you will have an n x q covariate matrix X that will be common
   # to each time point.  With m time points you will then create the full covariate matrix
   # by saying X.full = kronecker(X, matrix(1,m,1)).  That is, for each element in X you want to 
   # represent it multiple times in the X.full because there are multiple measurements for each individual 
   # and we assume that the first m measurements are for individual 1 and the second m for individual 2 and so on.
   
   # To create a covariate matrix such that covariates are different over time, it is a little more difficult.
   # You want the first line of your covariate matrix to represent the covariates for individual 1 at time point 1
   # and the second line to represent inividual 1 at time point 2 and so on.  For example, lets say we have two covariates: global mean and treatment (0 or 2). We have two time points and three individuals.  X.full might look like this.

   # | 1 0 |
   # | 1 2 |
   # | 1 0 |
   # | 1 2 |
   # | 1 0 |
   # | 1 2 |
   # An easy way to generate this matrix is to create the covariate matrix for one individual and then use Kronecker.
   # X.indiv = matrix(c(1,1,0,2),2,2) 
   # X.full = kronecker(matrix(1,n,1), X.indiv) # n=3
   
   # This does not easily generalize however and there are much more complicated situations.  It might best to draw
   # out what you want and then to figure out a simple way to construct it.  

  # THIS FUNCTION assumes that you are given a matrix X that is the shared covariate matrix for each time point.
  # That is each row represents an individual and the columns are the individuals covariates.
  # It then computes X.full using a kronecker product and returns.

  X.full = kronecker(X,matrix(1,m,1))
  return(X.full)
}


# ----------------- Analytical Power Section -------------------
#   
#  Here are a few functions that allow for the computation of analytical power in the mixed model.

# Computes the analytical power for a given setup.
# You must provide n,m, snp and the snpEffect size.
# If you have precomputed the covariance matrix and its inverse then provide it as Sigma.inv (see longGWAS.getSigma 
# Example: Sigma = longGWAS.getSigma(n,m,varComps,K,D); Sigma.inv = solve(Sigma)
# If not then you need to provide the varComps argument and the necessary matrices: K and D when m > 1
longGWAS.analyticalPower <- function(n, m, snp, snpEffect, K=NULL, D=NULL, varComps=NULL, Sigma.inv=NULL, X0=NULL, alpha=0.001){

  stopifnot(n > 0 && m > 0)

  if(var(snp) == 0){ return(NA) }

  N = n * m
  if(is.null(X0)){
    X = matrix(1,nrow=N,ncol=1)
  }
  else{ X = X0 }
  X = cbind(X, kronecker(snp, matrix(1,nrow=m,ncol=1)))

  if(is.null(Sigma.inv)){
    if(is.null(varComps) || is.null(K) || (is.null(D) && m > 1)){ 
       cat('ERROR: Must provide either Sigma.inv or variance components.\n')
       cat('When providing variance components you must also provide K and D (when m > 1).\n')
       return()
    }

    Sigma = longGWAS.getSigma(n,m,varComps,K,D)
    cat("Computing Sigma inverse...\n")
    Sigma.inv = solve(Sigma)
  }

  q = ncol(X) 
  ncp = longGWAS.powerNCP(snpEffect,X,Sigma.inv)

  v = qt(alpha/2,N-q,0)
  p = pt(v+sqrt(ncp),N-q) + 1 - pt(-v+sqrt(ncp),N-q)

  return(p)
}

longGWAS.getNCP <- function(delta,Sigma.inv,R,X){
  W = R %*% solve(t(X) %*% Sigma.inv %*% X) %*% t(R)
  return(t(delta) %*% solve(W) %*% delta)
}

longGWAS.powerNCP <- function(snpEffect,X,Sigma.inv){
  q = ncol(X)
  R = matrix(c(rep(0,q-1),1),nrow=1)
  ncp = longGWAS.getNCP(snpEffect,Sigma.inv,R,X); 
  return(ncp);
}

longGWAS.getSigma <- function(n,m,varComps,K,D=NULL){

    if(m == 1){
       Sigma = varComps$vg * K + diag(n) * (varComps$ve+varComps$vep)
    } 
    else{
       if(is.null(D)){cat('ERROR: Must provide D when m > 1\n'); return();}
       Sigma = varComps$vg * kronecker(K,matrix(1,m,m)) + varComps$ve * D + diag(n*m) * varComps$vep
    }
    return(Sigma)
}

# -------------- myFastLMM ----------------------------------
#   This code is hybrid between EMMA (Kang et al. 2008) and myFastLMM (Lippert et al. 2001) with my own spin.
#   This is a lightweight easy to understand version of code that fits a linear mixed model (LMM) with two variance components.
#   The basic idea is that when the variance of y has this structure V = g*K + e*I we can rotate the data by the eigenvectors
#   of K and then perform essentially basic linear regression where the residual i has variance g*lambda_i + e where lambda_i
#   is the ith eigenvalue for K.  All calculations are then based on this idea.
#
#   I have only made the REML version available here.
#   The optimization part () uses a heuristic to find the "best" solution - details in that function 

myFastLMM.cleanSNPs <- function(snps){

   nas = apply(snps,1,function(x){sum(is.na(x))})
   x = which( nas > 0 )
   cat("Number of SNPs with missing values: ", length(x),"\n")
   mafs = rowMeans(snps,na.rm=TRUE)
   if( length(x) > 0 ){
      for(i in x){ snps[i,which(is.na(snps[i,]))] <- mafs[i] } 
   }
   return(list(snps=snps,mafs=mafs))
}

# This calculation may take some time depending on the number of SNPs.  Therefore, it is advantageous
# to only do the calculation once and then save this result to a file.
# For example:
#  K = myFastLMM.kinship(snps)
#  write.table(K, 'output.K', row.names=F,col.names=F)
#  K.again = as.matrix(read.table('output.K')) # read it back in

myFastLMM.kinship <- function(snps){

   n = ncol(snps)
   K = emma.kinship(snps)
   S = diag(n) - 1/n*matrix(1,n,n)
   Ks = S %*% K %*% S
   diag(Ks) <- 1
   return(Ks)
}

myFastLMM.transform <- function(Y,X,n,q,K,K.eigen=NULL){

  if(is.null(K.eigen)){
    K.eigen = eigen(K)
  }

  H = K.eigen$vectors
  Ht = t(H)
  S = K.eigen$values

  Ys = Ht %*% Y
  Xs = Ht %*% X
  return(list(K.eigen=K.eigen, Ys=Ys, Xs=Xs, S=S, H=H))
}

myFastLMM.beta <- function(Ys,Xs,n,q,H,S,delta){
  #return(solve(t(Xs) %*% diag(1/(S+delta)) %*% Xs) %*% t(Xs) %*% diag(1/(S+delta)) %*% Ys)
  one = matrix(0,nrow=ncol(Xs),ncol=ncol(Xs))
  two = matrix(0,nrow=ncol(Xs),ncol=1) 
  for(i in 1:n){
    one = one + 1/(S[i]+delta)*matrix(Xs[i,],ncol=1) %*% Xs[i,]
    two = two + 1/(S[i]+delta) * matrix(Xs[i,],ncol=1) %*% Ys[i]
  }
  beta = solve(one) %*% two
  return(beta)
}

myFastLMM.sigma <- function(Ys,Xs,n,q,H,S,delta,beta=NULL){
  sigma = 0
  if(is.null(beta)){
    beta = myFastLMM.beta(Ys,Xs,n,q,H,S,delta)
  } 

  for(i in 1:n){
    sigma = sigma + 1/(S[i]+delta)*(Ys[i] - Xs[i,] %*% beta)**2
  }
  sigma = 1/n * sigma
  return(sigma)
}

myFastLMM.REML <- function(Ys,Xs,n,q,H,S,delta){

  R = myFastLMM.BetadBeta(Ys,Xs,n,q,H,S,delta)
  beta = R$beta
  sigma = myFastLMM.sigma(Ys,Xs,n,q,H,S,delta,beta)
  sigma = (sigma*n)/(n-q)

  one = n*log(2*pi*sigma) 
  three= q*log(2*pi*sigma) 

  two = 0
  for(i in 1:n){
    one = one + log(S[i]+delta) 
    two = two + 1/(S[i]+delta)*(Ys[i] - Xs[i,] %*% beta)**2
  }

  XX = H %*% Xs
  four = log(det(t(XX)%*%XX))
  LL = -1/2*(one + 1/sigma * two) + 1/2*(three + four - log(det(R$A)))
  return(list(delta=delta,beta=beta,LL=LL,sigma=sigma))
}

myFastLMM.dREML <- function(logdelta,Ys,Xs,n,q,H,S){
  delta = exp(logdelta)

  A = matrix(0,nrow=ncol(Xs),ncol=ncol(Xs))
  two = matrix(0,nrow=ncol(Xs),ncol=1) 
  dA = matrix(0,nrow=ncol(Xs),ncol=ncol(Xs))
  B = matrix(0,nrow=ncol(Xs),ncol=1) 
  dB = matrix(0,nrow=ncol(Xs),ncol=1) 
  Q = S + delta
  for(i in 1:n){
    PP = matrix(Xs[i,],ncol=1) %*% Xs[i,]
    MM = matrix(Xs[i,],ncol=1) %*% Ys[i]

    A = A + 1/Q[i] * PP
    dA = dA + (-1)/Q[i]**2 * PP
    B = B + 1/Q[i] * MM
    dB = dB + (-1)/Q[i]**2 * MM
  }
  A.inv = solve(A)
  beta = A.inv %*% B
  dBeta = -1*A.inv %*% dA %*% A.inv %*% B + A.inv %*% dB

  three = sum(diag(A.inv %*% dA)) # d(log |X'()X|) = Tr(X^{-1}dX)

  one = 0
  two = 0
  sigma = 0
  dSigma = 0

  for(i in 1:n){
    F = (Ys[i] - Xs[i,]%*%beta)**2
    one = one + 1/Q[i]
    #two = two - F/Q[i]**2
    two = two + (-2*(Ys[i] - Xs[i,]%*%beta) * (Xs[i,] %*% dBeta) * (S[i]+delta) - (Ys[i] - Xs[i,] %*% beta)**2)/(S[i]+delta)**2

    dSigma = dSigma + (-2*sqrt(F) * Q[i] * Xs[i,]%*%dBeta - F)/Q[i]**2
    sigma = sigma + F/Q[i]
  }
  sigma = 1/(n-q)*sigma
  dSigma = 1/(n-q)*dSigma

  dLL = -1/2*(one + n/(sigma*(n-q))*two)  + 1/2*(q/(sigma*(n-q)) * two - three)
  return(dLL)
}

myFastLMM.BetadBeta <- function(Ys,Xs,n,q,H,S,delta){
  A = matrix(0,nrow=ncol(Xs),ncol=ncol(Xs))
  B = matrix(0,nrow=ncol(Xs),ncol=1) 
  dA = matrix(0,nrow=ncol(Xs),ncol=ncol(Xs))
  dB = matrix(0,nrow=ncol(Xs),ncol=1) 
  for(i in 1:n){
    A = A + 1/(S[i]+delta) * matrix(Xs[i,],ncol=1) %*% Xs[i,]
    dA = dA + (-1)/(S[i]+delta)**2 * matrix(Xs[i,],ncol=1) %*% Xs[i,]
    B = B + 1/(S[i]+delta) * matrix(Xs[i,],ncol=1) %*% Ys[i]
    dB = dB + (-1)/(S[i]+delta)**2 * matrix(Xs[i,],ncol=1) %*% Ys[i]
  }

  #A.e = eigen(A)
  #cat("max = ",max(A.e$values)," min=",min(A.e$values)," ratio=",max(A.e$values)/min(A.e$values),"\n")

  if(length(A) > 1){ A.inv = solve(A) }
  else{ A.inv = 1/A }
  beta = A.inv %*% B
  dBeta = -1*A.inv %*% dA %*% A.inv %*% B + A.inv %*% dB
  return(list(beta=beta,dBeta=dBeta,stderrSQ=A.inv,A=A))
}

myFastLMM.optimizeDeltaREML <- function(Ys,Xs,n,q,H,S,K,Y,X, ngrids=100, llim=-10, ulim=10, esp=1e-10,verbose=FALSE){

  logdelta = (0:ngrids)/ngrids*(ulim-llim)+llim
  dLL = sapply(logdelta,function(x){myFastLMM.dREML(x,Ys,Xs,n,q,H,S)})

  optDeltas = c()
  N = length(dLL)

  # search for maximum using the derivative of the REML
  for(i in 1:(N-1)){
     if((dLL[i] * dLL[i+1] < 0 ) && (dLL[i] > 0) && (dLL[i+1] < 0)){
        optDeltas = append(optDeltas, uniroot(myFastLMM.dREML, lower=logdelta[i], upper=logdelta[i+1], Ys=Ys, Xs=Xs, n=n, q=q, H=H, S=S)$root) 
     } 
  }

  # Many times the function is monotonic so we add the boundary.
  optDeltas = append(optDeltas,llim)
  optDeltas = append(optDeltas,ulim)

  LLs = c()
  sigmas = c()
  if(verbose){
    cat("optDeltas:\n")
    cat(optDeltas,"\n")
  }

  # Calculate REML for each potential delta
  for(i in 1:length(optDeltas)){
    REMLE = myFastLMM.REML(Ys,Xs,n,q,H,S,exp(optDeltas[i]))
    LLs = append(LLs, REMLE$LL)
    sigmas = append(sigmas,REMLE$sigma)
  }

  x = which.max(LLs)
  maxLL = LLs[x]
  sigma = sigmas[x]
  ooptDelta = optDeltas[x]

  esp = (max(LLs) - min(LLs))/ngrids
  step = 0.1

  # If we are on the boundary this is likely not the best answer.
  # We use a heurstic approach to find a better solution.
  # We simply "backup" until we see that we have deviated significantly from the maxLL.
  # The rational is that on the boundary the function usually looks flat.  We want to identify
  # the position at which it leveled out.  This should be similar to what a hill-climbing algorithm would 
  # have found when the stopping condition was a small gradiant.

  # Therefore the next two conditions check the case when our best solution is on the boundary and 
  # attempt to find a better solution.  We have shown in simulation that this works much better
  # than the default in EMMA for example.
  if(optDeltas[x] == ulim){
  #if(FALSE){
    prevLL = maxLL 
    prevSigma = sigma
    for(i in seq(ulim-step,0,by=-step)){
      REMLE = myFastLMM.REML(Ys,Xs,n,q,H,S,exp(i))
      LL = REMLE$LL
      if(abs(LL - maxLL) > esp){
         ooptDelta = i + step; 
         maxLL = prevLL; 
         sigma = prevSigma;
         break
      }
      prevLL = LL
      prevSigma = REMLE$sigma
    } 
  }
  if(optDeltas[x] == llim){
  #if(FALSE){
    prevLL = maxLL
    for(i in seq(llim+step,0,by=step)){
      REMLE = myFastLMM.REML(Ys,Xs,n,q,H,S,exp(i))
      LL = REMLE$LL
      if(abs(LL - maxLL) > esp){
         ooptDelta = i - step; 
         maxLL = prevLL; 
         sigma = prevSigma;
         break
      }
      prevLL = LL
      prevSigma = REMLE$sigma
    } 
  }

  optDelta = exp(ooptDelta)
  beta = myFastLMM.beta(Ys,Xs,n,q,H,S,optDelta)

  if(verbose){
     cat("optDelta:\n")
     cat(ooptDelta,"\n")
     cat("max LL is ",maxLL,"\n")
  }

  return(list(LL=maxLL,delta=optDelta,vg=sigma,ve=optDelta*sigma,beta=beta))
}

# --------- Simple Utility Functions ------------------ # 
normMatrix <- function(M){
  return((M - matrix(rowMeans(M, na.rm=TRUE),nrow(M),ncol(M))) / matrix(apply(M,1,sd,na.rm=TRUE),nrow(M),ncol(M)));
}


#-------------- The following Functions were taken from the R version of EMMA. ----------# 
#-------------- EMMA can be found at http://mouse.cs.ucla.edu/emma/.
#-------------- I did not include the complete source because most functions are not needed.

emma.kinship <- function(snps, method="additive", use="all") {
  n0 <- sum(snps==0,na.rm=TRUE)
  nh <- sum(snps==0.5,na.rm=TRUE)
  n1 <- sum(snps==1,na.rm=TRUE)
  nNA <- sum(is.na(snps))

  stopifnot(n0+nh+n1+nNA == length(snps))

  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1

  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}
