rm(list=ls())
library(devtools)
# install_github("variani/lme4qtl")
library(lme4qtl)

# load synthetic data set `dat40` distributed within `lme4qtl`
# - table of phenotypes `dat40`
# - the double kinship matrix `kin2`
data(dat40)
dat40
str(kin2)

# (1) model continiuous trait `trait1`
mod <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), 
									dat40,
									relmat = list(ID = kin2))

# get the estimation of h2
(vf <- as.data.frame(VarCorr(mod))[, c("grp", "vcov")])
#       grp      vcov
#1       ID 5.2845001
#2    FAMID 0.0000000
#3 Residual 0.6172059

prop <- with(vf, vcov / sum(vcov))

(h2 <- prop[1]) 
#[1] `0.895419`

# (2) model binary trait `trait1bin`
gmod <- relmatGlmer(trait1bin ~ (1|ID), 
										dat40, 
										relmat = list(ID = kin2), 
										family = binomial)
gmod


