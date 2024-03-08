#include packages
library(tidyverse)
library(rjags)
library(mcmcplots)
library(MCMCvis)

#load data 
load("DataForAnalysis.RData")

#data
dat <- list("y" = y, "V" = V, "x" = x, "z" = z, "w" = w, "J"=J, "S" = S,
            "yind" = yind, "numFU" = numFU, "numcont" = numcont, "ref_arm" = ref_arm,
            "yend" = yend, "dim" = dim, "tind" = tind,
            "N" = N, "n" = n, "p" = p, "q" = q, "l"=l)


jags.m <- jags.model(file = "model_primary_RE.txt", data=dat, n.chains = 4, n.adapt = 10000) 
#jags.m <- jags.model(file = "model_secondary_FE.txt", data=dat, n.chains = 4, n.adapt = 10000) 
#jags.m <- jags.model(file = "model_SSVS.txt", data=dat, n.chains = 4, n.adapt = 10000) 

params <- c("alpha", "beta", "gamma", "phi", "eta", "tau") 

## run JAGS and save posterior samples
samps <- coda.samples( jags.m, params, n.iter=70000 )
## get DIC samples to check model fit
dicsamps <- dic.samples(jags.m, n.iter=20000) 
DIC <- sum(dicsamps$deviance) + sum(dicsamps$penalty) 

save.image(file = "Results.RData")
