library(MASS)
library(tidyverse)
library(BayesFactor)
library(HDInterval)
library(pkgcond)



simul = function(
  n,
  mean_hyp,
  mean_plb,
  mean_baseline,
  sd,
  r,
  SESOI
){
  tolerance = mvrnorm(n = n,
                      mu = c(0, 0, 0),
                      Sigma = matrix(c(1, r, r,
                                       r, 1, r,
                                       r, r, 1), nrow = 3))
  
  baseline = tolerance[,1]*sd+mean_baseline
  hyp = tolerance[,2]*sd+mean_hyp
  plb = tolerance[,3]*sd+mean_plb
  
  hyp_vs_baseline = hyp-baseline
  plb_vs_baseline = plb-baseline
  diff_hyp_plb = hyp_vs_baseline - plb_vs_baseline
  
  data = as.data.frame(c(hyp_vs_baseline, plb_vs_baseline))
  group = c(rep("hyp", n), rep("plb", n))
  data = cbind(data, group)
  names(data) = c("tolerance", "group")
  
  bf = ttestBF(x = diff_hyp_plb)
  posterior = posterior(bf, iterations = 1000)
  HDI = hdi(posterior, ci = 0.9)[,"mu"]
  
  decision = if(SESOI < HDI["lower"]){"H1"} else if(SESOI > HDI["upper"]){"H0"} else {"inconclusive"}
  

  
  return(decision)
}





iter = 1000



#### scenario 1: effect size = g = 0


results = replicate(iter, simul(  n = 120,
                                  mean_hyp = 100,
                                  mean_plb = 100,
                                  mean_baseline = 70,
                                  sd = 63,
                                  r = 0.7,
                                  SESOI = 15))

H1_detection_rate = mean(results=="H1")
H1_detection_rate

H0_detection_rate = mean(results=="H0")
H0_detection_rate

#### scenario 2: effect size = ...


results = replicate(iter, simul(  n = 120,
                                  mean_hyp = 115,
                                  mean_plb = 90,
                                  mean_baseline = 70,
                                  sd = 63,
                                  r = 0.7,
                                  SESOI = 10))

H1_detection_rate = mean(results=="H1")
H1_detection_rate

H0_detection_rate = mean(results=="H0")
H0_detection_rate



##################################################
### Same code as above, but with Bayesian MCMC ###
##################################################

### Set up your workspace

rm(list=ls())
bugsdir <- "C:/WinBUGS14" # Download WinBUGS from https://www.mrc-bsu.cam.ac.uk/software/bugs/
setwd("C:/ExcercisesforR/")

library(MASS)
library(tidyverse)
library(BayesFactor)
library(HDInterval)
library(pkgcond)
library(R2WinBUGS)

### Set up a simulation function with MCMC

simul = function(
  n,
  mean_hyp,
  mean_plb,
  mean_baseline,
  sd,
  r,
  SESOI
){
  tolerance = mvrnorm(n = n,
                      mu = c(0, 0, 0),
                      Sigma = matrix(c(1, r, r,
                                       r, 1, r,
                                       r, r, 1), nrow = 3))
  
  baseline = tolerance[,1]*sd+mean_baseline
  hyp = tolerance[,2]*sd+mean_hyp
  plb = tolerance[,3]*sd+mean_plb
  
  hyp_vs_baseline = hyp-baseline
  plb_vs_baseline = plb-baseline
  diff_hyp_plb = hyp_vs_baseline - plb_vs_baseline
  
  data = as.data.frame(c(hyp_vs_baseline, plb_vs_baseline))
  group = c(rep("hyp", n), rep("plb", n))
  data = cbind(data, group)
  names(data) = c("tolerance", "group")
  
  #Feed Data to winbugs
  x <- diff_hyp_plb
  x <- x/sd(x)

  #Set number of subjects
  ndata <- length(diff_hyp_plb)
  
  # Pass data to WinBUGS
  data  <- list("x", "ndata") # to be passed on to WinBUGS
  
  myinits <- list(
    list(deltatmp = runif(1), sigmatmp = runif(1)),
    list(deltatmp = runif(1), sigmatmp = runif(1)),
    list(deltatmp = runif(1), sigmatmp = runif(1)))
  
  # Parameters to be monitored
  parameters <- c("delta")
  
  # The following command calls WinBUGS with specific options.
  samples <- bugs(data, inits=myinits, parameters,
                  model.file="PLB-HYP-means.txt",
                  n.chains=3, n.iter=10000, n.burnin=1000, n.thin=1,
                  DIC=T, bugs.directory=bugsdir,
                  codaPkg=F, debug=FALSE)
  
  #Extract HDI of posterior
  HDI = hdi(samples$sims.list$delta, ci = 0.9)
  HDI["lower"]
  HDI["upper"]
  
  decision = if(SESOI < HDI["lower"]){"H1"} else if(SESOI > HDI["upper"]){"H0"} else {"inconclusive"}
  
  
  
  return(decision)
}


#Set n of iterations
iter = 10


#### scenario 1: effect size = g = 0


results = replicate(iter, simul(  n = 120,
                                  mean_hyp = 100,
                                  mean_plb = 100,
                                  mean_baseline = 70,
                                  sd = 63,
                                  r = 0.7,
                                  SESOI = 0.2))

H1_detection_rate = mean(results=="H1")
H1_detection_rate


H0_detection_rate = mean(results=="H0")
H0_detection_rate

#### scenario 2: effect size = ...


results = replicate(iter, simul(  n = 120,
                                  mean_hyp = 120,
                                  mean_plb = 100,
                                  mean_baseline = 70,
                                  sd = 63,
                                  r = 0.7,
                                  SESOI = 0.2))

H1_detection_rate = mean(results=="H1")
H1_detection_rate

H0_detection_rate = mean(results=="H0")
H0_detection_rate
