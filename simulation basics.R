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