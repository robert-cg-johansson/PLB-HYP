library(MASS)
library(tidyverse)





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
  
  data = as.data.frame(c(hyp_vs_baseline, plb_vs_baseline))
  group = c(rep("hyp", n), rep("plb", n))
  data = cbind(data, group)
  names(data) = c("tolerance", "group")

  ci = t.test(tolerance ~ group, data = data, paired = T)$conf.int
  decision = if(SESOI < ci[1]){"H1"} else if(SESOI > ci[2]){"H0"} else {"inconclusive"}
  

  
  return(decision)
}





iter = 1000



#### scenario 1: effect size = g = 0


results = replicate(iter, simul(  n = 100,
                                  mean_hyp = 100,
                                  mean_plb = 100,
                                  mean_baseline = 70,
                                  sd = 63,
                                  r = 0.7,
                                  SESOI = 10))

H1_detection_rate = mean(results=="H1")
H1_detection_rate

#### scenario 2: effect size = g = 0.15


results = replicate(iter, simul(  n = 200,
                                  mean_hyp = 120,
                                  mean_plb = 100,
                                  mean_baseline = 70,
                                  sd = 63,
                                  r = 0.7,
                                  SESOI = 10))

H1_detection_rate = mean(results=="H1")

H1_detection_rate
