rm(list=c(ls()))
#setwd("/Users/piresmm/git/Lego/Lego_Program/")
library(igraph)
library(plotrix)
library(RColorBrewer)

source("R/build_template.R")
source("R/test_compatability.R")

#sequence <- seq(10,2000,100)
#prop_active <- numeric(length(sequence))
sequence <- seq(10,100,1)
num_e <- numeric(length(sequence))
tic <- 0
for (i in sequence) {
  
  #for (i in sequence) {
  
  tic <- tic + 1
  
  num_play <- i
  
  #   pw_prob <- c(
  #     pr_ne = 0.025,
  #     pr_nn = 0.025,
  #     pr_ni = 0.05,
  #     pr_nm = 0.005,
  #     pr_ia = 0.05,
  #     pr_ie = 0.2,
  #     pr_ii = 0.5,
  #     pr_aa = 0.05,
  #     pr_ee = 0.05
  #   )
  
  #Defining probabilities of each type
  #Basic probs could also be based on num_play. e.g., We should expect p.n*num_play n's per column/row
  
  p.n=0.02 
  p.e=0.1
  p.i=0.72
  p.m=0.1
  p.a=0
  
  #Normalization [0,1]
  S_prob=sum(c(p.n,p.e,p.i,p.m,p.a))
  p.n=p.n/S_prob
  p.e=p.e/S_prob
  p.i=p.i/S_prob
  p.m=p.m/S_prob
  p.a=p.a/S_prob
  
  #Defining paiwise probabilities 
  pw_prob <- c(
    pr_ne = p.n*(p.e/(p.e+p.n+p.i+p.m)),
    pr_nn = p.n*(p.n/(p.e+p.n+p.i+p.m)),
    pr_ni = p.n*(p.i/(p.e+p.n+p.i+p.m)),
    pr_nm = p.n*(p.m/(p.e+p.n+p.i+p.m)),
    pr_ia = p.i*(p.a/(p.e+p.a+p.n+p.i)),
    pr_ie = p.i*(p.e/(p.e+p.a+p.n+p.i)),
    pr_ii = p.i*(p.i/(p.e+p.a+p.n+p.i)),
    pr_aa = p.a*(p.a/(p.a+p.i)),
    pr_ee = p.e*(p.e/(p.i+p.n+p.e))
  )
  
  
  #make sure this vector sums to one
  pw_prob <- pw_prob / sum(pw_prob)
  
  #Build the interaction template
  int_m <- build_template(num_play,pw_prob, 0.8)
  
  num_e[tic] <- length(which(int_m == "e"))
}
plot(sequence,num_e,xlab="Template size",ylab="Number of trophic interactions",pch=16)

