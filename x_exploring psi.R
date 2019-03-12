library(MCMCpack)
rdirichlet

set.seed(1)
n=10000
x=rgamma(n,0.001,1)
mean(x==0)

x=rgamma(n,0.01,1)
mean(x==0)