rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

setwd('U:\\GIT_models\\git_LDA_abundance\\phi given theta')
source('estimate phi given theta function.R')

#get data
setwd('U:\\GIT_models\\git_LDA_abundance')
dat=read.csv('fake data5.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)
nspp=ncol(y)
ncomm=5

#get phi
setwd('U:\\GIT_models\\git_LDA_abundance')
theta=data.matrix(read.csv('theta true.csv',as.is=T))

#run folding operation
ngibbs=1000
nburn=ngibbs/2
psi=0.1
res=phi.given.theta(y=y,ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,theta=theta,psi=psi)

plot(res$llk,type='l')
seq1=1:nrow(res$phi)
phi=colMeans(res$phi[seq1,])
phi1=matrix(phi,ncol(theta),ncol(y))

phi.true=data.matrix(read.csv('phi true.csv',as.is=T))
plot(phi1,phi.true)
hist(phi1-phi.true)