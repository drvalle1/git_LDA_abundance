rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

setwd('U:\\GIT_models\\git_LDA_abundance\\folding')
source('folding main function.R')

#get data
setwd('U:\\GIT_models\\git_LDA_abundance')
dat=read.csv('fake data5.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)
nspp=ncol(y)
ncomm=5

#get phi
setwd('U:\\GIT_models\\git_LDA_abundance')
phi=read.csv('phi true.csv',as.is=T)
npost=1000
phi.post=matrix(unlist(phi),npost,nspp*ncomm,byrow=T)

#run folding operation
ngibbs=1000
nburn=ngibbs/2
gamma=0.1
res=folding.in(y=y,ncomm=ncomm,ngibbs=ngibbs,nburn=nburn,phi.post=phi.post,gamma=gamma)

plot(res$llk,type='l')
seq1=1:nrow(res$theta)
theta=colMeans(res$theta[seq1,])
theta1=matrix(theta,nrow(y),ncomm)
boxplot(theta1)

plot(NA,xlim=c(0,nrow(y)),ylim=c(0,1))
for (i in 1:ncomm){
  lines(theta1[,i],col=i)
}