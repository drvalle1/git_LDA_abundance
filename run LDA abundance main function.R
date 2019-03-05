rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
setwd('U:\\GIT_models\\git_LDA_abundance')
source('gibbs functions.R')
source('LDA.abundance main function.R')
sourceCpp('aux1.cpp')

dat=read.csv('fake data5.csv',as.is=T)
ind=which(colnames(dat)=='X')
y=data.matrix(dat[,-ind]); dim(y)

ncomm=30
gamma=0.1
ngibbs=1000
nburn=ngibbs/2
res=LDA.abundance(y=y,ncomm=ncomm,gamma=gamma,
                  ngibbs=ngibbs,nburn=nburn)

plot(res$llk,type='l')
theta=res$theta[nrow(res$theta),]
theta1=matrix(theta,nrow(y),ncomm)
boxplot(theta1)

plot(NA,xlim=c(0,nrow(y)),ylim=c(0,1))
for (i in 1:ncomm){
  lines(theta1[,i],col=i)
}