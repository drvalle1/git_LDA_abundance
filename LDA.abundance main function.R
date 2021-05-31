LDA.abundance=function(y,ncomm,ngibbs,nburn,psi,gamma,mean.theta,mean.phi){
  #get data
  nspp=ncol(y)
  nloc=nrow(y)
  
  #useful stuff
  hi=0.999999
  lo=0.000001

  #initial values of parameters
  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(1/nspp,ncomm,nspp)

  #gibbs details
  if (!mean.theta) theta.out=matrix(NA,ngibbs,ncomm*nloc)
  if (mean.theta)  theta.out=matrix(0,nloc,ncomm)

  if (!mean.phi) phi.out=matrix(NA,ngibbs,ncomm*nspp)  
  if (mean.phi)  phi.out=matrix(0,ncomm,nspp)
  
  llk=rep(NA,ngibbs)
  # log.prior=rep(NA,ngibbs)
  
  options(warn=2)
  zeroes=array(0,dim=c(nloc,nspp,ncomm))
  for (i in 1:ngibbs){
    print(i)   
    
    #re-order z from time to time
    if (i<nburn & i%%50==0){
      med=apply(theta,2,median)
      ordem=order(med,decreasing=T)
      theta=theta[,ordem]
      phi=phi[ordem,]
    }
    
    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp, zeroes=zeroes)
    array.lsk=tmp$ArrayLSK
    nlk=tmp$nlk
    # nlk=nlk.true
    nks=tmp$nks
    # nks=nks.true
    
    #get parameters  
    tmp=get.theta(nlk=nlk,gamma=gamma,ncomm=ncomm,nloc=nloc,i=i,nburn=nburn,nks=nks,theta=theta,change1=20)
    nks=tmp$nks
    vmat=tmp$vmat
    theta=tmp$theta
    # theta[theta>hi]=hi; theta[theta<lo]=lo
    # theta=theta.true
    
    phi=rdirichlet1(alpha=nks+psi,ncomm=ncomm,nspp=nspp) 
    # phi[phi>hi]=hi; phi[phi<lo]=lo
    # phi=phi.true
    
    #calculate loglikelihood
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    
    #calculate log prior (I often get Inf!!)
    # vmat1=vmat[,-ncomm]
    # vmat1[vmat1>hi]=hi; vmat1[vmat1<lo]=lo
    # log.p.betas=sum(dbeta(vmat1,1,gamma,log=T))
    # log.p.phi=sum(log(ddirichlet(phi,rep(psi,nspp))))
    # print(c(log.p.betas,log.p.phi))
    
    #store results  
    llk[i]=sum(y*log(prob))
    # log.prior[i]=log.p.betas+log.p.phi
    
    if (mean.phi & i>=nburn)    phi.out=phi.out+phi
    if (!mean.phi)              phi.out[i,]=phi
    if (mean.theta & i>=nburn)  theta.out=theta.out+theta
    if (!mean.theta)            theta.out[i,]=theta
  }

  #output results
  seq1=nburn:ngibbs
  nseq=length(seq1)
  k=list(llk=llk,array.lsk=array.lsk)
  if (!mean.phi) k$phi=phi.out[seq1,]
  if (mean.phi)  k$phi=phi.out/nseq #just posterior mean
  if (!mean.theta) k$theta=theta.out[seq1,]
  if (mean.theta)  k$theta=theta.out/nseq #just posterior mean
  k
}
