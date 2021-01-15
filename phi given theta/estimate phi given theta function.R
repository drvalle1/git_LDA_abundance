phi.given.theta=function(y,ncomm,ngibbs,nburn,theta.post,psi){
  #get data
  nspp=ncol(y)
  nloc=nrow(y)
  
  #useful stuff
  hi=0.999999
  lo=0.000001

  #initial value
  phi=matrix(1/nspp,ncomm,nspp)
  npost=nrow(theta.post)
  ind=sample(1:npost,size=1)
  ncomm.init=ncol(theta.post)/nloc
  tmp=matrix(theta.post[ind,],nloc,ncomm.init)
  theta=tmp[,1:ncomm]
  
  #gibbs details
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  llk=rep(NA,ngibbs)
  # log.prior=rep(NA,ngibbs)
  
  options(warn=2)
  zeroes=array(0,dim=c(nloc,nspp,ncomm))
  for (i in 1:ngibbs){
    print(i)   
    
    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp, zeroes=zeroes)
    array.lsk=tmp$ArrayLSK
    nlk=tmp$nlk
    # nlk=nlk.true
    nks=tmp$nks
    # nks=nks.true
    
    #get parameters  
    phi=rdirichlet1(alpha=nks+psi,ncomm=ncomm,nspp=nspp) 

    #get theta
    ind=sample(1:npost,size=1)
    tmp=matrix(theta.post[ind,],nloc,ncomm.init)
    theta=tmp[,1:ncomm]
    
    #calculate loglikelihood
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    
    #store results  
    llk[i]=sum(y*log(prob))
    # log.prior[i]=log.p.betas+log.p.phi
    phi.out[i,]=phi
  }
  seq1=nburn:ngibbs
  list(llk=llk[seq1],
       # log.prior=log.prior[seq1],
       phi=phi.out[seq1,])
}
