LDA.abundance=function(y,ncomm,gamma,ngibbs,nburn){
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
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  llk=rep(NA,ngibbs)
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)   
    
    #sample z
    tmp=samplez(theta=theta, phi=phi, y=y, ncommun=ncomm, nloc=nloc, nspp=nspp)
    nlk=tmp$nlk
    # nlk=nlk.true
    nks=tmp$nks
    # nks=nks.true
    
    #get parameters  
    tmp=get.theta(nlk=nlk,gamma=gamma,ncomm=ncomm,nloc=nloc)
    vmat=tmp$vmat
    theta=tmp$theta
    # theta[theta>hi]=hi; theta[theta<lo]=lo
    # theta=theta.true
    
    phi=rdirichlet1(alpha=nks+1,ncomm=ncomm,nspp=nspp) 
    # phi[phi>hi]=hi; phi[phi<lo]=lo
    # phi=phi.true
    
    #calculate loglikelihood
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    
    #store results  
    llk[i]=sum(y*log(prob))
    theta.out[i,]=theta
    phi.out[i,]=phi
  }
  seq1=nburn:ngibbs
  list(llk=llk[seq1],
       theta=theta.out[seq1,],
       phi=phi.out[seq1,])
}
