theta.given.phi=function(y,ncomm,ngibbs,nburn,phi.post,gamma){
  #get data
  nspp=ncol(y)
  nloc=nrow(y)
  
  #useful stuff
  hi=0.999999
  lo=0.000001

  #initial values of parameters
  theta=matrix(1/ncomm,nloc,ncomm)
  npost=nrow(phi.post)
  ind=sample(1:npost,size=1)
  ncomm.init=ncol(phi.post)/nspp
  tmp=matrix(phi.post[ind,],ncomm.init,nspp)
  phi=tmp[1:ncomm,]

  #gibbs details
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
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
    theta=rdirichlet2(gamma+nlk,nloc=nloc,ncomm=ncomm)
    # theta[theta>hi]=hi; theta[theta<lo]=lo
    # theta=theta.true

    ind=sample(1:npost,size=1)
    tmp=matrix(phi.post[ind,],ncomm.init,nspp)
    phi=tmp[1:ncomm,]
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
    theta.out[i,]=theta
  }
  seq1=nburn:ngibbs
  list(llk=llk[seq1],
       # log.prior=log.prior[seq1],
       theta=theta.out[seq1,])
}
