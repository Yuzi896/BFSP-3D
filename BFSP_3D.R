library(fields);library(extraDistr);library(Matrix);library(sparseMVN);library(FastGP)

d2c<-function(x,y,z,center){
  d2c<-sqrt((x-center[1])^2+(y-center[2])^2+(z-center[3])^2)
  return(d2c)
}

a2c<-function(x,y,z,center){
  theta<-atan2(sqrt((x-center[1])^2+(y-center[2])^2),(z-center[3]))
  phi<-atan2((y-center[2]),(x-center[1]))
  phi<-ifelse(phi<0,phi+2*pi,phi)
  
  return(cbind(phi,theta))
}


calc.lik<-function(n,loc,value,params,cluster,spatial){
  
  if(n==0){
    ind0<-which(cluster==0)
    x<-loc[,1]
    y<-loc[,2]
    z<-loc[,3]
    x0<-x[ind0]
    y0<-y[ind0]
    z0<-z[ind0]
    loc0<-loc[ind0,]
    center0<-c(mean(x0),mean(y0),mean(z0))
    
    index<-list()
    index[[1]]<-ind0[which(x0<center0[1]&y0<center0[2]&z0<center0[3])]
    index[[2]]<-ind0[which(x0<center0[1]&y0>=center0[2]&z0<center0[3])]
    index[[3]]<-ind0[which(x0>=center0[1]&y0<center0[2]&z0<center0[3])]
    index[[4]]<-ind0[which(x0>=center0[1]&y0>=center0[2]&z0<center0[3])]
    index[[5]]<-ind0[which(x0<center0[1]&y0<center0[2]&z0>=center0[3])]
    index[[6]]<-ind0[which(x0<center0[1]&y0>=center0[2]&z0>=center0[3])]
    index[[7]]<-ind0[which(x0>=center0[1]&y0<center0[2]&z0>=center0[3])]
    index[[8]]<-ind0[which(x0>=center0[1]&y0>=center0[2]&z0>=center0[3])]
    
    
    lik8<-sapply(1:8,ind.lik,value,params,index,loc,spatial)
    log.lik<-sum(lik8)
    
  }else{
    ind1<-which(cluster==1)
    Y1<-value[ind1]
    n1<-length(Y1)
    if(n1<2) return(-Inf)
    
    if(spatial==TRUE){
      # Spatial
      cov<-params[6]*stationary.taper.cov(loc[ind1,],theta=1/params[4],Taper.args=list(k=2,aRange=.15,dimension=3))+params[2]*diag.spam(n1)
    }else{cov<-params[2]*diag.spam(n1)}

    cov<-params[2]*diag.spam(n1)
    log.lik<- dmvn.sparse(t(Y1), rep(params[5],n1), Matrix::Cholesky(as.dgCMatrix.spam(cov)), prec = F, log = TRUE)
  }
  return(log.lik)
}

ind.lik=function(indd,value,params,index,loc,spatial){
  ind<-index[[indd]]
  Y0<-value[ind]
  n0<-length(Y0)
  if(spatial==TRUE){
    cov<-params[3]*stationary.taper.cov(loc[ind,],theta=1/params[4],Taper.args=list(k=2,aRange=.15,dimension=3))+params[2]*diag.spam(n0)
  }else{cov<-params[2]*diag.spam(n0)}

  cov<-params[2]*diag.spam(n0)
  if(n0>0){log.lik<-dmvn.sparse(t(Y0), rep(params[1],n0), Cholesky(as.dgCMatrix.spam(cov)), prec = F, log = TRUE)}else{log.lik<-0}
  log.lik
}

likelihood=function(X_phi,Y_theta,intercept,R,R_minus,value,x,y,z,center,params,angles,M,spatial){
  
  f=(X_phi*Y_theta)%*%c(R,R_minus)+t((intercept)%*%(cos(outer(0:M,angles[,2]))))
  
  
  fx=f*sin(angles[,2])*cos(angles[,1])+center[1]
  fy=f*sin(angles[,2])*sin(angles[,1])+center[2]
  fz=f*cos(angles[,2])+center[3]
  
  cluster<-+(d2c(x,y,z,center)<=f)
  
  
  loc=cbind(x,y,z)
  
  log.lik<-sapply(0:1,calc.lik,loc=loc,value=value,params=params,cluster=cluster,spatial=spatial)
  return(list(sum(log.lik),cluster,f))
}

prior_cluster=function(X_phi,Y_theta,intercept,R,R_minus,angles,center,x,y,z,M){
  
  f=(X_phi*Y_theta)%*%c(R,R_minus)+t((intercept)%*%(cos(outer(0:M,angles[,2]))))
  
  fx=f*sin(angles[,2])*cos(angles[,1])+center[1]
  fy=f*sin(angles[,2])*sin(angles[,1])+center[2]
  fz=f*cos(angles[,2])+center[3]
  cluster<-+(d2c(x,y,z,center)<=f)
  
  # make sure the boundary does not exceed 
  if(any(fx>1.1)|any(fx<(-.1))|any(fy>1.1)|any(fy<(-.1))|any(fz>1.1)|any(fz<(-.1))){return(-Inf)}
  if(sum(cluster)<1){return(-Inf)}
  
  return(0)
}


propose_new<-function(param,cov){
  prop<-rcpp_rmvnorm(1,cov,param)
  return(prop)
}



prior_param=function(params,spatial){

  if(spatial==TRUE){
    # Let's try change priors to IG(2,1) instead of IG(0.1,0.1) for the variance 
    prior=dnorm(params[1],0,10,log=T)+sum(dgamma(1/params[2:3],2,1,log=T))+dgamma(params[4],3,.5,log=T)+dhnorm(params[5]-params[1], sigma = 10, log = TRUE)+dgamma(1/params[6],2,1,log=T)
  }else{prior=dnorm(params[1],0,10,log=T)+sum(dgamma(1/params[2],2,1,log=T))+dhnorm(params[5]-params[1], sigma = 10, log = TRUE)} 
   return(prior)
}

prior_intercept=function(intercept_params){
  prior=sum(dnorm(intercept_params[-1],0,1),log=T)+dunif(intercept_params[1],0,.5,log=T)
  return(prior)
}


update_R=function(X_phi,Y_theta,intercept,R,R_minus,value,x,y,z,center,i,ll,la0,lai,gam0,gami,laR,gamR,laR_minus,gamR_minus,params,angles,M,spatial){
  
  if(spatial==TRUE){
    if(i<=100){param_star<-propose_new(c(params),(.01/6)*diag(6))
    }else{
      param_star<-.99*propose_new(c(params),exp(la0)*gam0)+.01*propose_new(c(params),(.01/6)*diag(6))
      }
  }else{
    if(i <=100){
      param_star <- propose_new(c(params[c(1,2,5)]),(0.01/3)*diag(3))
      param_star <- c(param_star[1:2],0,0,param_star[3],0)
    }else{
      param_star<-.99*propose_new(c(params[c(1,2,5)]),exp(la0)*gam0)+.01* propose_new(c(params[c(1,2,5)]),(0.01/3)*diag(3))
      param_star <- c(param_star[1:2],0,0,param_star[3],0)
    }
  }

  prior_ratio=prior_param(param_star,spatial)-prior_param(params,spatial)
  
  if(prior_ratio>-Inf){
    ll_star<-likelihood(X_phi,Y_theta,intercept,R,R_minus,value,x,y,z,center,param_star,angles,M,spatial)
    like_ratio=ll_star[[1]]-ll[[1]]
    alpha=exp(prior_ratio+like_ratio)
    if(is.na(alpha)){alpha=0}
    if(runif(1)<alpha){
      params<-param_star
      ll<-ll_star
    }
  }
  
  if(i<=100) {
    param_star<-propose_new(c(intercept),(.01/(M+1))*diag((M+1)))
  }else{
    param_star<-.99*propose_new(c(intercept),exp(lai)*gami)+.01*propose_new(c(intercept),(.01/(M+1))*diag((M+1)))}
  
  
  prior_ratio=prior_intercept(param_star)-prior_intercept(intercept)+prior_cluster(X_phi,Y_theta,param_star,R,R_minus,angles,center,x,y,z,M)
  
  if(prior_ratio>-Inf){
    ll_star<-likelihood(X_phi,Y_theta,param_star,R,R_minus,value,x,y,z,center,params,angles,M,spatial)
    like_ratio=ll_star[[1]]-ll[[1]]
    alpha=exp(prior_ratio+like_ratio)
    if(is.na(alpha)){alpha=0}
    if(runif(1)<alpha){
      intercept<-param_star
      ll<-ll_star
    }
  }
  index<-rep(1:M,each=M)
  for(m in 1:M){
    if(i<=100) {
      param_star<-propose_new(c(R[index==m]),(.01/M)*diag(M))
    }else{
      param_star<-.99*propose_new(c(R[index==m]),exp(laR[[m]])*gamR[[m]])+.01*propose_new(c(R[index==m]),(.01/M)*diag(M))}
    
    R_star<-R
    R_star[index==m]<-param_star
    prior_ratio=sum(dnorm(param_star,0,1,log=T))-sum(dnorm(R[index==m],0,1,log=T))+prior_cluster(X_phi,Y_theta,intercept,R_star,R_minus,angles,center,x,y,z,M)
    
    if(prior_ratio>-Inf){
      ll_star<-likelihood(X_phi,Y_theta,intercept,R_star,R_minus,value,x,y,z,center,params,angles,M,spatial)
      like_ratio=ll_star[[1]]-ll[[1]]
      alpha=exp(prior_ratio+like_ratio)
      if(is.na(alpha)){alpha=0}
      if(runif(1)<alpha){
        R<-R_star
        ll<-ll_star
      }
    }
    
    if(i<=100) {
      param_star<-propose_new(c(R_minus[index==m]),(.01/M)*diag(M))
    }else{
      param_star<-.99*propose_new(c(R_minus[index==m]),exp(laR_minus[[m]])*gamR_minus[[m]])+.01*propose_new(c(R_minus[index==m]),(.01/M)*diag(M))}
    
    R_minus_star<-R_minus
    R_minus_star[index==m]<-param_star
    prior_ratio=sum(dnorm(param_star,0,1,log=T))-sum(dnorm(R_minus[index==m],0,1,log=T))+prior_cluster(X_phi,Y_theta,intercept,R,R_minus_star,angles,center,x,y,z,M)
    
    if(prior_ratio>-Inf){
      ll_star<-likelihood(X_phi,Y_theta,intercept,R,R_minus_star,value,x,y,z,center,params,angles,M,spatial)
      like_ratio=ll_star[[1]]-ll[[1]]
      alpha=exp(prior_ratio+like_ratio)
      if(is.na(alpha)){alpha=0}
      if(runif(1)<alpha){
        R_minus<-R_minus_star
        ll<-ll_star
      }
    }
  }
  
  
  return(list(R,R_minus,intercept,params,ll))
}

try.solve <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    if (!silent) {return(NA)}
    else{code}})}




range01 <- function(x){(x-min(x))/(max(x)-min(x))}



bfsp3d.partition<-function(data,iterations=10000,M=5,burn=1,vary_centroid=TRUE,center,spatial=FALSE,initial=NULL){
  
  
  x<-data[[1]]
  y<-data[[2]]
  z<-data[[3]]
  value<-data[[4]]
  
  if(length(initial)>0){
    intercept = initial[[1]]
    R = initial[[2]]
    R_minus = initial[[3]]
    params = initial[[4]]
  }else{
    # We want no spatial component so lambda = 0, sigma1^2 = 0, and sigma0^2 = 0
    if(spatial==TRUE){
      params<-c(0,var(value)/2,var(value)/2,.1,1,var(value)/2)
    }else{params<-c(0,var(value)/2,var(value)/2,0,1,0) }
  
    intercept<-c(.2,rep(0,M)) 
    R=rep(0,M^2) 
    R_minus=rep(0,M^2) 
  }

  angles<-a2c(x,y,z,center)
  thetas<-angles[,2]
  phis<-angles[,1]
  
  evens<-which(ceiling(1:(2*M*M)/M)%%2==0)
  odds<-c(1:(2*M*M))[-evens]
  Y<-array(dim=c(length(x),2*M*M)) 
  Y[,evens]<-t(sin(thetas)*sin(outer(1:M,thetas)))
  Y[,odds]<-t(sin(outer(1:M,thetas)))
  Y_theta<-Y
  
  X<-cbind(t(cos(outer(1:M,phis))),t(sin(outer(1:M,phis)))) 
  X_phi<-cbind(X[, rep(1:(2*M), each=M)]) 
  f=(X_phi*Y_theta)%*%c(R,R_minus)+t((intercept)%*%(cos(outer(0:M,thetas))))
 
  
 
  
  cluster<-+(d2c(x,y,z,center)<=f)
  ll<-likelihood(X_phi,Y_theta,intercept,R,R_minus,value,x,y,z,center,params,angles,M,spatial)
  
  
  R_keeps<-array(dim=c(iterations,M^2)) 
  R_minus_keeps<-array(dim=c(iterations,M^2))
  intercept_keeps<-array(dim=c(iterations,M+1))
  center_keeps<-array(dim=c(iterations,3))
  param_keeps<-array(dim=c(iterations,6)) 
  ll_keeps<-array(dim=c(iterations,1))
  
  index<-rep(1:M,each=M)
  muR<-list();muR_minus<-list();gamR<-list();gamR_minus<-list();laR<-list();laR_minus<-list()
  for(m in 1:M){
    muR[[m]]<-c(R[index==m])
    muR_minus[[m]]<-c(R_minus[index==m])
    gamR[[m]]<-diag(M)
    gamR_minus[[m]]<-diag(M)
    laR[[m]]<--1
    laR_minus[[m]]<--1
  }
  
  la0<--1
  lai<--1
  
  if(spatial==TRUE){
    mu0<-c(params)
    gam0<-diag(6)
  }else{
    mu0<-c(params[c(1,2,5)])
    gam0 <- diag(3)}
  
  mui<-c(intercept)
  gami<-diag(M+1)
  
  start.time <- Sys.time()
  for(i in 1:iterations){
    
    
    update=update_R(X_phi,Y_theta,intercept,R,R_minus,value,x,y,z,center,i,ll,la0,lai,gam0,gami,laR,gamR,laR_minus,gamR_minus,params,angles,M,spatial)
    R<-update[[1]]
    R_minus<-update[[2]]
    intercept<-update[[3]]
    params<-update[[4]]
    ll<-update[[5]]
    cluster<-ll[[2]]
    f<-ll[[3]]
    
    if(i<10){
      cat(paste(i,"\n","\n"))
    }
    
    if(i%%100==0){
      cat(paste(i,"\n","\n"))
    }
    R_keeps[i,]<-c(R)
    R_minus_keeps[i,]<-c(R_minus)
    intercept_keeps[i,]<-intercept
    param_keeps[i,]<-params
    center_keeps[i,]<-center
    ll_keeps[i,]<-ll[[1]]
    
    # variance adaption step for adaptive metropolis algorithm adjust term = 1/iteration+1
    index<-rep(1:M,each=M)
    if(i>=50){
      if(spatial==TRUE){
        mu0<-mu0+(1/(i+1))*(c(params)-mu0)
        gam0<-gam0+(1/(i+1))*((c(params)-mu0)%*%t(c(params)-mu0)-gam0)
      }else{
        mu0<-mu0+(1/(i+1))*(c(params[c(1,2,5)])-mu0)
        gam0<-gam0+(1/(i+1))*((c(params[c(1,2,5)])-mu0)%*%t(c(params[c(1,2,5)])-mu0)-gam0)
      }

      mui<-mui+(1/(i+1))*(c(intercept)-mui)
      gami<-gami+(1/(i+1))*((c(intercept)-mui)%*%t(c(intercept)-mui)-gami)
      
      for(m in 1:M){
        muR[[m]]<-muR[[m]]+(1/(i+1))*(c(R[index==m])-muR[[m]])
        gamR[[m]]<-gamR[[m]]+(1/(i+1))*((c(R[index==m])-muR[[m]])%*%t(c(R[index==m])-muR[[m]])-gamR[[m]])

        muR_minus[[m]]<-muR_minus[[m]]+(1/(i+1))*(c(R_minus[index==m])-muR_minus[[m]])
        gamR_minus[[m]]<-gamR_minus[[m]]+(1/(i+1))*((c(R_minus[index==m])-muR_minus[[m]])%*%t(c(R_minus[index==m])-muR_minus[[m]])-gamR_minus[[m]])
      }
      
      
    }
    
    # control step size
    if(i%%100==0 & i>=100){
      last20<-param_keeps[(i-20):i]
      accept<-1-mean(duplicated(last20))
      if(accept<.4) la0<-la0-.1
      if(accept>.4) la0<-la0+.1
      
      last20<-intercept_keeps[(i-20):i]
      accept<-1-mean(duplicated(last20))
      if(accept<.4) lai<-lai-.1
      if(accept>.4) lai<-lai+.1
      
      for(m in 1:M){
        last20<-R_keeps[(i-20):i,3*m]
        accept<-1-mean(duplicated(last20))
        if(accept<.4) laR[[m]]<-laR[[m]]-.1
        if(accept>.4) laR[[m]]<-laR[[m]]+.1

        last20<-R_minus_keeps[(i-20):i,3*m]
        accept<-1-mean(duplicated(last20))
        if(accept<.4) laR_minus[[m]]<-laR_minus[[m]]-.1
        if(accept>.4) laR_minus[[m]]<-laR_minus[[m]]+.1
      }
      
    }
    
    #update centroid 
    if(vary_centroid==T){center_star=c(mean(x[cluster==1]),mean(y[cluster==1]),mean(z[cluster==1]))}	
    if(vary_centroid==F){center_star<-center}
    if(sum(cluster)>1&(abs(center_star[1]-center[1])>.01|abs(center_star[2]-center[2])>.01|abs(center_star[3]-center[3])>.01)){
      
      bp<-cbind(fx=f*sin(angles[,2])*cos(angles[,1])+center[1],fy=f*sin(angles[,2])*sin(angles[,1])+center[2],fz=f*cos(angles[,2])+center[3])
      angle_star<-a2c(bp[,1],bp[,2],bp[,3],center_star)
      f_star<-d2c(bp[,1],bp[,2],bp[,3],center_star)
      thetas_star<-angle_star[,2]
      phis_star<-angle_star[,1]
      Y_star<-array(dim=c(length(x),2*M*M))
      Y_star[,evens]<-t(sin(thetas_star)*sin(outer(1:M,thetas_star)))
      Y_star[,odds]<-t(sin(outer(1:M,thetas_star)))
      Y_theta_star<-Y_star
      
      X_star<-cbind(t(cos(outer(1:M,phis_star))),t(sin(outer(1:M,phis_star))))
      X_phi_star<-cbind(X_star[, rep(1:(2*M), each=M)])
      X_star<-cbind(t(cos(outer(0:M,thetas_star))),X_phi_star*Y_theta_star)
      beta_star<-try.solve(solve(t(X_star)%*%X_star)%*%t(X_star)%*%f_star)
      if(beta_star[1]<0){
        beta_star[1]<-0
      }
      angle_star<-a2c(x,y,z,center_star)
      thetas_star<-angle_star[,2]
      phis_star<-angle_star[,1]
      Y_star<-array(dim=c(length(x),2*M*M))
      Y_star[,evens]<-t(sin(thetas_star)*sin(outer(1:M,thetas_star)))
      Y_star[,odds]<-t(sin(outer(1:M,thetas_star)))
      Y_theta_star<-Y_star
      
      X_star<-cbind(t(cos(outer(1:M,phis_star))),t(sin(outer(1:M,phis_star))))
      X_phi_star<-cbind(X_star[, rep(1:(2*M), each=M)])
      X_star<-cbind(t(cos(outer(0:M,thetas_star))),X_phi_star*Y_theta_star)
      prior_center<-prior_cluster(X_phi_star,Y_theta_star,beta_star[1:(M+1)],beta_star[(M+2):(M+1+M^2)],beta_star[(M+2+M^2):(M+M^2+1+M^2)],angle_star,center_star,x,y,z,M)
      if(prior_center>-Inf){
        intercept<-beta_star[1:(M+1)]
        angles<-angle_star
        R<-beta_star[(M+2):(M+1+M^2)]
        R_minus<-beta_star[(M+2+M^2):(M+M^2+1+M^2)]
        X_phi<-X_phi_star
        Y_theta<-Y_theta_star
        center<-center_star
        ll<-likelihood(X_phi,Y_theta,intercept,R,R_minus,value,x,y,z,center,params,angles,M,spatial)
      }
      
    }
    
    
  }
  end.time <- Sys.time()
  time.taken<-end.time-start.time
  
  print(laR)
  print(laR_minus)
  print(la0)
  print(lai)
  
  list1<-list(x,y,z,value,R_keeps,R_minus_keeps,param_keeps,intercept_keeps,center_keeps,time.taken,ll_keeps)
  list1
}

spher2cart <- function(r, theta, phi){
  
  x <- r * sin(theta) * cos(phi)
  y <- r * sin(theta) * sin(phi)
  z <- r * cos(theta)
  
  return(list(x = x, y = y, z = z))
}

update_fx <- function(fhat_star, theta_star,phi_star,evens, odds, M, x,y,z,center){
  Y_star <- array(dim=c(length(theta_star),2*M*M))
  Y_star[,evens]<-t(sin(theta_star)*sin(outer(1:M,theta_star)))
  Y_star[,odds]<-t(sin(outer(1:M,theta_star)))
  Y_theta_star<-Y_star
  
  X_star<-cbind(t(cos(outer(1:M,phi_star))),t(sin(outer(1:M,phi_star))))
  X_phi_star<-cbind(X_star[,rep(1:(2*M),each=M)])
  X_star<-cbind(t(cos(outer(0:M,theta_star))),X_phi_star*Y_theta_star)
  beta_star<-solve(t(X_star)%*%X_star)%*%t(X_star)%*%fhat_star
  
  # match the format 
  intercept<-beta_star[1:(M+1)]
  R<-beta_star[(M+2):(M+1+M^2)]
  R_minus<-beta_star[(M+2+M^2):(M+M^2+1+M^2)]
  angles<-a2c(x,y,z,center)
  thetas<-angles[,2]
  phis<-angles[,1]
  Y<-array(dim=c(length(thetas),2*M*M))
  Y[,evens]<-t(sin(thetas)*sin(outer(1:M,thetas)))
  Y[,odds]<-t(sin(outer(1:M,thetas)))
  Y_theta<-Y
  X<-cbind(1,t(cos(outer(rep(1:M,each=M),phis))),t(sin(outer(rep(1:M,each=M),phis))))
  X<-cbind(t(cos(outer(1:M,phis))),t(sin(outer(1:M,phis))))
  X_phi<-cbind(X[, rep(1:(2*M), each=M)])
  
  X_phi0<-X_phi
  Y_theta0<-Y_theta
  thetas0<-thetas
  phis0<-phis
  fhat<-(X_phi0*Y_theta0)%*%c(R,R_minus)+t(intercept%*%(cos(outer(0:M,thetas0))))
  return(fhat)
}

plot_3d <- function(out,burn = 50000,voxel_size = c(1.09375, 1.09375, 3.000000),data=NULL){
  
  iterations <- dim(out[[5]])[1]
  x<-out[[1]]
  y<-out[[2]]
  z<-out[[3]]
  value<-out[[4]]
  R_keeps=out[[5]][burn:iterations,]
  R_minus_keeps=out[[6]][burn:iterations,]
  center_keeps<-out[[9]][burn:iterations,]
  intercept_keeps<-out[[8]][burn:iterations,]
  M<-sqrt(dim(R_keeps)[2])
  center<-colMeans(center_keeps) 
  F_hat <- array(NA,dim = c(iterations-burn,10*50))
  
  for(i in 1:(iterations-burn)){
    R <- R_keeps[i,]
    R_minus <- R_minus_keeps[i,]
    intercept <- intercept_keeps[i,]
    
    theta <- seq(0, pi, length = 10)
    phi <- seq(0, 2*pi, length = 50)
    mesh <- mesh(theta, phi)
    thetas<-c(mesh$x)
    phis<-c(mesh$y)
    
    evens<-which(ceiling(1:(2*M*M)/M)%%2==0)
    odds<-c(1:(2*M*M))[-evens]
    Y<-array(dim=c(length(thetas),2*M*M))
    Y[,evens]<-t(sin(thetas)*sin(outer(1:M,thetas)))
    Y[,odds]<-t(sin(outer(1:M,thetas)))
    Y_theta<-Y
    X<-cbind(1,t(cos(outer(rep(1:M,each=M),phis))),t(sin(outer(rep(1:M,each=M),phis))))
    X<-cbind(t(cos(outer(1:M,phis))),t(sin(outer(1:M,phis))))
    X_phi<-cbind(X[, rep(1:(2*M), each=M)])
    
    X_phi0<-X_phi
    Y_theta0<-Y_theta
    thetas0<-thetas
    phis0<-phis
    
    fhat<-(X_phi0*Y_theta0)%*%c(R,R_minus)+t(intercept%*%(cos(outer(0:M,thetas0))))
    F_hat[i,] <- fhat
  }
  
  
  f_mean <- colMeans(F_hat)
  f_sd <- apply(F_hat, 2, sd)
  
  dev<-array(NA,dim=c(iterations-burn,10*50))
  for(i in 1:(iterations-burn)){
    dev[i,]<-abs(F_hat[i,]-f_mean)/f_sd
  }
  dev1<-rowMaxs(dev,value=T)
  l0<-quantile(dev1,probs=c(.95))
  
  f_u <- f_mean + l0*f_sd
  f_l <- f_mean - l0*f_sd

  f_mean_update <- update_fx(f_mean, c(mesh$x),c(mesh$y),evens,odds,M,x,y,z,center)
  f_u_update <- update_fx(f_u, c(mesh$x),c(mesh$y),evens,odds,M,x,y,z,center)
  f_l_update <- update_fx(f_l, c(mesh$x),c(mesh$y),evens,odds,M,x,y,z,center)
  
  cluster_upper = +(d2c(x,y,z,center)<=f_u_update)
  cluster_lower = +(d2c(x,y,z,center)<=f_l_update)
  cluster_mean = +(d2c(x,y,z,center)<=f_mean_update)
  
  voxel_volume_mm3 <- prod(voxel_size)
  
  volumn <- sum(cluster_mean)*voxel_volume_mm3
  volumn_ci <- c(sum(cluster_lower)*voxel_volume_mm3,sum(cluster_upper)*voxel_volume_mm3)
  
  prostate <- length(cluster_mean)*voxel_volume_mm3
  
  # now calculate the dimension 
  dim_x <- cbind.data.frame(x,y,z,cluster_mean)%>%
    filter(cluster_mean==1)%>%
    group_by(y,z)%>%
    summarise(n = n())%>%arrange(desc(n))
  
  dim_y <- cbind.data.frame(x,y,z,cluster_mean)%>%
    filter(cluster_mean==1)%>%
    group_by(x,z)%>%
    summarise(n = n())%>%arrange(desc(n))
  
  # calculating of the dimension should be modified by specific scenario
  dimension <- c(dim_x$n[1]*voxel_size[1],dim_y$n[1]*voxel_size[2],length(unique(z[cluster_mean==1]))*voxel_size[3])
  
  if(!is.null(data)){
    value_table <- list()
    names <-c('T2', 'AUGC', 'KTRANS', 'ADC', 'KEP')
    data <- data[,names]
    for(i in 1:5){
      if(names[i]=='ADC'){
        quantile_l25_upper <- quantile(as.vector(data[,i][cluster_upper==1]),0.25)
        quantile_l25_lower <- quantile(as.vector(data[,i][cluster_lower==1]),0.25)
        quantile_l25_mean <- quantile(as.vector(data[,i][cluster_mean==1]),0.25)
      }

      value_upper <- mean(as.vector(data[,i][cluster_upper==1]))
      value_lower <- mean(as.vector(data[,i][cluster_lower==1]))
      value_mean <- mean(as.vector(data[,i][cluster_mean==1]))

      if(names[i]=='ADC'){
        value_table[[names[i]]][[1]] <- data.frame(
          upper = value_upper,
          lower = value_lower,
          mean  = value_mean)
        value_table[[names[i]]][[2]] <- data.frame(
          quantile_upper = quantile_l25_upper,
          quantile_lower = quantile_l25_lower,
          quantile_mean  = quantile_l25_mean)
      }else{
        value_table[[names[i]]] <- data.frame(
          upper = value_upper,
          lower = value_lower,
          mean  = value_mean)
      }
      
    }
  }
  if(is.null(data)){
    return(list(cluster = list(post_mean=cluster_mean,lower=cluster_lower,upper=cluster_upper),surface = list(upper=f_u,mean=f_mean,lower=f_l), measurements = list(volumn=volumn,volumn_ci=volumn_ci,dimension=dimension,prostate_volumn=prostate)))
  }else{
    return(list(cluster = list(post_mean=cluster_mean,lower=cluster_lower,upper=cluster_upper),surface = list(upper=f_u,mean=f_mean,lower=f_l), measurements = list(volumn=volumn,volumn_ci=volumn_ci,dimension=dimension,prostate_volumn=prostate,apprx_value = value_table)))
  }
}
