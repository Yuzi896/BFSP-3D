library(fields);library(extraDistr);library(Matrix);library(sparseMVN);library(FastGP)

generate_simulation_data<-function(shape,smoothness1,tau){
  loc<-as.matrix(expand.grid(seq(0,1,length=15),seq(0,1,length=15),seq(0,1,length=15)))
  x<-loc[,1]
  y<-loc[,2]
  z<-loc[,3]
  
  if(shape=="ellipse"){
    clus<-which((x-.5)^2+(y-.5)^2+.4*(z-.5)^2<.06)
  }
  
  if(shape=="rectangle"){
    clus<-which(x<.8&x>.2&y<.6&y>.3&z>.3&z<.6)
  }

  if(shape=="X"){
    clus<-which(x<(-z+1.1)&x>(-z+.9)&z>.2&z<.8&y>.3&y<.7|x<(z+.1)&x>(z-.1)&z>.2&z<.8&y>.3&y<.7)
  }
  if(shape=="cone"){
    clus<-which((x-.5)^2/.5^2+(y-.5)^2/.5^2<(z-.5)^2/.5^2&z<.5&z>.1)
  }
  if(shape=="cross"){
    clus<-which(x<.6&x>.4&y>.2&y<.5&z<.8&z>.2|x>.2&x<.8&y>.2&y<.5&z>.4&z<.6)
  }
  
  
  loc0<-loc[-clus,]
  loc0<-loc0[order(loc0[,1],loc0[,2],loc0[,3],decreasing=FALSE),]
  n0=dim(loc0)[1]
  loc1<-loc[clus,]
  n1<-dim(loc1)[1]
  
  
  Cov0<-stationary.cov(loc0, x2=NULL, Covariance = "Matern", theta=.05,smoothness=1)
  Y0     <- rcpp_rmvnorm(1,Cov0+tau*diag(n0),rep(0,n0))
  
  Cov1<-stationary.cov(loc1, x2=NULL, Covariance = "Matern", theta=.05,smoothness=smoothness1)
  Y1    <- rcpp_rmvnorm(1,Cov1+tau*diag(n1),rep(2,n1))
  
  value<-array(dim=c(length(x),1))
  
  value[-clus]<-Y0
  value[clus]<-Y1
  output <- list(x,y,z,value,clus) #,c(smoothness1,tau))
  return(output)
}

# run the simulation
shape = "ellipse"
smooth = 1.5
tau = 0.5
sim <- generate_simulation_data(shape,smooth,tau)
