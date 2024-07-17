library(IndepTest)
library(expm)

#set.seed(313)

#X=data


n=nrow(X)
p=ncol(X)

X1=X[1:floor(n/2),]
X2=X[floor(n/2+1):n,]

test.stat.comp=function(X1,X2){
  # use X1 to compute mu, Sigma
  # use X2 to compute T
  n=nrow(X1)
  p=ncol(X1)
  
  mu.hat=colMeans(X1)
  Sigma.hat=cov(X1)*(n-1)/n
  
  eps=1e-6
  
  Sigma.hat.sqrt=sqrtm(Sigma.hat+diag(eps,p))
  Sigma.hat.inv=solve(Sigma.hat+diag(eps,p))
  Sigma.hat.sqrt.inv=sqrtm(Sigma.hat.inv+diag(eps,p))
  
  n=nrow(X2)
  Y=X2
  for(i in 1:n){
    Y[i,]=Sigma.hat.sqrt.inv%*%(X2[i,]-mu.hat)
  }
  U=sqrt(rowSums(Y^2))
  V=Y
  for(i in 1:n){
    V[i,]=Y[i,]/U[i]
  }
  
  
  taup=min(c(2/5,4/(4+3*p),p/4/(1+floor(p/4))))
  tau1=1/4
  
  Kp=Cp*ceiling(p*n^taup)
  K1=C1*ceiling(n^tau1)
  if(p>3){
    Wp=L2OptW(Kp,p)
  }else{
    Wp=rep(1/Kp,Kp)
  }
  W1=rep(1/K1,K1)
  
  
  U.sq=U^2
  D1=matrix(U.sq,n,n,byrow=FALSE)
  D2=matrix(U.sq,n,n,byrow=TRUE)
  D3=Y%*%t(Y)
  DY.sqr=D1+D2-2*D3
  diag(DY.sqr)=0
  DY=sqrt(DY.sqr)
  
  KNN.Y=matrix(NA,n,Kp)
  for(i in 1:n){
    RKYi=rep(0,n)
    RKYi[-i]=rank(DY[i,-i],ties.method="random")
    for(j in 1:Kp){
      KNN.Y[i,j]=which(RKYi==j)
    }
  }
  
  
  
  D4=U%*%t(U)
  DU.sqr=D1+D2-2*D4
  diag(DU.sqr)=0
  DU=sqrt(DU.sqr)
  
  
  KNN.U=matrix(NA,n,K1)
  for(i in 1:n){
    RKUi=rep(0,n)
    RKUi[-i]=rank(DU[i,-i],ties.method="random")
    for(j in 1:K1){
      KNN.U[i,j]=which(RKUi==j)
    }
  }
  
  
  Vp=pi^(p/2)/gamma(1+p/2)
  log.xiY=matrix(NA,n,Kp)
  for(i in 1:n){
    for(j in 1:Kp){
      log.xiY[i,j]=log((n-1)*max(DY[i,KNN.Y[i,j]],eps)^p*Vp/exp(digamma(j)))
    }
  }
  
  HY=0
  for(i in 1:n){
    for(j in 1:Kp){
      HY=HY+Wp[j]*log.xiY[i,j]
    }
  }
  HY=HY/n
  HY
  
  
  V1=pi^(1/2)/gamma(1+1/2)
  log.xiU=matrix(NA,n,K1)
  for(i in 1:n){
    for(j in 1:K1){
      log.xiU[i,j]=log((n-1)*max(DU[i,KNN.U[i,j]],eps)*V1/exp(digamma(j)))
    }
  }
  
  HU=0
  for(i in 1:n){
    for(j in 1:K1){
      HU=HU+W1[j]*log.xiU[i,j]
    }
  }
  HU=HU/n
  HU
  
  
  
  cp=gamma(p/2)/(2*pi^(p/2))
  
  test.stat=-HY+(p-1)*mean(log(U))+HU-log(cp)
  test.stat
  
}



test.stat=(test.stat.comp(X1,X2)+test.stat.comp(X2,X1))/2*sqrt(n)
test.stat




#KLentropy(Y,Kp,weights=TRUE)
#KLentropy(U,K1,weights=TRUE)
