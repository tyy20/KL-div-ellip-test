library(IndepTest)
library(expm)

#set.seed(313)




n=nrow(X)
p=ncol(X)

eps=1e-6

# Sigma.sqrt=sqrtm(Sigma)
# Sigma.inv=solve(Sigma)
# Sigma.sqrt.inv=sqrtm(Sigma.inv)


Y=X
# for(i in 1:n){
#   Y[i,]=Sigma.sqrt.inv%*%(X[i,]-mu)
# }

### Compute U, V, Theta
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
DY.sqr[DY.sqr<0]=0
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
DU.sqr[DU.sqr<0]=0
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


KLentropy(Y,Kp,weights=TRUE)
KLentropy(U,K1,weights=TRUE)

cp=gamma(p/2)/(2*pi^(p/2))

test.stat=sqrt(n)*(-HY+(p-1)*mean(log(U))+HU-log(cp))

test.stat



