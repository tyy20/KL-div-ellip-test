alpha=0.05




n=nrow(X)
p=ncol(X)

mu.hat=colMeans(X)
Sigma.hat=cov(X)*(n-1)/n

eps=1e-6

Sigma.hat.sqrt=sqrtm(Sigma.hat+diag(eps,p))
Sigma.hat.inv=solve(Sigma.hat+diag(eps,p))
Sigma.hat.sqrt.inv=sqrtm(Sigma.hat.inv+diag(eps,p))

Y=X
for(i in 1:n){
  Y[i,]=Sigma.hat.sqrt.inv%*%(X[i,]-mu.hat)
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











Ehat.logU=mean(log(U))


Sigma.influence=solve(Sigma.hat.sqrt%x%Sigma.hat+Sigma.hat%x%Sigma.hat.sqrt+diag(eps,p^2))
tr=function(A){
  sum(diag(A))
}


psi.Sigma=array(NA,dim=c(n,p,p))
psi.Sigma.invsqrt=array(NA,dim=c(n,p,p))

for(i in 1:n){
  psi.Sigma[i,,]=(X[i,]-mu.hat)%*%t(X[i,]-mu.hat)-Sigma.hat
  psi.Sigma.invsqrt[i,,]=-matrix(Sigma.influence%*%matrix(psi.Sigma[i,,],p^2,1),p,p)
}



VY=0
for(i in 1:n){
  VY.add=0
  for(j in 1:Kp){
    VY.add=VY.add+Wp[j]*log.xiY[i,j]
  }
  VY=VY+(VY.add-HY+tr(psi.Sigma.invsqrt[i,,]%*%Sigma.hat.sqrt))^2
}
VY=VY/n





term1=matrix(0,p,p)
for(j in 1:n){
  term1=term1+V[j,]%*%t(V[j,])
}
term1=term1/n

term4=matrix(0,1,p)
for(j in 1:n){
  term4=term4+t(V[j,])/U[j]
}
term4=term4/n

VlogU=0
for(i in 1:n){
  VlogU=VlogU+(tr(term1%*%psi.Sigma.invsqrt[i,,]%*%Sigma.hat.sqrt)-term4%*%Y[i,]+
                 log(U[i])-Ehat.logU)^2
}
VlogU=VlogU/n





h=n^(-1/5)

K.mat=matrix(NA,n,n)
K.der=matrix(NA,n,n)
for(i in 1:n){
  for(j in i:n){
    DUij=U[i]-U[j]
    K.mat[i,j]=1/sqrt(2*pi)*(1/h)*exp(-DUij^2/(2*h^2))
    K.mat[j,i]=K.mat[i,j]
    K.der[i,j]=-1/sqrt(2*pi)*(1/h^2)*exp(-DUij^2/(2*h^2))*DUij
    K.der[j,i]=-K.der[i,j]
  }
}

fhat.val=rowMeans(K.mat)
fhat.der=rowMeans(K.der)


term3=matrix(0,p,p)
for(j in 1:n){
  term3=term3+fhat.der[j]/fhat.val[j]*U[j]*V[j,]%*%t(V[j,])
}
term3=term3/n

term5=matrix(0,1,p)
for(j in 1:n){
  term5=term5+fhat.der[j]/fhat.val[j]*t(V[j,])
}
term5=term5/n


VU=0
for(i in 1:n){
  VU.add=0
  for(j in 1:K1){
    VU.add=VU.add+W1[j]*log.xiU[i,j]
  }
  VU=VU+(VU.add-HU-tr(term3%*%psi.Sigma.invsqrt[i,,]%*%Sigma.hat.sqrt)+term5%*%Y[i,])^2
}
VU=VU/n


critical.val=test.stat/sqrt(n)-qnorm(alpha,lower.tail=FALSE)*
  (sqrt(3)*sqrt(VY+VU+(p-1)^2*VlogU))/sqrt(n)
critical.val
decision=ifelse(0<critical.val,1,0)
decision
p.val=pnorm(test.stat/(sqrt(3)*sqrt(VY+VU+(p-1)^2*VlogU)),
            lower.tail=FALSE)
p.val

result=list(test.stat=test.stat,VY=VY,VU=VU,VlogU=VlogU,
            critical.val=critical.val,decision=decision,p.val=p.val)
