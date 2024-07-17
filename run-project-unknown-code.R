C1=1
Cp=1
B=100


X=data

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




test.stat.resamp=rep(NA,B)
for(b in 1:B){
  Yb=matrix(rnorm(n*p),n,p)
  Xb=matrix(NA,n,p)
  for(i in 1:n){
    Yb[i,]=Yb[i,]/sqrt(sum(Yb[i,]^2))*U[i]
    Xb[i,]=Sigma.hat.sqrt%*%Yb[i,]+mu.hat
  }
  X=Xb
  source("test-stat-unknown.R")
  
  test.stat.resamp[b]=test.stat
}


Tb.bar=mean(test.stat.resamp)

X=data
source("test-stat-unknown.R")
T.hat=test.stat
test.stat=T.hat-Tb.bar
source("asymp-unknown-cs3-2.R")