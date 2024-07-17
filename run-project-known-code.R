C1=1
Cp=1
B=100



X=data

Y=X
U=sqrt(rowSums(Y^2))


test.stat.resamp=rep(NA,B)
for(b in 1:B){
  Yb=matrix(rnorm(n*p),n,p)
  for(i in 1:n){
    Yb[i,]=Yb[i,]/sqrt(sum(Yb[i,]^2))*U[i]
  }
  X=Yb
  source("test-stat-known.R")
  
  test.stat.resamp[b]=test.stat
}


Tb.bar=mean(test.stat.resamp)

X=data
source("test-stat-known.R")
T.hat=test.stat
test.stat=T.hat-Tb.bar
source("asymp-known-cs3-2.R")