alpha=0.05

Ehat.logU=mean(log(U))

VY=0
for(i in 1:n){
  VY.add=0
  for(j in 1:Kp){
    VY.add=VY.add+Wp[j]*log.xiY[i,j]
  }
  VY=VY+(VY.add-HY)^2
}
VY=VY/n

VU=0
for(i in 1:n){
  VU.add=0
  for(j in 1:K1){
    VU.add=VU.add+W1[j]*log.xiU[i,j]
  }
  VU=VU+(VU.add-HU)^2
}
VU=VU/n


VlogU=(n-1)/n*var(log(U))


critical.val=test.stat/sqrt(n)-qnorm(alpha,lower.tail=FALSE)*sqrt(3)*sqrt(VY+VU+(p-1)^2*VlogU)/sqrt(n)
critical.val
decision=ifelse(0<critical.val,1,0)
decision
p.val=pnorm(test.stat/(sqrt(3)*sqrt(VY+VU+(p-1)^2*VlogU)),lower.tail=FALSE)
p.val

result=list(test.stat=test.stat,VY=VY,VU=VU,critical.val=critical.val,decision=decision,p.val=p.val)
