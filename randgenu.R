randgenu=function(n,k){
#randomly generate a membership matrix U
U=matrix(runif(n*k),nrow=n,ncol=k)
U=diag(1/colSums(t(U)))%*%U;
ind=max.col(U)
U=matrix(0,nrow=n,ncol=k)
for(i in 1:n){
  U[i,ind[i]]=1
}
return(U)
}