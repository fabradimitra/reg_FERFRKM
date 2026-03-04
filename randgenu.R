randgenu<-function(n,k){
#randomly generate a membership matrix U
U=matrix(runif(n*k),nrow=n,ncol=k)
U=diag(1/colSums(t(U)))%*%U;
return(U)
}
