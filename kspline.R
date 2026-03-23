kspline = function(t){ 
 #t<-seq(1,10,1)
  T<-length(t)
  h<-matrix(0,T,1)
  Q<-matrix(0,T,T-2)
  R<-matrix(0,T-2,T-2)

  for (i in 1:(T-1)) {
    h[i]=t[i+1]-t[i]
    }
for (j in 2:(T-1)) {
  Q[j-1,j-1]=1/h[j-1]
  Q[j,j-1]=-1/h[j-1]-1/h[j]
  Q[j+1,j-1]=1/h[j]
}
    for (i in 2:(T-2)) {
      R[i-1,i-1]=(h[i-1]+h[i])/3
      R[i-1,i]=h[i]/6
      R[i,i-1]=h[i]/6 
    }
  R[T-2,T-2]=(h[T-2]+h[T-1])/3;
  K<-Q%*%solve(R)%*%t(Q)    
  return(K)
}