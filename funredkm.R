funredkm = function(X,U,A,B,K,Pk,Lk,pen,lam,start){ 
  fungcv=function (lamt,lslk,snormxt,pxt2,I,J){
    llk=1+lamt*lslk
    #num = snormxt + sum(pxt2%*%(1/(llk^2)-2/llk))
    num = snormxt + sum(pxt2*(1/(llk^2)-2/llk))
    den =((I%*%J-sum(1/llk)))^2
    y=num/den;
    return(y)
  }
  laminf=10
  lamup=1e+5
  #lk=diag(Lk)
  lk<-matrix(diag(Lk), nrow=nrow(Lk), ncol=1)
  I<-nrow(U)
  G<-ncol(U)
  J<-nrow(B)
  Q<-ncol(B)
  Ij<-diag(J)
  Iq<-diag(Q)
  uI<-rep(1, length.out=I)
  uG<-rep(1, length.out=G)
  #dist<-matrix(0, G, I)
  snormx<-sum(X^2)
  if ((start==1)|(start==3)){
  KMEANS<-kmeans(X,U,0,0) 
    U<-KMEANS$U
    Xmean<-KMEANS$Xmean
    res<-KMEANS$res
  }
  su<-colSums(U)
  # start PCA
  if ((start==2)|(start==3)){
    SVD=svd(diag(1/su)%*%t(U)%*%X)
    Ps<-SVD$u
    Ls<-SVD$d
    Qs<-SVD$v
    A<-Ps[,1:Q]
    B<-Qs[,1:Q]%*%Ls[1:Q, 1:Q]
  }
  st<-sum(colSums(X^2))
  KrIqK=kronecker(Iq,K)
  lfo=Inf
  dif=st
  it=0
  eps=0.000001
  while ((dif > eps)&(it<999)){
    it=it+1
    if (G!=Q){
      SVD<-svd(B)
      l=SVD$d
      gam=max(su)*l[1]^2
      SVD=svd(A+t(U)%*%(X-U%*%A%*%t(B))%*%(B/gam))
      Ps<-SVD$u
      Ls<-SVD$d
      Qs<-SVD$v
      A=Ps%*%t(Qs)
    }
    Uold=U
    rob=colSums(t((kronecker(X,uG)-kronecker(uI,A%*%t(B)))^2))
    dist<-matrix(rob, G, I)
    ind=max.col(-t(dist))
    U=matrix(0,nrow=I,ncol=G)
    for(i in 1:I){
      U[i,ind[i]]=1
    }
   if (sum(sum(U))!=I){
     U<-matrix(0,I,G)
     for (i in 1:I) {
       m<-min(dist[,i])[,]
       p<-min(dist[,i])[,]
      U(i,p)=1
     }
     }
    su=colSums(U);
    while (sum(su==0)>0){
     indg<-which(su==0)
     indi<-which(Uold[,indg[1]]>0)
     U[indi[1],]<-Uold[indi[1],]
     su=colSums(U)
    }
    if (pen>0){
    SVD=svd (U%*%A,(nu=(min(nrow(U%*%A),ncol(U%*%A)))),nv=(min(nrow(U%*%A),ncol(U%*%A)))) 
    Pua<-SVD$u
    Lua<-SVD$d
    Qua<-SVD$v
    pxt2=matrix(t(Pk)%*%t(X)%*%Pua,nrow=Q*J)^2
    lslk<-kronecker((Lua^-2),lk)
    lam=optimise(fungcv,c(laminf,lamup),lslk=lslk,snormx=snormx,pxt2=pxt2,I=I,J=J, maximum = FALSE) 
lam=lam$minimum
        }
    iXXr=solve(kronecker(t(A)%*%diag(su)%*%A,Ij)+lam*KrIqK)
        B<-matrix(iXXr%*%matrix(t(X)%*%U%*%A, nrow =J%*%Q, ncol=1), nrow = J, ncol = Q)
    lf=sum(colSums((X-U%*%A%*%t(B))^2))+lam%*%sum(diag(t(B)%*%K%*%B));
    dif=lfo-lf
    if (pen>0) {
      dif=abs(dif)
          }
    lfo=lf      
  }
  SVD=svd (U%*%A,(nu=(min(nrow(U%*%A),ncol(U%*%A)))),nv=(min(nrow(U%*%A),ncol(U%*%A)))) 
  Pua<-SVD$u
  Lua<-SVD$d #in matlab è 2x2
  Qua<-SVD$v
  pxt2=(matrix(t(Pk)%*%t(X)%*%Pua,Q%*%J,1))^2
  lslk=kronecker((Lua)^(-2),lk)
  gcv=fungcv(lam,lslk,snormx,pxt2,I,J)
  M=A%*%t(B)
  dw=sum(colSums((X-U%*%M)^2))
  db=sum(colSums((U%*%M)^2))
  llk=1+lam*lslk
  trH=sum(1/llk)
  ssr=100%*%sum(colSums((X-U%*%A%*%t(B))^2))/st
  
  return(list(gcv=gcv,U=U,B=B,A=A,ssr=ssr,lam=lam,db=db,dw=dw,trH=trH))
  }