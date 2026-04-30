fungcv <- function (lambda,lmsk,snormCbarpD,CbarpDP,G,J){
    lLK<-1+lambda*lmsk
    num <- snormCbarpD + sum((CbarpDP^2)*(1/(lLK^2)-2/lLK))
    den <-((G*J-sum(1/lLK)))^2
    return(num/den)
  }