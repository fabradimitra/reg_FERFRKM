loss_function <- function(U, C, Cbar, D, A, B, K, lambda, gamma){
  cnorm2 <- rowSums(C*C) # vector of lenght I having squared norm of c_i as element
  cbnorm2 <- rowSums(Cbar*Cbar) # vector of lenght G having squared norm of cbar_g as element
  Dist2 <- outer(cnorm2,cbnorm2, "+") - 2 * (C %*% t(Cbar)) # matrix of dimension IxG with squared norm of c_i - cbar_g 
  wdev <- sum(U*Dist2) + sum((D%*%Cbar - D%*%A%*%t(B))^2)
  loss <- wdev+lambda*sum(diag(A%*%t(B)%*%K%*%B%*%t(A)))
  lossp <- loss
  if(gamma>0){
    ulu <- c(U*log(U))
    ulu[is.na(ulu)]<-0
    lossp <- loss + gamma*sum(ulu)
  }
  return(list(loss=loss,lossp=lossp, wdev = wdev))
}