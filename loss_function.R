loss_function <- function(U, C, Cbar, D, A, B, K, lambda, gamma){
  cnorm2 <- rowSums(C*C) # vector of lenght I having squared norm of c_i as element
  cbnorm2 <- rowSums(Cbar*Cbar) # vector of lenght G having squared norm of cbar_g as element
  D2 <- outer(cnorm2,cbnorm2, "+") - 2 * (C %*% t(Cbar)) # matrix of dimension IxG with squared norm of c_i - cbar_g 
  loss <- sum(U*D2) + sum((D%*%Cbar - D%*%A%*%t(B))*(D%*%Cbar - D%*%A%*%t(B))) + gamma*sum(U*log(U))
  lossp <- loss + lambda*sum(diag(A%*%t(B)%*%K%*%B%*%t(A)))
  return(list(loss=loss,lossp=lossp))
}