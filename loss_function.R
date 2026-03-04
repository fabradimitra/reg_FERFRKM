loss_function <- function(U, C, Cbar, D, A, B, K, lambda, gamma){
  c2 <- rowSums(C*C) # vector of lenght I having squared norm of c_i as element
  cb2 <- rowSums(Cbar*Cbar) # vector of lenght G having squared norm of cbar_g as element
  D2 <- outer(c2,cb2, "+") - 2 * (C %*% t(Cbar)) # matrix of dimension IxG with squared norm of c_i - cbar_g 
  return(sum(U*D2) + 
    sum((D%*%Cbar - D%*%A%*%t(B))*(D%*%Cbar - D%*%A%*%t(B))) + 
    lambda*sum(diag(t(B)%*%K%*%B)) + gamma*sum(U*log(U)) 
)
}