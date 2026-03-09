FEFRKM <- function(C,K,U,A,B,lambda,gamma,max_iter = Inf,tol = 1e-6){
  # This function implements the FEFRKM algorithm for fuzzy functional data clustering.
  # Inputs:
  #  C: (I x J) coefficients for the natural cubic spline basis functions for each individual i
  #  K: (J x J) inner products of second derivative of the natural cubic spline basis functions 
  #  U: (I x G) initial fuzzy membership matrix
  #  A: (G x Q) initial coefficients of the cluster centroid expansion in the reduced subspace
  #  B: (J x Q) initial coefficients of the basis of the reduced subspace wrt the natural cubic 
  #  spline basis functions
  #  lambda: roughness penalty parameter for the cluster centroids
  #  gamma: degree of fuzziness
  #  max_iter: maximum number of iterations
  #  tol: tolerance for convergence
  #######################################################################################
  J <- nrow(B)
  Q <- ncol(B)
  fab <- colSums(U) 
  D <- diag(fab^(-1/2))
  Cbar <- (D^(-2)) %*% t(U) %*% C
  loss_function_curr <- loss_function(U, C, Cbar, D, A, B, K, lambda, gamma)
  cat("Initial loss function value: ", loss_function_curr, "\n", file = "log/log.txt", append = TRUE)
  dif <- Inf
  iter <- 0
  while(dif > tol & iter < max_iter){
    iter <- iter + 1
    # Update U
    cnorm2 <- rowSums(C*C) # vector of lenght I having squared norm of c_i as element
    aBnorm2 <- rowSums((B %*% t(A))*(B %*% t(A))) # vector of lenght G having squared norm of B a_g as element
    distxmg2 <- outer(cnorm2,aBnorm2, "+") - 2 * (C %*% B %*% t(A)) # matrix of dimension IxG ||c_i - B a_g||^2
    fab <- exp(-distxmg2/gamma)
    U <- fab/colSums(fab)
    #######################################################################################
    # Update B
    fab <- solve(kronecker(t(A) %*% (D^2) %*% A, diag(numeric(J)+1)) + 
      lambda*kronecker(numeric(Q)+1,K)) %*% as.vector(t(Cbar) %*% (D^2) %*% A) # Ridge regression type update for B
    B <- matrix(fab, nrow = J, ncol = Q)
    #######################################################################################
    # Update A
    FAB <- kronecker(D^2,t(B) %*% B)
    alpha_max <- max(eigen(FAB, symmetric = TRUE, only.values = TRUE)$values) # maximum eigenvalue to be inserted in Kiers maj. method
    matr_kiers <- A + ((D^2) %*% Cbar %*% B - (D^2) %*% A %*% t(B) %*% B)/alpha_max # matrix of Kiers maj. method
    S <- svd(matr_kiers)
    A <- S$u %*% t(S$v)
    #######################################################################################
    # Compute loss function and check convergence
    loss_function_new <- loss_function(U, C, Cbar, D, A, B, K, lambda, gamma)
    dif <-  loss_function_new - loss_function_curr
    loss_function_curr <- loss_function_new
    cat("Iteration: ", iter, " Loss function value: ", loss_function_curr, " Difference: ", dif, "\n", file = "log/log.txt", append = TRUE)
  }
} 