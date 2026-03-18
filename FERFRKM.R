FERFRKM <- function(C,K,U,A,B,lambda,gamma,max_iter = Inf,tol = 1e-6){
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
  I <- nrow(C) 
  J <- ncol(C) 
  G <- ncol(U) 
  Q <- ncol(B)
  stopifnot(is.matrix(C), is.matrix(K), is.matrix(U), is.matrix(A), is.matrix(B))
  if (!all(dim(K) == c(J, J))) stop("K must be J x J with J = ncol(C).")
  if (!all(dim(U) == c(I, G))) stop("U must be I x G with I = nrow(C).")
  if (!all(dim(A) == c(G, Q))) stop("A must be G x Q with G = ncol(U).")
  if (!all(dim(B) == c(J, Q))) stop("B must be J x Q with J = ncol(C) and Q = ncol(A).")
  if (!isTRUE(all.equal(K, t(K), tol = 1e-10))) stop("K must be symmetric.")
  if (lambda < 0 || gamma <= 0 || tol <= 0 || max_iter < 1) stop("Invalid hyperparameters.")
  if (any(!is.finite(C)) || any(!is.finite(K)) || any(!is.finite(U)) || any(!is.finite(A)) || any(!is.finite(B))) {
  stop("Inputs contain NA/NaN/Inf.")
  }
  cnorm2 <- rowSums(C*C) # vector of lenght J having squared norm of c_i as element
  IJ <- diag(J)
  IQkK <- kronecker(diag(Q),K)
  loss_function_curr <- Inf
  dif <- Inf
  iter <- 0
  while(dif > tol & iter < max_iter){
    iter <- iter + 1
    # Update U
    BAp <- B %*% t(A)
    aBnorm2 <- colSums(BAp*BAp) # vector of lenght G having squared norm of B a_g as element
    distxmg2 <- outer(cnorm2,aBnorm2, "+") - 2 * (C %*% BAp) # matrix of dimension IxG ||c_i - B a_g||^2
    U <- exp(-distxmg2/gamma)
    U <- U/matrix(rowSums(U),I,G)
    #######################################################################################
    # Update B
    D <- diag(sqrt(colSums(U)))
    D2 <- D^2
    Cbar <- diag(1/diag(D2)) %*% t(U) %*% C
    B <- solve(kronecker(t(A) %*% (D2) %*% A, IJ) + lambda*IQkK) %*% c(t(Cbar) %*% D2 %*% A) # Ridge regression type update for B
    B <- matrix(B, nrow = J, ncol = Q)
    #######################################################################################
    # Update A
    Fab <- kronecker(D2,t(B) %*% B)
    alpha_max <- max(eigen(Fab, symmetric = TRUE, only.values = TRUE)$values) # maximum eigenvalue to be inserted in Kiers maj. method
    matr_kiers <- A + (D2 %*% Cbar %*% B - D2 %*% A %*% t(B) %*% B)/alpha_max # matrix of Kiers maj. method
    S <- svd(matr_kiers)
    A <- S$u %*% t(S$v)
    #######################################################################################
    # Compute loss function and check convergence
    loss_function_new <- loss_function(U, C, Cbar, D, A, B, K, lambda, gamma)
    dif <-  loss_function_curr - loss_function_new
    loss_function_curr <- loss_function_new
    cat("Iteration: ", iter, " Loss function value: ", loss_function_curr, " Difference: ", dif, "\n") #, file = "log/log.txt", append = TRUE)
  }
  return(list(U = U, A = A, B = B, loss_function = loss_function_curr, iterations = iter))
} 