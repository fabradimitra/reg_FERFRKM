source("kspline.R")
source("randgenuf.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FEFRKM.R")
# Simulated data example
set.seed(200)
I <-100; J <- 30; Q <- 2; G <- 3
X <- matrix(rnorm(I*J),I,J)
K <- kspline(1:J)
U <- randgenuf(I,G)
A <- rand_orthogonal(G,Q)
B <- t(t(A)%*%solve(t(U)%*%U)%*%t(U)%*%X)
#
FEFRKM(X,K,U,A,B,lambda=10,gamma=1,max_iter = Inf,tol = 1e-6)
# Simulated data to check if the model reconstruct the data generating process
set.seed(42)
I <- 100
J <- 40
Q <- 3
G <- 4
var.err <- 0.5
# Generate U crisp
U <- randgenuc(I, G)
# Generate A orthogonal
A <- rand_orthogonal(G, Q)
# Generate B
B <- matrix(rnorm(J * Q), nrow = J, ncol = Q)
# Generate the error
E <- matrix(rnorm(I * J, sd = var.err), nrow = I, ncol = J)
# Compute the data matrix X
X <- U %*% A %*% t(B) + E
# Set-up for FEFRKM algorithm
set.seed(123)
K <- kspline(1:J)
# Starting values
U_init <- randgenuf(I,G)
A_init <- rand_orthogonal(G,Q)
B_init <- t(t(A)%*%solve(t(U)%*%U)%*%t(U)%*%X)
lambda <- 10
gamma <- 1
max_iter <- Inf
tol <- 1e-6
# Run FEFRKM algorithm
res <- FEFRKM(X, K, U_init, A_init, B_init, lambda, gamma, max_iter, tol)
# Check results 
# A%*%t(B) should be close to res$A%*%t(res$B)
mean((A %*% t(B) - res$A %*% t(res$B))^2) # should be small