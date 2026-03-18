source("kspline.R")
source("randgenuf.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FEFRKM.R")
# Simulated data to check if the model reconstruct the data generating process
set.seed(42)
I <- 100
J <- 30
Q <- 3
G <- 4
var.err <- 0.1
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
K <- kspline(1:J)
lambda <- 1
gamma <- 1
max_iter <- Inf
tol <- 1e-6
# K-means initialization of U
kmeans_res <- kmeans(X, G, nstart = 10)
U_init <- matrix(0, nrow = I, ncol = G)
U_init[cbind(seq_len(I), kmeans_res$cluster)] <- 1L
# PCA initialization of A and B
SVD <- svd(diag(1/colSums(U_init)) %*% t(U_init) %*% X)
A_init <- SVD$u[, 1:Q]
B_init <- SVD$v[, 1:Q] %*% diag(SVD$d[1:Q])
# Run FEFRKM algorithm
res <- FEFRKM(X, K, U_init, A_init, B_init, lambda, gamma, max_iter, tol)
# Check results 
mean(abs(A %*% t(B) - res$A %*% t(res$B))) # A B'= A_hat B_hat' should be small
# Random starts
n_starts <- 10
results <- vector("list", n_starts)
discrepancies <- numeric(n_starts)
for (i in seq_len(n_starts)) {
    set.seed(42 + i) # Different seed for each start
    # Random initialization of U
    U_init <- randgenuc(I, G)
    # Initialization of A and B
    A_init <- rand_orthogonal(G, Q)
    B_init <- t(t(A_init)%*%solve(t(U_init)%*%U_init)%*%t(U_init)%*%X)
    # Run FEFRKM algorithm
    results[[i]] <- FEFRKM(X, K, U_init, A_init, B_init, lambda, gamma, max_iter, tol)
    # Compute discrepancy
    discrepancies[i] <- mean(abs(A %*% t(B) - results[[i]]$A %*% t(results[[i]]$B)))
}
discrepancies

  