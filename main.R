library(mclust)
source("kspline.R")
source("randgenuf.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FERFRKM.R")
# Simulated data to check if the model reconstruct the data generating process
set.seed(42)
I <- 100
J <- 30
Q <- 3
G <- 4
var.err <- 0.0001
# Generate U crisp
U <- randgenuc(I, G)
# U: I x G fuzzy memberships
cluster_id <- max.col(U, ties.method = "first")
# Generate A orthogonal
A <- rand_orthogonal(G, Q)
# Generate B
B <- matrix(rnorm(J * Q), nrow = J, ncol = Q) 
# Generate the error
E <- matrix(rnorm(I * J, sd = var.err), nrow = I, ncol = J)
# Compute the data matrix X
X <- U %*% A %*% t(B) + E
# Set-up for FERFRKM algorithm
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
# Run FERFRKM algorithm
res <- FERFRKM(X, K, U_init, A_init, B_init, lambda, gamma, max_iter, tol)
# Check res:
# A B'= A_hat B_hat'
ABp <- A %*% t(B)
sum((ABp - res$A %*% t(res$B))^2)/sum(ABp^2) 
# Adjusted Rand Index between true cluster and estimated cluster
cluster_id_est <- max.col(res$U, ties.method = "first")
adjustedRandIndex(cluster_id, cluster_id_est)
# Random starts
n_starts <- 10
res <- vector("list", n_starts)
errors <- numeric(n_starts)
adjusted_rand_indices <- numeric(n_starts)
for (i in seq_len(n_starts)) {
    set.seed(42 + i) # Different seed for each start
    # Random initialization of U
    U_init <- randgenuc(I, G)
    # Initialization of A and B
    A_init <- rand_orthogonal(G, Q)
    B_init <- t(t(A_init)%*%solve(t(U_init)%*%U_init)%*%t(U_init)%*%X)
    # Run FERFRKM algorithm
    res[[i]] <- FERFRKM(X, K, U_init, A_init, B_init, lambda, gamma, max_iter, tol)
    # Compute discrepancy
    errors[i] <- sum((A %*% t(B) - res[[i]]$A %*% t(res[[i]]$B))^2)/sum((A %*% t(B))^2)
    # Compute Adjusted Rand Index
    cluster_id_est <- max.col(res[[i]]$U, ties.method = "first")
    adjusted_rand_indices[i] <- adjustedRandIndex(cluster_id, cluster_id_est)
}
errors
adjusted_rand_indices

  