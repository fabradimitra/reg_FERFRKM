library(mclust)
library(splines)
library(clue)
source("kspline.R")
source("randgenuf.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FERFRKM.R")
source("perm_hungarian_fast.R")
# Simulation parameters
I <- 150
J <- 100
Q <- 2
G <- 4
var.err <- 0.2
# smooth smooth
psi1_smooth <- function(t) {
  t + sin(pi * t) * exp(-t)
}
psi2_smooth <- function(t) {
  cos(3 + pi * t)
}
psi1_wiggly <- function(t) {
  cos(20 * t)
}
psi2_wiggly <- function(t) {
  sin(20 * t)
}
# True A matrix (orthogonal)
A <- matrix(c(1,0,1,-1,0,1,1,1), nrow= G, ncol = Q)
# Evaluate the curves at a grid of observed points
t_grid <- seq(0.1, 1, length.out = J)
f1 <- psi1_smooth(t_grid)
f2 <- psi2_smooth(t_grid)
# Cluster centroids
curves <- apply(A, 1, function(a) a[1] * f1 + a[2] * f2)
res <- kspline(t_grid)
K <- res$K
Pk <- res$Pk
Lk <- res$Lk
# Hyperparameters for FERFRKM
lambda <- 0.01
gamma <- 1
max_iter <- Inf
tol <- 1e-8
random_init <- FALSE  
adjustedRandIndices <- numeric(250)
sSqErrors <- numeric(250)
# Simulation loop
for(iter in c(1:250)){
  set.seed(iter)
  # Data generation
  U <- randgenuc(I, G)
  cluster_labels <- max.col(U, ties.method = "first")
  # Generate the error
  E <- matrix(rnorm(I * J, sd = var.err), nrow = I, ncol = J)
  # Compute the data matrix X
  X <- U %*% t(curves) + E
  # Run the FERFRKM algorithm
  # K-means initialization of U
  kmeans_res <- kmeans(X, G, nstart = 10)
  U_init <- matrix(0, nrow = I, ncol = G)
  U_init[cbind(seq_len(I), kmeans_res$cluster)] <- 1L
  # PCA initialization of A and B
  SVD <- svd(diag(1/colSums(U_init)) %*% t(U_init) %*% X)
  A_init <- SVD$u[, 1:Q]
  B_init <- SVD$v[, 1:Q] %*% diag(SVD$d[1:Q])
  # Run FERFRKM algorithm
  res <- FERFRKM(C=X,
                 K=K,
                 Pk=Pk,
                 Lk=Lk,
                 U=U_init,
                 A=A_init,
                 B=B_init,
                 lambda=lambda,
                 gamma = gamma,
                 max_iter = max_iter,
                 tol = tol) 
  cluster_labels_est <- max.col(res$U, ties.method = "first")
  adjustedRandIndices[iter] <- adjustedRandIndex(cluster_labels,cluster_labels_est)
  ABp <- res$A %*% t(res$B)
  tcurves <- t(curves)
  # Permute the estimated curves to match the true curves
  perm <- perm_hungarian_fast(tcurves, ABp, J)
  sSqErrors[iter] <- sum((ABp[perm,] - t(curves))^2)
}
mean(adjustedRandIndices)
mean(sSqErrors)
# Plotting utilities
# Plot the centroids and their reconstruction
tt <- seq(min(t_grid), max(t_grid), length.out = 400)
Ym <- apply(t(curves), 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
par(mfrow = c(1,2))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
Ymr <- apply(ABp[perm,], 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
legend(
  "bottomright",
  legend = paste0("cluster ", 1:4),
  col = c("red","blue","darkgreen","orange"),
  lwd = 2, bty = "n"
)
# Plot the data points colored by their true cluster
cols <- c("red", "blue", "darkgreen", "orange")[cluster_labels]
Y <- t(apply(X, 1, function(y) splinefun(t_grid, y, method = "natural")(tt)))
matplot(tt, t(Y), type = "l", lty = 1, col = cols, lwd = 2,
        xlab = "t", ylab = "spline value")
legend("bottomright", legend = paste("label", 1:4),
       col = c("red","blue","darkgreen","orange"), lwd = 2, bty = "n")


  

