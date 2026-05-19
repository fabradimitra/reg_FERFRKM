require(doParallel)
require(parallel)
require(foreach)
require(mclust)
require(clue)
#
source("kspline.R")
source("randgenuf.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FERFRKM.R")
source("perm_hungarian_fast.R")
source("cv_ferfrkm_XB.R")
#
seed <- 123
set.seed(seed)
# Data generation -----
# Simulation parameters
I <- 150
J <- 100
Q <- 2
G <- 4
var.err <- 0.3
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
f1 <- psi1_wiggly(t_grid)
f2 <- psi2_smooth(t_grid)
# Cluster centroids
curves <- apply(A, 1, function(a) a[1] * f1 + a[2] * f2)
res <- kspline(t_grid)
K <- res$K
Pk <- res$Pk
Lk <- res$Lk
U <- randgenuc(I, G)
cluster_labels <- max.col(U, ties.method = "first")
# Generate the error
E <- matrix(rnorm(I * J, sd = var.err), nrow = I, ncol = J)
# Compute the data matrix X
X <- U %*% t(curves) + E
# Cross validation -----
default_gamma_grid <- seq(0.5, 10, by = 0.5)
default_lambda_grid <- 10^seq(-5, 0, by = 0.25)
# Make grid funtion
make_ferfrkm_grid <- function(
  gamma_grid = default_gamma_grid,
  lambda_grid = default_lambda_grid
) {
  expand.grid(
    gamma = gamma_grid,
    lambda = lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}
# Make fold function
make_folds <- function(n, k = 10, seed = 123) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}
# Initialize starting values for the algorithm
init_ferfrkm <- function(X, G, Q, seed = NULL, nstart_kmeans = 10) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  kmeans_res <- kmeans(X, G, nstart = nstart_kmeans)
  U_init <- matrix(0, nrow = nrow(X), ncol = G)
  U_init[cbind(seq_len(nrow(X)), kmeans_res$cluster)] <- 1

  svd_init <- svd(diag(1 / colSums(U_init)) %*% t(U_init) %*% X)
  A_init <- svd_init$u[, seq_len(Q), drop = FALSE] # drop = FALSE forces R to keep a matrix
  B_init <- svd_init$v[, seq_len(Q), drop = FALSE] %*% diag(svd_init$d[seq_len(Q)])

  list(U = U_init, A = A_init, B = B_init)
}
#
cv_res <- cv_ferfrkm_XB(
  Xtr = X,
  G = 4,
  Q = 2,
  K = K,
  Pk = Pk,
  Lk = Lk,
  folds = 5,
  gamma_grid = default_gamma_grid,
  lambda_grid = default_lambda_grid,
  max_iter = Inf,
  tol = 1e-8,
  nstart_kmeans = 10,
  parallel = TRUE,
  seed = seed
)
cv_res$best
# Fit the best combination:
init <- init_ferfrkm(X, G, Q, seed = seed, nstart_kmeans = 10)
# Run FERFRKM algorithm
res <- FERFRKM(C=X,
               K=K,
               Pk=Pk,
               Lk=Lk,
               U=init$U,
               A=init$A,
               B=init$B,
               lambda= cv_res$best$lambda,
               gamma = cv_res$best$gamma,
               max_iter = Inf,
               tol = 1e-8) 
cluster_labels_est <- max.col(res$U, ties.method = "first")
ARI <- adjustedRandIndex(cluster_labels,cluster_labels_est)
ABp <- res$A %*% t(res$B)
tcurves <- t(curves)
# Permute the estimated curves to match the true curves
perm <- perm_hungarian_fast(tcurves, ABp, J)
sSq <- sum((ABp[perm,] - t(curves))^2)
# Plot the centroids and their reconstruction
tt <- seq(min(t_grid), max(t_grid), length.out = 400)
Ym <- apply(t(curves), 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
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

