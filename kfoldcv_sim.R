require(doParallel)
require(parallel)
require(foreach)
require(mclust)
require(clue)
#
source("kspline.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FERFRKM.R")
source("perm_hungarian_fast.R")
source("cv_ferfrkm_XB.R")
source("make_ferfrkm_grid.R")
source("init_ferfrkm.R")
source("make_folds.R")
# Simulation preparation -----
randomstarts <- 5
kmeans_starts <- 20
# Set up dimensions and centroids
I <- 50
J <- 101
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
t_grid <- seq(-1, 1, length.out = J)
f1 <- psi1_wiggly(t_grid)
f2 <- psi2_wiggly(t_grid)
# Cluster centroids
curves <- apply(A, 1, function(a) a[1] * f1 + a[2] * f2)
res <- kspline(t_grid)
K <- res$K
Pk <- res$Pk
Lk <- res$Lk
# Grid of lambda and gamma
default_gamma_grid <- seq(0.2, 3, by = 0.20)
default_lambda_grid <- 10^seq(-7, -3, by = 0.25)
#
simulation_results <- data.frame(
  lambda_best = numeric(250),
  gamma_best = numeric(250),
  ARI = numeric(250),
  SSQerr = numeric(250)
)
# Monte Carlo simulations
for(iter in c(3:4)){
  set.seed(iter)
  U <- randgenuc(I, G)
  cluster_labels <- max.col(U, ties.method = "first")
  # Generate the error
  E <- matrix(rnorm(I * J, sd = var.err), nrow = I, ncol = J)
  # Compute the data matrix X
  X <- U %*% t(curves) + E
  # Cross validation
  cv_res <- cv_ferfrkm_XB(
    Xtr = X,
    G = G,
    Q = Q,
    K = K,
    Pk = Pk,
    Lk = Lk,
    folds = 5,
    gamma_grid = default_gamma_grid,
    lambda_grid = default_lambda_grid,
    max_iter = Inf,
    tol = 1e-8,
    nstart_kmeans = kmeans_starts,
    parallel = TRUE,
    seed = iter,
    randomstarts = randomstarts
  )
  simulation_results$lambda_best[iter] <- cv_res$best$lambda
  simulation_results$gamma_best[iter] <- cv_res$best$gamma
  # Fit the best combination:
  cur_loss <- Inf
  for(start in seq_len(randomstarts)){
    if(start == 1){
      init <- init_ferfrkm(X, G, Q, seed = iter, nstart_kmeans = kmeans_starts) 
    }else{
      U_init <- randgenuc(I, G)
      A_init <- rand_orthogonal(G, Q)
      B_init <- t(t(A_init)%*%solve(t(U_init)%*%U_init)%*%t(U_init)%*%X)
      init <- list(U=U_init, A=A_init, B=B_init)
    }
    # Run FERFRKM algorithm
    res_cur <- tryCatch(
      FERFRKM(C=X,
              K=K,
              Pk=Pk,
              Lk=Lk,
              U=init$U,
              A=init$A,
              B=init$B,
              lambda= cv_res$best$lambda,
              gamma = cv_res$best$gamma,
              max_iter = Inf,
              tol = 1e-8),
          error = function(e) NULL
        )
      if (is.null(res_cur)) {
          next
        }
    if(cur_loss>res_cur$loss_function){
      res <- res_cur
      cur_loss <- res$loss_function
    }
  }
  cluster_labels_est <- max.col(res$U, ties.method = "first")
  simulation_results$ARI[iter] <- adjustedRandIndex(cluster_labels,cluster_labels_est)
  ABp <- res$A %*% t(res$B)
  tcurves <- t(curves)
  # Permute the estimated curves to match the true curves
  perm <- perm_hungarian_fast(tcurves, ABp, J)
  simulation_results$SSQerr[iter] <- sum((ABp[perm,] - t(curves))^2)
  cat("End Monte Carlo Simulation: ", iter, " Lambda* ", simulation_results$lambda_best[iter], " Gamma* ", simulation_results$gamma_best[iter],
    " ARI: ", simulation_results$ARI[iter], " sSqErr: ", simulation_results$SSQerr[iter], "\n")
}
# Results
mean(simulation_results$ARI)
sd(simulation_results$ARI)
mean(simulation_results$SSQerr)
sd(simulation_results$SSQerr)
# Plot the centroids and their reconstruction for one iteration ----
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
cols <- c("red", "blue", "darkgreen", "orange")[cluster_labels]
Y <- t(apply(X, 1, function(y) splinefun(t_grid, y, method = "natural")(tt)))
matplot(tt, t(Y), type = "l", lty = 1, col = cols, lwd = 2,
        xlab = "t", ylab = "spline value")
legend("bottomright", legend = paste("label", 1:4),
       col = c("red","blue","darkgreen","orange"), lwd = 2, bty = "n")
