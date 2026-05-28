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
source("CV_FERFRKM.R")
source("init_FERFRKM.R")
source("make_folds.R")
# Simulation preparation -----
randomstarts <- 5
kmeans_starts <- 20
# Set up dimensions and centroids
I <- 150
J <- 101
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
t_grid <- seq(-1, 1, length.out = J)
f1 <- psi1_wiggly(t_grid)
f2 <- psi2_wiggly(t_grid)
# Cluster centroids
curves <- apply(A, 1, function(a) a[1] * f1 + a[2] * f2)
res <- kspline(t_grid)
K <- res$K
Pk <- res$Pk
Lk <- res$Lk
#
simulation_results <- data.frame(
  lambda_best = numeric(250),
  gamma_best = numeric(250),
  ARI = numeric(250),
  SSQerr = numeric(250)
)
# Monte Carlo simulations
for(iter in c(1)){
  set.seed(iter)
  U <- randgenuc(I, G)
  cluster_labels <- max.col(U, ties.method = "first")
  # Generate the error
  E <- matrix(rnorm(I * J, sd = var.err), nrow = I, ncol = J)
  # Compute the data matrix X
  X <- U %*% t(curves) + E
  # Cross validation
  invisible(capture.output(
    cv_res <- CV_FERFRKM(
      Xtr = X,
      G = G,
      Q = Q,
      K = K,
      Pk = Pk,
      Lk = Lk,
      lambda_init = 0.001,
      gamma_init = 1,
      folds = 5,
      max_iter = Inf,
      tol = 1e-8,
      nstart_kmeans = kmeans_starts,
      seed = iter,
      randomstarts = randomstarts
    ),
    type = "output"
  ))
  simulation_results$lambda_best[iter] <- cv_res$par[1]
  simulation_results$gamma_best[iter] <- cv_res$par[2]
  # Fit the best combination:
  cur_loss <- Inf
  for(start in seq_len(randomstarts)){
    if(start == 1){
      init <- init_FERFRKM(X, G, Q, seed = iter, nstart_kmeans = kmeans_starts) 
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
              lambda= cv_res$par[1],
              gamma = cv_res$par[2],
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
