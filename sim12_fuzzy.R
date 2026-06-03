require(doParallel)
require(parallel)
require(foreach)
require(fclust)
require(clue)
require(MASS)
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
randomstarts_cv <- 3
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
#
simulation_results <- data.frame(
  lambda_best = numeric(250),
  gamma_best = numeric(250),
  FARI = numeric(250),
  SSQerr = numeric(250)
)
IJ <- diag(J)
# Monte Carlo simulations
for(iter in c(1)){
  set.seed(iter)
  # Simulate labels 
  dummy_labels <- t(rmultinom(
  n = I, size = 1, 
  prob = rep(1/G,G)
  ))
  cluster_labels <- max.col(dummy_labels, ties.method = "first")
  # Draw data from multivariate normal distriubution
  X <- t(sapply(cluster_labels, function(lbl){
    mvrnorm(1, mu = curves[, lbl], Sigma = IJ)
  }
  ))
  # Generate the error
  E <- matrix(rnorm(I * J, sd = var.err), nrow = I, ncol = J)
  # Compute the data matrix X 
  X <- X + E
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
      randomstarts = randomstarts_cv
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
  simulation_results$FARI[iter] <- ARI.F(cluster_labels,res$U)
  ABp <- res$A %*% t(res$B)
  tcurves <- t(curves)
  # Permute the estimated curves to match the true curves
  perm <- perm_hungarian_fast(tcurves, ABp, J)
  simulation_results$SSQerr[iter] <- sum((ABp[perm,] - t(curves))^2)
  cat("End Monte Carlo Simulation: ", iter, " Lambda* ", simulation_results$lambda_best[iter], " Gamma* ", simulation_results$gamma_best[iter],
    "FARI: ", simulation_results$FARI[iter], " sSqErr: ", simulation_results$SSQerr[iter], "\n")
}
# Results
mean(simulation_results$FARI)
sd(simulation_results$FARI)
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
