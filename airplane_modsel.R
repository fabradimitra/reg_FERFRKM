require(fclust)
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
source("preggq_int.R")
# load dataset
load("data/plane.RData")
X <- scale(X,center=TRUE,scale=FALSE)
I <- nrow(X)
J <- ncol(X)
# Hyperparams
randomstarts <- 5
randomstarts_cv <- 3
kmeans_starts <- 20
folds <- 5
seed <- 123
lambda_init <- 1 
gamma_init <- 1
#
t_grid <- seq(1, ncol(X))
# Cluster centroids
res <- kspline(t_grid)
K <- res$K
Pk <- res$Pk
Lk <- res$Lk
# G and Q selection
G_max <- 9
gridGQ <- do.call(rbind, lapply(2:G_max, function(G) {
  data.frame(
    G = G,
    Q = 2:(G - 1),
    lambda = 0,
    gamma = 0,
    wdev = 0,
    ARI = 0,
    loss = 0
  )
}))
for(i in nrow(gridGQ)){
  G <- gridGQ[i,1]
  Q <-  gridGQ[i,2]
  # Cross validation
  invisible(capture.output(
    cv_res <- CV_FERFRKM(
      Xtr = X,
      G = G,
      Q = Q,
      K = K,
      Pk = Pk,
      Lk = Lk,
      lambda_init = lambda_init,
      gamma_init = gamma_init,
      folds = folds,
      max_iter = Inf,
      tol = 1e-8,
      nstart_kmeans = kmeans_starts,
      seed = seed,
      randomstarts = randomstarts_cv
    ),
    type = "output"
  ))
  # Fit the best combination:
  cur_loss <- Inf
  for(start in seq_len(randomstarts)){
    if(start == 1){
      init <- init_FERFRKM(X, G, Q, seed = seed, nstart_kmeans = kmeans_starts) 
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
  gridGQ[i,3:7] <- c(cv_res$par[1],cv_res$par[2],res$wdev, ARI.F(tcm,res$U), cur_loss) 
  cat("G: ", G, "Q: ", Q, " Lambda* ", gridGQ[i,3], " Gamma* ", gridGQ[i,4],
    "within-cluster deviance: ", gridGQ[i,5], " FARI: ", gridGQ[i,6] ," loss: ", gridGQ[i,7],"\n")
  save.image("data/fit_plane.RData")
}
load("data/fit_plane.RData")
# Computing the Xie-Beni index
compute_xie_beni <- function(X, centers, gamma) {
  cnorm2 <- rowSums(X^2)
  vnorm2 <- rowSums(centers^2)
  dist2 <- outer(cnorm2, vnorm2, "+") - 2 * (X %*% t(centers))
  dist2[!is.finite(dist2)] <- Inf
  U <- exp(-dist2 / gamma)
  row_sums <- rowSums(U)
  row_sums[!is.finite(row_sums) | row_sums == 0] <- 1
  U <- U / row_sums
  U[!is.finite(U)] <- 0
  U[U < 1e-12] <- 1e-12
  sep2 <- as.matrix(dist(centers))^2
  sep2[sep2 == 0] <- Inf
  min_sep2 <- min(sep2)
  if (!is.finite(min_sep2) || min_sep2 <= 0) {
    return(Inf)
  }
  score <- sum(U * dist2) / (nrow(X) * min_sep2)
  if (!is.finite(score)) {
    Inf
  } else {
    score
  }
}

modelsel$XB <- NA_real_
for (i in seq_len(nrow(modelsel))) {
  G <- modelsel$G[i]
  Q <- modelsel$Q[i]
  cur_loss <- Inf
  res <- NULL
  for (start in seq_len(randomstarts)) {
    if (start == 1) {
      init <- init_FERFRKM(X, G, Q, seed = seed, nstart_kmeans = kmeans_starts)
    } else {
      U_init <- randgenuc(I, G)
      A_init <- rand_orthogonal(G, Q)
      B_init <- t(t(A_init) %*% solve(t(U_init) %*% U_init) %*% t(U_init) %*% X)
      init <- list(U = U_init, A = A_init, B = B_init)
    }
    res_cur <- tryCatch(
      FERFRKM(
        C = X,
        K = K,
        Pk = Pk,
        Lk = Lk,
        U = init$U,
        A = init$A,
        B = init$B,
        lambda = modelsel$lambda[i],
        gamma = modelsel$gamma[i],
        max_iter = Inf,
        tol = 1e-8
      ),
      error = function(e) NULL
    )
    if (is.null(res_cur)) {
      next
    }
    if (cur_loss > res_cur$loss_function) {
      res <- res_cur
      cur_loss <- res$loss_function
    }
  }
  if (!is.null(res)) {
    modelsel$XB[i] <- compute_xie_beni(X, res$A %*% t(res$B), modelsel$gamma[i])
  } else {
    modelsel$XB[i] <- Inf
  }
  cat("Row ", i, " G: ", G, " Q: ", Q, " XB: ", modelsel$XB[i], "\n")
}
save(modelsel, file = "data/modelsel_plane_XB.RData")