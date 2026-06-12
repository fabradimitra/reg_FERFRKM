require(clue)
source("perm_hungarian_fast.R")
source("kspline.R")
source("randgenuc.R")
source("loss_function.R")
source("FERFRKM.R")
source("CV_FERFRKM.R")
source("init_FERFRKM.R")
source("make_folds.R")
source("rand_orthogonal.R")
source("preggq_int.R")
# Hyper-parameters
randomstarts <- 5
seed <- 123
kmeans_starts <- 20

# load data
load("data/Plane.RData")
X <- scale(X, center = TRUE, scale = FALSE)
res_kspline <- kspline(seq_len(ncol(X)))
K <- res_kspline$K
Pk <- res_kspline$Pk
Lk <- res_kspline$Lk
I <- nrow(X)
J <- ncol(X)

# load model selection
load("data/modelsel_plane.RData")
modelsel[] <- lapply(modelsel, as.numeric)
modelsel$loss[2] <- max(modelsel$loss[-2])

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

modelsel$XB[2] <- max(modelsel$XB[-2])
best_idx <- which.min(modelsel$XB)
G <- modelsel$G[best_idx]
Q <- modelsel$Q[best_idx]
# Refit best model
cur_loss <- Inf
for (start in seq_len(randomstarts)) {
  if (start == 1) {
    init <- init_FERFRKM(X, G, Q, seed = seed, nstart_kmeans = kmeans_starts)
  } else {
    U_init <- randgenuc(I, G)
    A_init <- rand_orthogonal(G, Q)
    B_init <- t(t(A_init) %*% solve(t(U_init) %*% U_init) %*% t(U_init) %*% X)
    init <- list(U = U_init, A = A_init, B = B_init)
  }
  # Run FERFRKM algorithm
  res_cur <- tryCatch(
    FERFRKM(
      C = X,
      K = K,
      Pk = Pk,
      Lk = Lk,
      U = init$U,
      A = init$A,
      B = init$B,
      lambda = modelsel$lambda[best_idx],
      gamma = modelsel$gamma[best_idx],
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
est_centroids <- res$A %*% t(res$B)
centroid_1 <- t(as.matrix(colMeans(X[which(tcm == 1),])))
centroid_2 <- t(as.matrix(colMeans(X[which(tcm == 2),])))
centroid_3 <- t(as.matrix(colMeans(X[which(tcm == 3),])))
centroid_4 <- t(as.matrix(colMeans(X[which(tcm == 4),])))
centroid_5 <- t(as.matrix(colMeans(X[which(tcm == 5),])))
centroid_6 <- t(as.matrix(colMeans(X[which(tcm == 6),])))
centroid_7 <- t(as.matrix(colMeans(X[which(tcm == 7),])))
true_centroids <- rbind(centroid_1,centroid_2,centroid_3,centroid_4,centroid_5,centroid_6,centroid_7)
# Plot the centroids and their reconstruction for one iteration ----
tt <- seq(1, 144, length.out = 400)
t_grid <- c(1:J)
# Centroid 1 ----
Ym <- apply(centroid_1, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
est_centroid_1 <- t(as.matrix(est_centroids[2,]))
Ymr <- apply(est_centroid_1, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
# Centroid 2 ----
Ym <- apply(centroid_2, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
est_centroid_2 <- t(as.matrix(est_centroids[6,]))
Ymr <- apply(est_centroid_2, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
# Centroid 3 ----
Ym <- apply(centroid_3, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
est_centroid_3 <- t(as.matrix(est_centroids[4,]))
Ymr <- apply(est_centroid_3, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
# Centroid 4 ----
Ym <- apply(centroid_4, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
est_centroid_4 <- t(as.matrix(est_centroids[5,]))
Ymr <- apply(est_centroid_4, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
# Centroid 5 ----
Ym <- apply(centroid_5, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
est_centroid_5 <- t(as.matrix(est_centroids[1,]))
Ymr <- apply(est_centroid_5, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
# Centroid 6 ----
Ym <- apply(centroid_6, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
est_centroid_6 <- t(as.matrix(est_centroids[6,]))
Ymr <- apply(est_centroid_6, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
# Centroid 7 ----
Ym <- apply(centroid_7, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 1, lty = 1,
  col = c("red","blue","darkgreen","orange"),
  xlab = "", ylab = ""
)
est_centroid_7 <- t(as.matrix(est_centroids[7,]))
Ymr <- apply(est_centroid_7, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 1, lty = 2,
  col = c("red","blue","darkgreen","orange")
)
