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
load("data/modelsel_plane_XB.RData")
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
pred_class <- max.col(res$U[,c(2,4,6,1,3,5)], ties.method = "first")
conf_mat_raw <- table(True = tcm, Predicted = pred_class)

plot_confusion_matrix <- function(mat, main) {
  mat <- as.matrix(mat)
  nr <- nrow(mat)
  nc <- ncol(mat)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 5, 4, 2) + 0.1)
  image(
    x = seq_len(nc),
    y = seq_len(nr),
    z = t(mat[nr:1, , drop = FALSE]),
    col = heat.colors(30),
    axes = FALSE,
    xlab = "Predicted",
    ylab = "True",
    main = main
  )
  axis(1, at = seq_len(nc), labels = colnames(mat))
  axis(2, at = seq_len(nr), labels = rev(rownames(mat)), las = 1)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      text(j, nr - i + 1, labels = mat[i, j])
    }
  }
}
plot_confusion_matrix(conf_mat_raw, "Confusion matrix")
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

# Plot curves of Mirage and Eurofighter classes:
Mirage <- X[which(tcm == 1),]
Ymirage <- apply(Mirage, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ymirage, type = "l", lwd = 1, lty = 1,
  col = c("red"),
  xlab = "", ylab = "", ylim = c(-1.5,2),
  main = "Mirages & Eurofighters"
)
Eurofighters <- X[which(tcm == 2),]
Yeurofighter <- apply(Eurofighters, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Yeurofighter, lwd = 1, lty = 1,
  col = c("blue")
)
est_centroid_2 <- t(as.matrix(est_centroids[2,]))
Ymr <- apply(est_centroid_2, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 3, lty = 2,
  col = c("black")
)
# Centroids and estimated together: 
matplot(
  tt, Ymr, type = "l", lwd = 1, lty = 1,
  col = c("black"),
  xlab = "", ylab = "", ylim = c(-1.5,2),
  main = "Mirages & Eurofighters"
)
centroid_mir <- apply(centroid_1, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, centroid_euro, lwd = 1, lty = 1,
  col = c("red")
)
centroid_euro <- apply(centroid_2, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, centroid_euro, lwd = 1, lty = 1,
  col = c("blue")
)

# Centroid 3 ----
Ym <- apply(centroid_3, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 2, lty = 1,
  col = c("gold"),
  xlab = "", ylab = "", ylim = c(-1.5,2),
  main = "F14 wing closed"
)
est_centroid_3 <- t(as.matrix(est_centroids[6,]))
Ymr <- apply(est_centroid_3, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 2, lty = 2,
  col = c("black")
)
# Centroid 4 ----
Ym <- apply(centroid_4, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Ym, type = "l", lwd = 2, lty = 1,
  col = c("darkgreen"),
  xlab = "", ylab = "", ylim = c(-1.5,2),
  main = "F14 wing open"
)
est_centroid_4 <- t(as.matrix(est_centroids[1,]))
Ymr <- apply(est_centroid_4, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 2, lty = 2,
  col = c("black")
)
# Centroid 5 ----
# plot harrier curves:
harriers_1 <- X[which((tcm == 5)&(pred_class==5)),]
ycurvesharriers_1 <- apply(harriers_1, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, ycurvesharriers_1, type = "l", lwd = 1, lty = 1,
  col = c("blue"),
  xlab = "", ylab = "", ylim = c(-2,3),
  main = "Harrier 1"
)
harriers_2 <- X[which((tcm == 5)&!(pred_class==5)),]
ycurvesharriers_2 <- apply(harriers_2, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, ycurvesharriers_2, lwd = 1, lty = 1,
  col = c("red")
)

# Centroids Harrier 1
centroid_harrier_1 <- t(as.matrix(colMeans(X[which((tcm == 5)&(pred_class==5)),])))
Yharriers_1 <- apply(centroid_harrier_1, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Yharriers_1, type = "l", lwd = 2, lty = 1,
  col = c("blue"),
  xlab = "", ylab = "",
  main = "Harrier 1"
)
est_centroid_harrier_1 <- t(as.matrix(est_centroids[3,]))
Ymr <- apply(est_centroid_harrier_1, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 2, lty = 2,
  col = c("black")
)

# Centroids Harrier 2
centroid_harrier_2<- t(as.matrix(colMeans(X[which((tcm == 5)&!(pred_class==5)),])))
Yharriers_2 <- apply(centroid_harrier_2, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, Yharriers_2, type = "l", lwd = 2, lty = 1,
  col = c("red"),
  xlab = "", ylab = "",
  main = "Harrier 2"
)
est_centroid_harrier_2 <- t(as.matrix(est_centroids[5,]))
Ymr <- apply(est_centroid_harrier_2, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, Ymr, lwd = 2, lty = 2,
  col = c("black")
)

# Plot curves of F22 and F15 classes:
F22 <- X[which(tcm == 6),]
YF22 <- apply(F22, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matplot(
  tt, YF22, type = "l", lwd = 1, lty = 1,
  col = c("red"),
  xlab = "", ylab = "", ylim = c(-1.5,2),
  main = "F-22 & F-15"
)
F15 <- X[which(tcm == 7),]
YF15 <- apply(F15, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, YF15, lwd = 1, lty = 1,
  col = c("blue")
)

est_centroid_F22_15 <- t(as.matrix(est_centroids[4,]))
Ymr <- apply(est_centroid_F22_15, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
# Centroids and estimated together: 
matplot(
  tt, Ymr, type = "l", lwd = 1, lty = 1,
  col = c("black"),
  xlab = "", ylab = "", ylim = c(-2,3),
  main = "F-22 & F-15"
)
centroid_F22 <- apply(centroid_6, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, centroid_F22, lwd = 1, lty = 1,
  col = c("red")
)
centroid_F15 <- apply(centroid_7, 1, function(y) splinefun(t_grid, y, method = "natural")(tt))
matlines(
  tt, centroid_F15, lwd = 1, lty = 1,
  col = c("blue")
)

