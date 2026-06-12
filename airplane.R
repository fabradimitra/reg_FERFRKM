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
# load model selection
load("data/modelsel_plane.RData")
modelsel[] <- lapply(modelsel, as.numeric)
modelsel$loss[2] <- max(modelsel$loss[-2])
# Select the number of dimensions of the subspace with objective elbow method.
preggq_int(modelsel$wdev[13:17], modelsel$G[13:17], modelsel$Q[13:17])
G <- 7
Q <- 4
I <- nrow(X)
J <- ncol(X)
# Refit best model
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
      lambda= modelsel$lambda[8],
      gamma = modelsel$gamma[8],
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
ABp <- res$A%*%t(res$B)
centroid_1 <- t(as.matrix(colMeans(X[which(tcm == 1),])))
centroid_2 <- t(as.matrix(colMeans(X[which(tcm == 2),])))
centroid_3 <- t(as.matrix(colMeans(X[which(tcm == 3),])))
centroid_4 <- t(as.matrix(colMeans(X[which(tcm == 4),])))
centroid_5 <- t(as.matrix(colMeans(X[which(tcm == 5),])))
centroid_6 <- t(as.matrix(colMeans(X[which(tcm == 6),])))
centroid_7 <- t(as.matrix(colMeans(X[which(tcm == 7),])))
true_centroids <- rbind(centroid_1,centroid_2,centroid_3,centroid_4,centroid_5,centroid_6,centroid_7)
perm <- perm_hungarian_fast(true_centroids, ABp, J)
est_centroids <- ABp[perm,]
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
est_centroid_2 <- t(as.matrix(est_centroids[2,]))
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
est_centroid_3 <- t(as.matrix(est_centroids[3,]))
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
est_centroid_4 <- t(as.matrix(est_centroids[4,]))
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
est_centroid_5 <- t(as.matrix(est_centroids[5,]))
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
