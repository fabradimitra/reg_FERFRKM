init_FERFRKM <- function(X, G, Q, seed = NULL, nstart_kmeans = 10) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  #
  kmeans_res <- kmeans(X, G, nstart = nstart_kmeans)
  U_init <- matrix(0, nrow = nrow(X), ncol = G)
  U_init[cbind(seq_len(nrow(X)), kmeans_res$cluster)] <- 1
  #
  svd_init <- svd(diag(1 / colSums(U_init)) %*% t(U_init) %*% X)
  A_init <- svd_init$u[, seq_len(Q), drop = FALSE] # drop = FALSE forces R to keep a matrix
  B_init <- svd_init$v[, seq_len(Q), drop = FALSE] %*% diag(svd_init$d[seq_len(Q)])
  #
  list(U = U_init, A = A_init, B = B_init)
}