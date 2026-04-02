perm_hungarian_fast <- function(M_true, M_est, J) {
  # This function computes the optimal permutation of the estimated cluster centroids
  # to match the true cluster centroids, using the Hungarian algorithm.
  #######
  # Efficient function to compute the cost matrix
  # To plug in the Hungarian algorithm for label switching correction in clustering
  # row-wise squared norms
  norm_true <- rowSums(M_true^2)   # length G
  norm_est  <- rowSums(M_est^2)    # length G
  # Cross-product
  cross <- M_true %*% t(M_est)     # G x G matrix
  # Squared distances
  cost <- outer(norm_true, norm_est, "+") - 2 * cross
  # Convert to MSE (mean over rows)
  cost / J
  ########
  perm <- solve_LSAP(cost)
  as.vector(perm)
}