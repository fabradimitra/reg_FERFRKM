randgenuc <- function(n,k){
  M <- matrix(0L, nrow = n, ncol = k)
  M[cbind(seq_len(n), sample(c(1:k), n, replace = TRUE))] <- 1L
  M
}
