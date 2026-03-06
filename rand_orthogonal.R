rand_orthogonal <- function(G,Q) {
  stopifnot(Q <= G)
  Z <- matrix(rnorm(G * Q), nrow = G, ncol = Q)
  qrZ <- qr(Z)
  A <- qr.Q(qrZ)     
  return(A)
}