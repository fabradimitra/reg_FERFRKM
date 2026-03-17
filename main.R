source("kspline.R")
source("randgenu.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FEFRKM.R")
# Simulated data example
set.seed(5)
I <-100; J <- 30; Q <- 2; G <- 3
X <- matrix(rnorm(I*J),I,J)
K <- kspline(1:J)
U <- randgenu(I,G)
A <- qr.Q(qr(matrix(rnorm(G*Q),G,Q)))
B <- t(t(A)%*%solve(t(U)%*%U)%*%t(U)%*%X)
#
FEFRKM(X,K,U,A,B,lambda=10,gamma=1,max_iter = Inf,tol = 1e-6)
