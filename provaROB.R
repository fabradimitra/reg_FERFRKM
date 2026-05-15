source("kspline.R")
source("randgenuf.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FERFRKM_BcolwiseNconstr.R")
source("fungcv.R")
source("perm_hungarian_fast.R")
I<-150; J<-50; G <- 4; Q <- 2
#
set.seed(123)
U_true <- randgenuc(I, G)
M <- matrix(rnorm(G*J),G,J)
X <- U_true%*%M
U_init <- randgenuc(I, G)
A_init <- rand_orthogonal(G, Q)
# A_init <- matrix(rnorm(G * Q), nrow = G, ncol = Q)
B_init <- t(t(A_init)%*%solve(t(U_init)%*%U_init)%*%t(U_init)%*%X)
res <- kspline(1:J)
K <- res$K
Pk_f <- res$Pk
Lk_f <- res$Lk
idx <- diag(Lk_f > 1e-6)
Pk <- Pk_f[,idx] # For reconstruction of K
Lk <- Lk_f[idx,idx] # For reconstruction of K
lambda <- 0.01; gamma <- 1; max_iter <- 1000; tol <- 1e-6
res <- FERFRKM_BcolwiseNconstr(C=X,
                               K=K,
                               Pk=Pk,
                               Lk=Lk,
                               U=U_init,
                               A=A_init,
                               B=B_init,
                               lambda=lambda,
                               gamma = gamma,
                               max_iter = max_iter,
                               tol = tol)
