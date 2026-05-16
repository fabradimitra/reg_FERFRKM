source("kspline.R")
source("randgenuf.R")
source("randgenuc.R")
source("rand_orthogonal.R")
source("loss_function.R")
source("FERFRKM.R")
source("plotX.R")
source("fungcv.R")
source("perm_hungarian_fast.R")
I<-100; J<-12; G <- 4; Q <- 2
res <- kspline(1:J)
K <- res$K
Pk_f <- res$Pk
Lk_f <- res$Lk
idx <- diag(Lk_f > 1e-6)
Pk <- Pk_f[,idx] # For reconstruction of K
Lk <- Lk_f[idx,idx] # For reconstruction of K
#
set.seed(123)
U_true <- randgenuc(I, G)
A_true <- matrix(rnorm(G*Q),G,Q)
B <- matrix(rnorm(J*Q),J,Q)
plotX(t(B),1:2)
B_true <- solve(diag(J)+3*K)%*%B
plotX(t(B_true),1:2)
M <- A_true%*%t(B_true)
X <- U_true%*%M + 0.3*matrix(rnorm(I*J),I,J)
U_init <- randgenuc(I, G)
A_init <- rand_orthogonal(G, Q)
# A_init <- matrix(rnorm(G * Q), nrow = G, ncol = Q)
B_init <- t(t(A_init)%*%solve(t(U_init)%*%U_init)%*%t(U_init)%*%X)

lambda <- 2; gamma <- 1; max_iter <- Inf; tol <- 1e-10
res <- FERFRKM(C=X,K=K,Pk=Pk,Lk=Lk,U=U_init,A=A_init,B=B_init,
               lambda=lambda,gamma = gamma,max_iter = max_iter,tol = tol)
t(U_true)%*%res$U

plotX(rbind(M,res$A%*%t(res$B)),c(1,1,1,1,2,2,2,2))

