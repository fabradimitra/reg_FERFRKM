library(fda)
library(splines2)
data(growth)
with(growth, matplot(age, hgtf[, 1:10], type="b"))

# Natural cubic spline basis functions
x_obs  <- growth$age
iknots <- growth$age[-c(1, length(growth$age))]
ns_obs  <- ns(
        x = x_obs,
        knots = iknots,
        intercept = TRUE)
qrNS <- qr(ns_obs)
ns_orth <- qr.Q(qrNS) # Orthogonalized version of the natural spline basis functions
R <- qr.R(qrNS) # Upper triangular matrix from the QR decomposition

# Focus on the height of the 39 male individuals
# transforme into meters instead of centimeters
Y      <- t(growth$hgtm)/100           

# Matrix of coefficients for the natural spline basis functions for each of the 39 individuals
C <- solve(crossprod(ns_orth))%*%crossprod(ns_orth, t(Y)) 

# Plotting the natural spline basis functions and the fitted curves
x_grid <- seq.int(1, 18, .05)
ns_grid <- ns(
        x = x_grid, 
        knots = iknots,
        intercept = TRUE)
ns_orth_grid <- ns_grid %*% solve(R) 
matplot(
  x_grid, ns_orth_grid,
  type = "l", lty = 1, lwd = 2,
  xlab = "Age", ylab = "Basis value",
  main = "Natural Spline Basis Functions")
yhat   <- ns_orth_grid %*% C
matplot(x_obs, Y[1,], ty="b", lty=1, col = "grey",
        xlab = "Age",
        ylab = "Height (cm)")
matlines(x_grid, yhat[,1], col = "darkred")

# Hyperparameters of the FEFRKM algorithm
lambda <- 1
gamma  <- 1
G <- 3
Q <- 2
I <- nrow(Y)
J <- ncol(ns_obs)
# Building matrix K
K <- kspline(x_obs) # internal knots or full grid?
Korth <- solve(R) %*% K %*% t(solve(R)) # Orthogonalized version of K
# Trasposing matrix C
C <- t(C)
# Random generate initial U matrix
U_init <- randgenuf(I, G)
# Random generate initial A matrix
A_init <- rand_orthogonal(G, Q)
# Check A'A=I_Q constraint
t(A_init) %*% A_init 
# Random generate initial B matrix
B_init <- matrix(rnorm(J * Q), nrow = J, ncol = Q)
# Running FEFRKM algorithm
res <- FERFRKM(
  C = C,
  K = Korth,
  U = U_init,
  A = A_init,
  B = B_init,
  lambda = lambda,
  gamma = gamma,
  max_iter = 1000,
  tol = 1e-6
)