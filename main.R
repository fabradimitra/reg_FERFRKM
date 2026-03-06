library(fda)
library(splines2)
source("kspline.R")
source("randgenu.R")
source("rand_orthogonal.R")

# Loading data growth
data(growth)
with(growth, matplot(age, hgtf[, 1:10], type="b"))

# Natural cubic spline basis functions
x_obs  <- growth$age
iknots <- growth$age[-c(1, length(growth$age))]
nc_obs  <- naturalSpline(
        x = x_obs, 
        knots = iknots, 
        derivs = 0)  

# Focus on the height of the 39 male individuals
Y      <- t(growth$hgtm)                                     

# Matrix of coefficients for the natural spline basis functions for each of the 39 individuals
C <- solve(crossprod(nc_obs))%*%crossprod(nc_obs, t(Y)) 

# Plotting the natural spline basis functions and the fitted curves
x_grid <- seq.int(1, 18, .05)
nc_grid <- naturalSpline(
        x = x_grid, 
        knots = iknots,
        Boundary.knots = range(x_obs))   # 73 x 32
matplot(
  x_grid, nc_grid,
  type = "l", lty = 1, lwd = 2,
  xlab = "Age", ylab = "Basis value",
  main = "Natural Spline Basis Functions")

yhat   <- nc_grid %*% C

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
J <- ncol(C)

# Building matrix K
K <- kspline(x_obs)
# Trasposing matrix C
C <- t(C)
# Random generate initial U matrix
U_init <- randgenu(I, G)
# Random generate initial A matrix
A <- rand_orthogonal(G, Q)
# Check A'A=I_Q constraint
t(A) %*% A 
# Random generate initial B matrix
B <- matrix(rnorm(G * Q), nrow = G, ncol = Q)