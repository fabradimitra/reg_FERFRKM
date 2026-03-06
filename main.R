library(fda)
library(splines2)

# Loading data growth
data(growth)
with(growth, matplot(age, hgtf[, 1:10], type="b"))

# Natural spline basis functions
x_obs  <- growth$age
iknots <- growth$age[-c(1, length(growth$age))]
nc_obs  <- naturalSpline(
        x = x_obs, 
        knots = iknots, 
        derivs = 0)  

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

# Building matrix K

