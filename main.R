library(fda)
library(splines2)
data(growth)
with(growth, matplot(age, hgtf[, 1:10], type="b"))

x_obs  <- growth$age
iknots <- growth$age[-c(1, length(growth$age))]
bs_obs  <- naturalSpline(
        x = x_obs, 
        knots = iknots, derivs = 0)   # 31 x 30

Y      <- t(growth$hgtm)                                     # 39 x 31

coef_basis <- solve(crossprod(bs_obs))%*%crossprod(bs_obs, t(Y))  # 32 x 39

x_grid <- seq.int(1, 18, .05)
bs_grid <- naturalSpline(
        x = x_grid, 
        knots = iknots,
        Boundary.knots = range(x_obs))   # 73 x 32

matplot(
  x_grid, bs_grid,
  type = "l", lty = 1, lwd = 2,
  xlab = "Age", ylab = "Basis value",
  main = "Natural Spline Basis Functions")

yhat   <- bs_grid %*% coef_basis

matplot(x_obs, Y[1,], ty="b", lty=1, col = "grey",
        xlab = "Age",
        ylab = "Height (cm)")
matlines(x_grid, yhat[,1], col = "darkred")

