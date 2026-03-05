library(fda)
library(splines2)
data(growth)
with(growth, matplot(age, hgtf[, 1:10], type="b"))

argvals <- seq.int(0.75, 18.25, .25)
nsbase <- naturalSpline(x = argvals, knots = growth$age)
matplot(nsbase, ty="l")

Y<- t(growth$hgtm)

c    <- solve(crossprod(nsbase)) %*% crossprod(nsbase, t(Y))
yhat <- nsbase %*% c

matplot(argvals, Y[1,], ty="l", lty=1, col = "lightgrey", ylim = c(0.35, 0.77),
        xlab = "Location within corpus collosum",
        ylab = "Fractional Anisotropy")
matlines(argvals, yhat[,1], col = "darkred")