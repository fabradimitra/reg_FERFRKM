preggq_int<-function(y,x1,x2){
  #
  # piecewise regression
  #
  nm <- length(y)
  res2 <- rep(0,nm)
  for(i in 1:nm){
    g <- x1[i]
    q <- x2[i]
    x1g <- (x1-g)*(x1>=g) 
    x2q <- (x2-q)*(x2>=q)
    op <- lm(y~x1+x2+x1g+x2q+x1*x2)
    res2[i] <- sum(op$residuals^2)
  }
  print(res2)
  print(which.min(res2))
}
