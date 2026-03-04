# This is a test file for the equality of the Fuzzy Entropic Functional Reduced K Means
i <- 20
j <- 15
q <- 4
g <- 7

A <- matrix(rnorm(g*q), nrow = g, ncol = q)
B <- matrix(rnorm(j*q), nrow = j, ncol = q)
C <- matrix(rnorm(g*j), nrow = g, ncol = j)
K <- matrix(rchisq(j*j, 1), nrow = j, ncol = j)
D <- diag(rchisq(g, 1))

# smoothing penalty in compact matrix form
value <- 0
for(h in c(1:g)){
    value <- value + t(A[h,])%*%t(B)%*%K%*%B%*%A[h,]
}
print(value)
print(
    sum(
        diag(
            t(B)%*%K%*%B%*%t(A)%*%A
        )
    )
)

# mean and centroids in compact matrix form
value <- 0
for(h in c(1:g)){
    value <- value + D[h,h]^2 * sum((matrix(C[h,],ncol = 1) - B%*%A[h,])^2)
}
print(value)
print(
    sum(
        diag(
            t(D%*%C - D%*%A%*%t(B))%*%(D%*%C - D%*%A%*%t(B))
        )
    )
)
