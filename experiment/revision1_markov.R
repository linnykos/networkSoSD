rm(list=ls())
# ?quadprog::solve.QP

n <- 6
delta <- 0.2

D <- matrix(0, n, n)
diag(D) <- rep(c(0.5^2, 0.4^2, 0.1^2), each = 2)
D[1,6] <- D[6,1] <- 0.1*0.5
D[2,4] <- D[4,2] <- 0.5*0.4
D[3,5] <- D[5,3] <- 0.4*0.1
D <- 2*D

tmp <- c(0.5*0.4, 0.5*0.1, 0.5*0.4, 0.4*0.1, 0.5*0.1, 0.4*0.1)
d <- 2*tmp - 2*tmp*(1-delta)

A <- matrix(0, nrow = 3, ncol = n)
A[1,1:2] <- 1
A[2,3:4] <- 1
A[3,5:6] <- 1
A <- t(A)

b0 <- rep(delta, 3)

tmp <- eigen(D)
D2 <- tmp$vectors[,1:3] %*% diag(tmp$values[1:3]) %*% t(tmp$vectors[,1:3])

# res <- quadprog::solve.QP(Dmat = D2, dvec = d,
#                 Amat = A, bvec = b0, meq = 3)

res <- kernlab::ipop(c = matrix(-d, ncol = 1), 
                     H = D, A = t(A), 
                     b = matrix(b0, ncol = 1), 
                     l = matrix(rep(0,6), ncol = 1), 
                     u = matrix(rep(1,6), ncol = 1), 
                     r = matrix(rep(0, 3), ncol = 1))

#################
P <- matrix(0, 3, 3)
diag(P) <- 1-delta
P[1,2:3] <- res@primal[1:2]
P[2,c(1,3)] <- res@primal[3:4]
P[3,1:2] <- res@primal[5:6]
x <- c(0.5, 0.4, 0.1)
round(P, 3)

x %*% P
x %*% round(P, 4)






