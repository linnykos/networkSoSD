##
# this generates B1, B2 with equal colSum, but then no connection between (12) and (3)
.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(0.6, 0.6, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(0.6, 0.6, -0.4)) %*% t(eigen_mat)
round(B1, 2); round(B2, 2); colSums(B1); colSums(B2)
K <- 3

#############

# this doesn't seem to make a difference b/w ss-debias and ss ... (using rho=0.015)
.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1, 0.3, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1, 0.3, -0.4)) %*% t(eigen_mat)
round(B1, 2); round(B2, 2); colSums(B1); colSums(B2)
K <- 3
zz <- solve(B1, rep(1,3)); zz <- zz/sum(zz)


#######################

# even though degree are the same (also matrices), for paramMat[1,], ss-debias
#  still is better than ss
.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,2,2)/3
vec2 <- c(-2,-1,2)/3
vec3 <- c(2,-2,1)/3
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1, 0.3, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1, 0.3, 0.4)) %*% t(eigen_mat)
round(B1, 2); round(B2, 2); colSums(B1); colSums(B2)
K <- 3

zz <- solve(B1, rep(1,3)); zz <- zz/sum(zz)
paramMat <- cbind(500, 100, seq(0.015, 0.2, length.out = 15), zz[1], zz[2], zz[3])
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")
