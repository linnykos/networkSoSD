rm(list=ls())
library(simulation)
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho.R")

paramMat <- cbind(500, 100, seq(0.1, 0.15, length.out = 11), 0.45, 0.1, 0.45)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

.l2norm <- function(x){sqrt(sum(x^2))}
# vec1 <- c(1,1,sqrt(2))
# vec2 <- c(1,1,-sqrt(2))
# vec3 <- c(-1,1,0)
vec1 <- c(1,2,2)/3
vec2 <- c(-2,-1,2)/3
vec3 <- c(2,-2,1)/3
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1, 0.5, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, -0.1, 0.6)) %*% t(eigen_mat)
# B2 <- eigen_mat %*% diag(c(1.5, 0.1, 0.6)) %*% t(eigen_mat)
# it seems like we don't need positive/negative eigenvalues to really show debiasing is needed
round(B1,2); round(B2,2); colSums(B1+B2)
K <- 3

vec <- paramMat[1,]
n <- vec["n"]; L <- vec["L"]; # rho <- vec["rho"]
# rho <- 0.05
rho <- 1
mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]

membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))

prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})

#########################################

set.seed(10)
adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})

agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss_debias")
res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss")
res3 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

res1; res3

#######################################
# now compute the supposed biases
# first convert to array
probsq_array <- .form_array(prob_list)
# now square each layer
for(i in 1:dim(probsq_array)[1]){
  probsq_array[i,,] <- crossprod(probsq_array[i,,])
}
yy <- apply(probsq_array, c(2,3), sum)

diag_vec <- rep(0, dim(probsq_array)[2])
for(i in 1:dim(probsq_array)[1]){
  diag_vec <- diag_vec + colSums(prob_list[[i]])
}

basis1 <- eigen(yy)$vectors[,1:K]
basis2 <- eigen(yy + diag(diag_vec))$vectors[,1:K]

par(mfrow = c(1,2))
image(basis1); plot(diag(yy))
image(basis2); plot(diag(yy + diag(diag_vec)))
basis2_orth <- eigen(yy + diag(diag_vec))$vectors[,-c(1:3)]
dist_vec <- svd(t(basis1) %*% basis2_orth)$d 
dist_vec

########################################

# now compute the supposed biases
# first convert to array
probsq_array <- .form_array(prob_list)
# now square each layer
for(i in 1:dim(probsq_array)[1]){
  probsq_array[i,,] <- crossprod(probsq_array[i,,])
}
# now take the eigen-decomp of the sum
yy <- apply(probsq_array, c(2,3), sum)
yy2 <- yy
max_val <- max(probsq_array)
for(i in 1:n){
  yy2[i,i] <- yy2[i,i] + L*n*max_val
}
yy3 <- yy
diag_vec <- rep(0, dim(probsq_array)[2])
for(i in 1:dim(probsq_array)[1]){
  diag_vec <- diag_vec + colSums(prob_list[[i]])
}
for(i in 1:n){
  yy3[i,i] <- yy3[i,i] + L*n*max_val + diag_vec[i]
}

image(yy)

basis1 <- eigen(yy2)$vectors[,1:3]
basis2 <- eigen(yy3)$vectors[,1:3]
basis2_orth <- eigen(yy3)$vectors[,-c(1:3)]
dist_vec <- svd(t(basis1) %*% basis2_orth)$d 
dist_vec

##########################################

# previous attempt

# now compute the supposed biases
# first convert to array
probsq_array <- .form_array(prob_list)
# now square each layer
for(i in 1:dim(probsq_array)[1]){
  probsq_array[i,,] <- crossprod(probsq_array[i,,])
}
# now take the eigen-decomp of the sum
yy <- apply(probsq_array, c(2,3), sum)
Matrix::rankMatrix(yy)
basis1 <- eigen(yy)$vectors[,1:3]
head(eigen(yy)$values)
image(t(basis1)) # this is our target population quantity

# now take the eigen-decomp of the debiased
diag_vec <- rep(0, dim(probsq_array)[2])
for(i in 1:dim(probsq_array)[1]){
  diag_vec <- diag_vec + colSums(prob_list[[i]])
}
zz <- yy + diag(diag_vec) # this is what the "denoised" version of the variant w/ debiasing looks like
head(eigen(zz)$values)
basis2 <- eigen(zz)$vectors[,1:3]
plot(diag_vec)
plot(colSums(prob_list[[1]]))
plot(colSums(prob_list[[100]]))
plot(basis2[,2])
image(t(basis2))
Matrix::rankMatrix(zz)
plot(eigen(zz)$values[1:10], ylim = c(0, max(eigen(zz)$values)))

basis2_orth <- eigen(zz)$vectors[,-c(1:3)]
dist_vec <- svd(t(basis1) %*% basis2_orth)$d 
dist_vec


##########################################

# below are just random stuff, looking at the minimum distances between clusters

.min_distance <- function(mat){
  idx <- sort(unique(unlist(lapply(1:ncol(mat), function(j){which(abs(diff(mat[,j])) >= 1e-6)}))))
  idx <- c(1,idx+1)
  mat2 <- mat[idx,]
  zz <- as.matrix(dist(mat2))
  for(j in 1:ncol(mat)){
    zz[j,] <- zz[j,]/.l2norm(mat2[,j])
  }
  diag(zz) <- Inf
  min(zz)
}

.min_distance(basis1)
.min_distance(basis2)

