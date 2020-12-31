rm(list=ls())
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho.R")

paramMat <- cbind(500, 100, seq(0.1, 0.15, length.out = 11), 0.45, 0.1, 0.45)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

#######################################

# against sum

vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1.5, 0.4, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, 0.4, -0.4)) %*% t(eigen_mat)
round(B1,2); round(B2,2); colSums(B1+B2)
K <- 3

vec <- paramMat[1,]
n <- vec["n"]; L <- vec["L"]; # rho <- vec["rho"]
rho <- 1
mem_prop1 <- 1/3; mem_prop2 <- 1/3; mem_prop3 <- 1/3

membership_vec <- c(rep(1, floor(mem_prop1*n)), rep(2, floor(mem_prop2*n)), rep(3, floor(mem_prop3*n)))
if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))

prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})

set.seed(10)
adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss_debias")
res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "sum")
res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

res1; res2

#######################################

vec1 <- c(1,2,2)/3
vec2 <- c(-2,-1,2)/3
vec3 <- c(2,-2,1)/3
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1, 0.5, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, 0.1, 0.6)) %*% t(eigen_mat)
round(B1,2); round(B2,2); colSums(B1+B2)
K <- 3

vec <- paramMat[1,]
n <- vec["n"]; L <- vec["L"]; # rho <- vec["rho"]
rho <- 0.1 # can work even at rho <- 0.04, but not further
# mem_prop1 <- 0.45; mem_prop2 <- 0.1; mem_prop3 <- 0.45
mem_prop1 <- 0.4; mem_prop2 <- 0.1; mem_prop3 <- 0.5

membership_vec <- c(rep(1, floor(mem_prop1*n)), rep(2, floor(mem_prop2*n)), rep(3, floor(mem_prop3*n)))
if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))

prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})

set.seed(10)
adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss_debias")
res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss")
res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

res1; res2

probsq_array <- .form_array(prob_list)
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
basis2_orth <- eigen(yy + diag(diag_vec))$vectors[,-c(1:K)]
dist_vec <- svd(t(basis1) %*% basis2_orth)$d 
dist_vec # this is a useless metric it seems. the more important metric is the eigengap, 
# i.e., eigen(yy + diag(diag_vec))$values[3:4]

par(mfrow = c(1,2)); plot(diag(yy)); plot(diag(yy)+diag_vec)

par(mfrow = c(1,3))
tmp <- eigen(yy + diag(diag_vec))$values; tmp[3:12]; tmp[3]/tmp[4]
image(t(basis1)); image(t(basis2)); plot(tmp[3:12], ylim = c(0, tmp[3]))

zz <- eigen(yy + diag(diag_vec))
tmp2 <- zz$vectors[,-c(1:2)] %*% diag(zz$values[-c(1:2)]) %*% t(zz$vectors[,-c(1:2)])
par(mfrow = c(1,1)); image(tmp2, asp = T)
quantile(tmp2)

###########################################

# setting to also ruin reweighting
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)
K <- 3

vec <- paramMat[1,]
n <- vec["n"]; L <- vec["L"]; # rho <- vec["rho"]
rho <- 0.15 # let's go from rho = 0.12 to 0.18 (roughly)
# mem_prop1 <- 0.45; mem_prop2 <- 0.1; mem_prop3 <- 0.45
mem_prop1 <- 0.4; mem_prop2 <- 0.1; mem_prop3 <- 0.5

membership_vec <- c(rep(1, floor(mem_prop1*n)), rep(2, floor(mem_prop2*n)), rep(3, floor(mem_prop3*n)))
if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))

prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})

set.seed(10)
adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss_debias")
res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss")
res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)

res1; res2

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "sum")
res <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)

res

set.seed(10)
agg_network <- networkSoSD::aggregate_networks(adj_list, method = "ss")
res <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)

res

probsq_array <- .form_array(prob_list)
for(i in 1:dim(probsq_array)[1]){
  probsq_array[i,,] <- crossprod(probsq_array[i,,])
}
yy <- apply(probsq_array, c(2,3), sum)

diag_vec <- rep(0, dim(probsq_array)[2])
for(i in 1:dim(probsq_array)[1]){
  diag_vec <- diag_vec + colSums(prob_list[[i]])
}

par(mfrow = c(1,3))
tmp <- eigen(yy + diag(diag_vec))$values; tmp[3:12]; tmp[3]/tmp[4]
image(t(basis1)); image(t(basis2)); plot(tmp[3:12], ylim = c(0, tmp[3]))

