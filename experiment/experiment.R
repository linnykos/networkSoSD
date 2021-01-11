rm(list=ls())
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()

.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,2,2)/3
vec2 <- c(-2,-1,2)/3
vec3 <- c(2,-2,1)/3
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(0.8, 0.2, 1.5)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, -0.1, 0.6)) %*% t(eigen_mat)
round(B1, 2); round(B2, 2); colSums(B1); colSums(B2)
K <- 3

val3 <- 0.4; val1 <- (1-val3)/2

paramMat <- cbind(500, 100, seq(0.025, 0.2, length.out = 15), val1, val1, val3)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

rule <- function(vec){
  n <- vec["n"]; L <- vec["L"]; rho <- vec["rho"]
  mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list, prob_list = prob_list)
}

set.seed(10)
dat <- rule(paramMat[10,]); zz = colSums(dat$prob_list[[1]]); plot(zz)
# 
# agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
# res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
# 
# agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
# res3 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
# 
# 
# res1
# res3
# 
