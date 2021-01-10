rm(list=ls())
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho_dcsbm.R")
run_suffix <- ""

paramMat <- cbind(500, 100, seq(0.1, 1, length.out = 15), 0.4, 0.1, 0.5)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)
K <- 3

trials <- 50
ncores <- 10
doMC::registerDoMC(cores = ncores)
rm(list = c("vec1", "vec2", "vec3", "eigen_mat"))

###############################

rule <- function(vec){
  n <- vec["n"]; L <- vec["L"]; rho <- vec["rho"]
  mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
  
  degree_vec <- c(seq(0.2, 0.3, length.out = 0.2*n), seq(0.7, 0.8, length.out = 0.2*n), 
                  rep(0.5, 0.1*n), 
                  seq(0.1, 0.2, length.out = 0.2*n),  seq(0.9, 1, length.out = 0.3*n))
   
  prob_mat1 <- networkSoSD:::.mult_mat_vec(networkSoSD:::.mult_vec_mat(degree_vec, prob_mat1), degree_vec)
  prob_mat2 <- networkSoSD:::.mult_mat_vec(networkSoSD:::.mult_vec_mat(degree_vec, prob_mat2), degree_vec)
  
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list, prob_list = prob_list)
}

i <- 10; y <- 1; set.seed(y); dat <- rule(paramMat[i,])
prob_agg_network <- networkSoSD::aggregate_networks(dat$prob_list, method = "ss")
svd_mat <- .svd_projection(prob_agg_network, K = 3, weighted = F)
svd_mat <- t(apply(svd_mat, 1, function(x){x/.l2norm(x)}))
# image(t(svd_mat))
tmp <- prob_agg_network; diag(tmp) <- 0; quantile(tmp)
quantile(diag(prob_agg_network)); plot(diag(prob_agg_network)); plot(colSums(prob_agg_network))

####
agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
# image(agg_network)
# plot(diag(agg_network))
tmp <- agg_network; diag(tmp) <- 0; quantile(tmp)
quantile(diag(agg_network))
svd_mat <- .svd_projection(agg_network, K = 3, weighted = F)
svd_mat <- t(apply(svd_mat, 1, function(x){x/.l2norm(x)}))
# image(t(svd_mat))
networkSoSD::spectral_clustering(agg_network, K = K, weighted = F, row_normalize = T)


agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
# image(agg_network)
# plot(diag(agg_network))
tmp <- agg_network; diag(tmp) <- 0; quantile(tmp)
quantile(diag(agg_network))
svd_mat <- .svd_projection(agg_network, K = 3, weighted = F)
svd_mat <- t(apply(svd_mat, 1, function(x){x/.l2norm(x)}))
# image(t(svd_mat))
# zz <- RSpectra::svds(agg_network, 10); plot(zz$d)

tmp_mat <- agg_network
diag(tmp_mat) <- colSums(tmp_mat)/nrow(tmp_mat)
# plot(diag(tmp_mat))
svd_mat <- .svd_projection(tmp_mat, K = 3, weighted = F)
svd_mat <- t(apply(svd_mat, 1, function(x){x/.l2norm(x)}))
image(t(svd_mat))







