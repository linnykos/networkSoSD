rm(list=ls())
library(simulation)
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho.R")

paramMat <- cbind(500, 100, seq(0.1, 0.15, length.out = 11), 0.45, 0.1, 0.45)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)
K <- 3

trials <- 100
ncores <- 15
doMC::registerDoMC(cores = ncores)
rm(list = c("vec1", "vec2", "vec3", "eigen_mat"))

###############################

rule <- function(vec){
  n <- vec["n"]; L <- vec["L"]; rho <- vec["rho"]
  mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list)
}

criterion <- function(dat, vec, y){
  agg_network <- aggregate_networks(dat$adj_list, method = "ss_debias")
  res1 <- spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of adding
  agg_network <- aggregate_networks(dat$adj_list, method = "sum")
  res2 <- spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of ss
  agg_network <- aggregate_networks(dat$adj_list, method = "ss")
  res3 <- spectral_clustering(agg_network, K = K, weighted = F)
  
  ### now all the weighted versions
  # try naive method of adding
  agg_network <- aggregate_networks(dat$adj_list, method = "sum")
  res4 <- spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of ss
  agg_network <- aggregate_networks(dat$adj_list, method = "ss")
  res5 <- spectral_clustering(agg_network, K = K, weighted = T)
  
  list(res_ss_debias_F = res1, 
       res_sum_F = res2, res_ss_F = res3, 
       res_sum_T = res4, res_ss_T = res5)
}

## i <- 6; y <- 1; set.seed(y); zz1 <- criterion(rule(paramMat[i,]), paramMat[i,], y); zz1

#########################


res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = ncores, as_list = T,
                                        filepath = "../results/simulation_rho_tmp.RData",
                                        verbose = T)

save.image("../results/simulation_rho.RData")
