rm(list=ls())
library(simulation)
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho_same.R")
run_suffix <- "_same"

.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,2,2)/3
vec2 <- c(-2,-1,2)/3
vec3 <- c(2,-2,1)/3
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1, 0.3, 0.4)) %*% t(eigen_mat)
K <- 3

zz <- solve(B1, rep(1,3)); zz <- zz/sum(zz)
paramMat <- cbind(500, 100, seq(0.005, 0.025, length.out = 11), zz[1], zz[2], zz[3])
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

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
  prob_list <- lapply(1:L, function(i){prob_mat1})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list, prob_list = prob_list)
}

criterion <- function(dat, vec, y){
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of adding
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of ss
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of flattening
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4 <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = F)
  
  # try less aggressive debiasing
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias2")
  res5 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F, row_normalize = T)
  
  ### now all the weighted versions
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of adding
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of ss
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of flattening
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4b <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = T)
  
  # try less aggressive debiasing
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias2")
  res5b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T, row_normalize = T)
  
  ### now the greedy method
  res6 <- networkSoSD::greedy_clustering(dat$adj_list, K = K)$cluster
  
  list(res_ss_debias_F = res1, res_ss_debias_T = res1b, 
       res_sum_F = res2, res_sum_T = res2b, 
       res_ss_F = res3, res_ss_T = res3b,
       res_flat_F = res4, res_flat_T = res4b,
       res_ss_debias2_F = res5, res_ss_debias2_T = res5b,
       res_greedy = res6)
}

## i <- 1; y <- 1; set.seed(y); zz <- criterion(rule(paramMat[i,]), paramMat[i,], y); zz

#########################


res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = ncores, as_list = T,
                                        filepath = "../results/simulation_rho_tmp.RData",
                                        verbose = T)

save.image(paste0("../results/simulation_rho", run_suffix, ".RData"))
