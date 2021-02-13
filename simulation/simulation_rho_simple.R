rm(list=ls())
library(simulation)
library(networkSoSD)

session_info <- devtools::session_info()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho_simple.R")
run_suffix <- ""

paramMat <- cbind(200, 30, seq(0.02, 0.06, length.out = 9), 0.5, 0.5)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2")

.l2norm <- function(x){sqrt(sum(x^2))}
B1 <- matrix(c(3/4, sqrt(3)/8, sqrt(3)/8, 1/2), 2, 2)
B2 <- matrix(c(7/8, 3*sqrt(3)/8, 3*sqrt(3)/8, 1/8), 2, 2)
K <- ncol(B1)

trials <- 100
ncores <- 10
doMC::registerDoMC(cores = ncores)

###############################

rule <- function(vec){
  n <- vec["n"]; L <- vec["L"]; rho <- vec["rho"]
  mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(2, n-length(membership_vec)))
  
  prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list)
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
  
  list(res_ss_debias_F = res1, res_sum_F = res2, 
       res_ss_F = res3, res_flat_F = res4)
}

## i <- 1; y <- 1; set.seed(y); zz <- criterion(rule(paramMat[i,]), paramMat[i,], y); zz

#########################


res <- simulation::simulation_generator(rule = rule, criterion = criterion,
                                        paramMat = paramMat, trials = trials,
                                        cores = ncores, as_list = T,
                                        filepath = "../results/simulation_rho_simple_tmp.RData",
                                        verbose = T)

save.image(paste0("../results/simulation_rho_simple", run_suffix, ".RData"))
