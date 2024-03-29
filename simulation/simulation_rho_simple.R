rm(list=ls())
library(customSimulator)
library(networkSoSD)

session_info <- devtools::session_info()
date_of_run <- Sys.time()
source("../simulation/Codes_Spectral_Matrix_Paul_Chen_AOS_2020.r")
source_code_info <- readLines("../simulation/simulation_rho_simple.R")
run_suffix <- ""

df_param <- cbind(200, 30, seq(0.02, 0.06, length.out = 9), 0.5, 0.5)
colnames(df_param) <- c("n", "L", "rho", "mem_prop1", "mem_prop2")
df_param <- as.data.frame(df_param)

.l2norm <- function(x){sqrt(sum(x^2))}
B1 <- matrix(c(3/4, sqrt(3)/8, sqrt(3)/8, 1/2), 2, 2)
B2 <- matrix(c(7/8, 3*sqrt(3)/8, 3*sqrt(3)/8, 1/8), 2, 2)
K <- ncol(B1)

ntrials <- 100
ncores <- 1

###############################

rule <- function(vec){
  n <- as.numeric(vec["n"]); L <- as.numeric(vec["L"]); rho <- as.numeric(vec["rho"])
  mem_prop1 <- as.numeric(vec["mem_prop1"]); mem_prop2 <- as.numeric(vec["mem_prop2"])
  
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
  set.seed(10)
  res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of adding
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  set.seed(10)
  res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of ss
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  set.seed(10)
  res3 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of flattening
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  set.seed(10)
  res4 <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = F)
  
  # yuguo's methods
  set.seed(10)
  res5 <- networkSoSD::lmfo(dat$adj_list, k = K)
  
  set.seed(10)
  reg_val <- max(sapply(dat$adj_list, function(x){4*max(abs(RSpectra::eigs_sym(x, k = min(K,5))$values))}))
  tmp <- networkSoSD::coreg(dat$adj_list, k = K, beta = reg_val,
                            verbose = F, max_iter = 50)
  res6 <- tmp[[length(dat$adj_list)+1]]
  
  list(res_ss_debias_F = res1, res_sum_F = res2, 
       res_ss_F = res3, res_flat_F = res4,
       chen_linked = res5, chen_coreg = res6)
}

## i <- 1; y <- 1; set.seed(y); zz <- criterion(rule(df_param[i,]), df_param[i,], y); zz

#########################


res <- customSimulator::simulator(rule = rule, criterion = criterion,
                                  df_param = df_param, ntrials = ntrials,
                                  cores = ncores,
                                  filepath = "../results/simulation_rho_simple_tmp.RData",
                                  required_packages = c("networkSoSD", "irlba", "clue", "stats", "Matrix"),
                                  verbose = T)

save.image(paste0("../results/simulation_rho_simple", run_suffix, ".RData"))
