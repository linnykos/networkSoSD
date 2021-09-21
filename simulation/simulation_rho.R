rm(list=ls())
library(customSimulator)
library(networkSoSD)

session_info <- devtools::session_info()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho.R")
run_suffix <- ""

df_param <- cbind(3, 500, 100, seq(0.025, 0.2, length.out = 15), 0.4, 0.1, 0.5)
colnames(df_param) <- c("K", "n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")
df_param <- as.data.frame(df_param)

ntrials <- 100
ncores <- 4
worker_variables <- NA

###############################

generator <- function(vec, worker_variables){
  .l2norm <- function(x){sqrt(sum(x^2))}
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
  B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)

  n <- as.numeric(vec["n"]); L <- as.numeric(vec["L"]); rho <- as.numeric(vec["rho"])
  K <- as.numeric(vec["K"])
  mem_prop1 <- as.numeric(vec["mem_prop1"]); mem_prop2 <- as.numeric(vec["mem_prop2"])
  mem_prop3 <- as.numeric(vec["mem_prop3"])
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list)
}

executor <- function(dat, vec, y, worker_variables){
  K <- as.numeric(vec["K"])

  time_start1 <- Sys.time()
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  time_end1 <- Sys.time()
  
  # try naive method of adding
  time_start2 <- Sys.time()
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  time_end2 <- Sys.time()
  
  # try naive method of ss
  time_start3 <- Sys.time()
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  time_end3 <- Sys.time()
  
  # try naive method of flattening
  time_start4 <- Sys.time()
  set.seed(10)
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4 <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = F)
  time_end4 <- Sys.time()
  
  ### now all the weighted versions
  time_start1b <- Sys.time()
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  time_end1b <- Sys.time()
  
  # try naive method of adding
  time_start2b <- Sys.time()
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  time_end2b <- Sys.time()
  
  # try naive method of ss
  time_start3b <- Sys.time()
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  time_end3b <- Sys.time()
  
  # try naive method of flattening
  time_start4b <- Sys.time()
  set.seed(10)
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4b <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = T)
  time_end4b <- Sys.time()
  
  # yuguo's methods
  time_start5 <- Sys.time()
  set.seed(10)
  res5 <- networkSoSD::lmfo(dat$adj_list, k = K)
  time_end5 <- Sys.time()
  
  # set.seed(10)
  time_start6 <- Sys.time()
  reg_val <- max(sapply(dat$adj_list, function(x){4*max(abs(RSpectra::eigs_sym(x, k = min(K,5))$values))}))
  tmp <- networkSoSD::coreg(dat$adj_list, k = K, beta = reg_val,
               verbose = F, max_iter = 50)
  res6 <- tmp[[length(dat$adj_list)+1]]
  time_end6 <- Sys.time()
  
  list(res_ss_debias_F = res1, time_ss_debias_F = time_end1 - time_start1,
       res_ss_debias_T = res1b, time_ss_debias_T = time_end1b - time_start1b,
       res_sum_F = res2, time_sum_F = time_end2 - time_start2,
       res_sum_T = res2b, time_sum_T = time_end2b - time_start2b,
       res_ss_F = res3, time_ss_F = time_end3 - time_start3,
       res_ss_T = res3b, time_ss_T = time_end3b - time_start3b,
       res_flat_F = res4, time_flat_F = time_end4 - time_start4,
       res_flat_T = res4b, time_flat_T = time_end4b - time_start4b,
       chen_linked = res5, time_chen_linked = time_end5 - time_start5,
       chen_coreg = res6, time_chen_coreg = time_end6 - time_start6)
}

## i <- 1; y <- 1; set.seed(y); zz <- executor(generator(df_param[i,], worker_variables), df_param[i,], y, worker_variables); zz

#########################


res <- customSimulator::simulator(generator = generator, executor = executor,
                                  df_param = df_param, ntrials = ntrials,
                                  cores = ncores, 
                                  filepath = "../results/simulation_rho_tmp.RData",
                                  required_packages = "networkSoSD",
                                  verbose = T)

save.image(paste0("../results/simulation_rho", run_suffix, ".RData"))
