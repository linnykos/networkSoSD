rm(list=ls())
library(customSimulator)
library(networkSoSD)

session_info <- devtools::session_info()
date_of_run <- Sys.time()

df_param <- cbind(3, 500, 100, 0.2, seq(0, 0.5, length.out = 11), 
                  0.4, 0.1, 0.5)
colnames(df_param) <- c("K", "n", "L", "rho", "switch_prob",
                        "mem_prop1", "mem_prop2", "mem_prop3")
df_param <- as.data.frame(df_param)

ntrials <- 50
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
  K <- as.numeric(vec["K"]); switch_prob <- as.numeric(vec["switch_prob"])
  stationary_vec <- as.numeric(c(vec["mem_prop1"], vec["mem_prop2"], vec["mem_prop3"]))
  mem_prop1 <- stationary_vec[1]; mem_prop2 <- stationary_vec[2]
  mem_prop3 <- stationary_vec[3]
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  transition <- networkSoSD::compute_markov_transition(stationary_vec, switch_prob)
  
  prob_list <- lapply(1:L, function(i){
    membership_vec2 <- membership_vec
    for(k in 1:3){
      idx <- which(membership_vec == k)
      membership_vec2[idx] <- sample(c(1:3), size = length(idx), 
                                     replace = T, prob = transition[k,])
    }
    
    if(i <= L/2){
      networkSoSD::compute_prob_mat(rho*B1, membership_vec2)
    } else {
      networkSoSD::compute_prob_mat(rho*B2, membership_vec2)
    }
  })
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list)
}

executor <- function(dat, vec, y, worker_variables){
  K <- as.numeric(vec["K"])
  
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of adding
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of ss
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3 <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = F)
  
  # try naive method of flattening
  set.seed(10)
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4 <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = F)
  
  ### now all the weighted versions
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
  res1b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of adding
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "sum")
  res2b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of ss
  set.seed(10)
  agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
  res3b <- networkSoSD::spectral_clustering(agg_network, K = K, weighted = T)
  
  # try naive method of flattening
  set.seed(10)
  flat_mat <- networkSoSD::flatten(dat$adj_list)
  res4b <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = T)
  
  # yuguo's methods
  set.seed(10)
  res5 <- networkSoSD::lmfo(dat$adj_list, k = K)
  
  # set.seed(10)
  reg_val <- max(sapply(dat$adj_list, function(x){4*max(abs(RSpectra::eigs_sym(x, k = min(K,5))$values))}))
  tmp <- networkSoSD::coreg(dat$adj_list, k = K, beta = reg_val,
                            verbose = F, max_iter = 50)
  res6 <- tmp[[length(dat$adj_list)+1]]
  
  list(res_ss_debias_F = res1, res_ss_debias_T = res1b,
       res_sum_F = res2, res_sum_T = res2b,
       res_ss_F = res3, res_ss_T = res3b,
       res_flat_F = res4, res_flat_T = res4b,
       chen_linked = res5, 
       chen_coreg = res6)
}

## i <- 1; y <- 1; set.seed(y); zz <- executor(generator(df_param[i,], worker_variables), df_param[i,], y, worker_variables); zz

#########################


res <- customSimulator::simulator(generator = generator, executor = executor,
                                  df_param = df_param, ntrials = ntrials,
                                  cores = ncores, 
                                  filepath = "../results/simulation_rho_switching_tmp2.RData",
                                  required_packages = "networkSoSD",
                                  verbose = T)

save.image(paste0("../results/simulation_rho_switching2", run_suffix, ".RData"))
