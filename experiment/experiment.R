rm(list=ls())
library(customSimulator)
library(networkSoSD)

session_info <- devtools::session_info()
date_of_run <- Sys.time()
source("../simulation/Codes_Spectral_Matrix_Paul_Chen_AOS_2020.r")
source_code_info <- readLines("../simulation/simulation_rho.R")
run_suffix <- ""

# df_param <- cbind(3, 500, 100, seq(0.025, 0.2, length.out = 15), 0.4, 0.1, 0.5)
df_param <- cbind(3, 30, 10, seq(0.025, 0.2, length.out = 2), 0.4, 0.1, 0.5)
colnames(df_param) <- c("K", "n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")
df_param <- as.data.frame(df_param)

# ntrials <- 50
ntrials <- 5
ncores <- 2
worker_variables <- list(lmfo = lmfo, coreg = coreg, vec = 1:10)

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
  
  ### now the greedy method
  set.seed(10)
  res5 <- networkSoSD::greedy_refinement(dat$adj_list, K = K)$cluster
  
  # yuguo's methods
  set.seed(10)
  res6 <- worker_variables$lmfo(dat$adj_list, n = nrow(dat$adj_list[[1]]), k = K)
  # res5 <- mean(worker_variables$vec)
  
  # set.seed(10)
  # reg_val <- max(sapply(dat$adj_list, function(x){4*max(abs(RSpectra::eigs_sym(x, k = min(K,5))$values))}))
  # tmp <- worker_variables$coreg(dat$adj_list, n = nrow(dat$adj_list[[1]]), k = K, beta = reg_val,
  #              verbose = F, max_iter = 50)
  # res7 <- tmp[[length(dat$adj_list)+1]]
  
  list(res_ss_debias_F = res1, res_ss_debias_T = res1b, 
       res_sum_F = res2, res_sum_T = res2b, 
       res_ss_F = res3, res_ss_T = res3b,
       res_flat_F = res4, res_flat_T = res4b,
       res_greedy = res5, chen_linked = res6) #, chen_coreg = res7)
}

############

ntrials = 5
specific_trials = NA
cores = 2
shuffle_group = NA
chunking_num = nrow(df_param)
required_packages = c("networkSoSD", "irlba", "clue", "stats", "Matrix", "RSpectra", "psych")
filepath = NA
verbose = F

# construct the scheduler
df_schedule <- customSimulator:::.construct_scheduler(nrow(df_param), ntrials, specific_trials)
pb <- function(){invisible()}

# if shuffling is used, shuffle the scheduler
if(!all(is.na(shuffle_group))){
  df_schedule <- customSimulator:::.shuffle(df_schedule, shuffle_group)
}

# create the empty list that we will populate with the results
res_all <- lapply(1:nrow(df_param), function(i){
  if(!is.na(ntrials)){
    tmp <- vector("list", length = ntrials)
    names(tmp) <- paste0("trial_", 1:ntrials)
    tmp
  } else {
    tmp <- vector("list", length = length(specific_trials))
    names(tmp) <- paste0("trial_", specific_trials)
    tmp
  }
})

# create the chunking
chunking_list <- customSimulator:::.split_rows(nrow(df_schedule), chunking_num)

fun <- function(i, df_schedule, generator, executor, worker_variables, pb){
  x <- df_schedule$row[i]
  y <- df_schedule$trial[i]
  set.seed(y)
  
  dat <- generator(df_param[x,], worker_variables)
  start_time <- proc.time()
  res <- executor(dat, df_param[x,], y, worker_variables)
  end_time <- proc.time()
  
  res <- list(result = res)
  res$elapsed_time <- end_time-start_time
  res
}

if(cores > 1) future::plan(future::multisession, workers = cores)

for(k in 1:length(chunking_list)){
  if(verbose) cat(paste0("\n", Sys.time(), ": Chunk ", k, " of ", length(chunking_list), " started!\n"))
  
  if(cores > 1){
    # set up progress bar
    if(verbose & cores > 1){
      progressr::handlers(global = T)
      pb <- progressr::progressor(along = chunking_list[[k]])
    }
    
    # parallel version
    res_tmp <- future.apply::future_lapply(chunking_list[[k]], function(i){
      fun(i, df_schedule, generator, executor, worker_variables, pb)
    }, future.globals = list(generator = generator, executor = executor,
                             df_param = df_param, df_schedule = df_schedule,
                             worker_variables = worker_variables,
                             pb = pb, verbose = verbose, fun = fun),
    future.packages = required_packages, future.seed = TRUE)
    
  } else {
    # sequential version
    if(verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
    res_tmp <- pbapply::pblapply(chunking_list[[k]], function(i){
      fun(i, df_schedule, generator, executor, worker_variables, pb)
    })
  }
  
  # copy the results in
  for(i in 1:length(chunking_list[[k]])){
    x <- df_schedule$row[chunking_list[[k]][i]]
    y <- df_schedule$trial[chunking_list[[k]][i]]
    
    if(!is.na(ntrials)){
      res_all[[x]][[y]] <- res_tmp[[i]]
    } else {
      res_all[[x]][[which(specific_trials == y)]] <- res_tmp[[i]]
    }
  }
  
  if(!is.na(filepath)) save(res_all, file = filepath)
}

if(cores > 1) future::plan(future::sequential)
if(verbose){
  progressr::handlers(global = F)
}

names(res_all) <- paste0("row_", 1:length(res_all))
res_all

