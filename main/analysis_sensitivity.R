rm(list=ls()); set.seed(10)

library(networkSoSD); library(org.Mmu.eg.db); library(RSpectra); library(irlba)

load("../../../data/bakken_pnas/pnas.RData")
session_info <- devtools::session_info()
date_of_run <- Sys.time()
source_code_info <- readLines("../main/analysis_sensitivity.R")

# manage gene names
gene_name <- read.csv("../../../data/bakken_pnas/All_human_genes.txt", header = F)

library(org.Mmu.eg.db)
rhesus_all <- as.list(org.Mmu.eg.db::org.Mmu.egALIAS2EG)
rhesus_all <- rhesus_all[!is.na(rhesus_all)]
entrez_id <- sapply(1:nrow(gene_name), function(i){
  if(i %% floor(nrow(gene_name)/10) == 0) cat('*')
  
  id <- gene_name[i,1]
  idx <- which(names(rhesus_all) == id)
  
  if(length(idx) == 0){
    return(NA)
  } else {
    rhesus_all[[idx[1]]][1]
  }
})
keep_idx <- which(!is.na(entrez_id))
entrez_id <- entrez_id[keep_idx]
gene_name <- gene_name[keep_idx,1]

##################################

analysis_function <- function(dat_list, 
                              keep_idx = 1:ncol(dat_list[[1]]),
                              cor_threshold = 0.15, 
                              degree_threshold = 90,
                              num_clusters = 8){
  
  print("Converting")
  adj_list <- lapply(dat_list, function(mat){
    adj_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    adj_mat[which(abs(mat) >= cor_threshold)] <- 1
    adj_mat[keep_idx, keep_idx]
  })
  
  print("Filtering")
  degree_list <- lapply(adj_list, function(adj_mat){
    colSums(adj_mat)
  })
  
  degree_mat <- do.call(rbind, degree_list)
  degree_vec <- colSums(degree_mat)
  keep_idx2 <- which(degree_vec >= degree_threshold) 
  
  adj_list <- lapply(adj_list, function(adj_mat){
    Matrix::Matrix(adj_mat[keep_idx2, keep_idx2], sparse = T)
  })
  entrez_id2 <- entrez_id[keep_idx2]
  gene_name2 <- gene_name[keep_idx2]
  
  print("Aggregating networks")
  total_network <- networkSoSD::aggregate_networks(adj_list, method = "ss_debias", verbose = T)

  print("Clustering")
  set.seed(10)
  clustering_res <- networkSoSD::spectral_clustering(total_network, K = num_clusters, weighted = F, row_normalize = F)
  
  list(clustering_res = clustering_res,
       keep_idx2 = keep_idx2,
       entrez_id2 = entrez_id2, 
       gene_name2 = gene_name2,
       cor_threshold = cor_threshold, 
       degree_threshold = degree_threshold,
       num_clusters = num_clusters)
}

##############################

# we choose these values since 0.15^(1/6) = 0.72
# so we try 0.65 and 0.8, which correspond to 0.65^6 = 0.08
# and 0.8^6 = 0.26
df_param <- expand.grid(c(0.08, 0.15, 0.26), c(7:9))
colnames(df_param) <- c("cor_threshold", "num_clusters")

res_list <- lapply(1:nrow(df_param), function(x){
  print(paste0("On row ", x))
  set.seed(10)
  analysis_function(dat_list, 
                    keep_idx = keep_idx, 
                    cor_threshold = as.numeric(df_param$cor_threshold[x]), 
                    degree_threshold = 90,
                    num_clusters = as.numeric(df_param$num_clusters[x]))
})

save(dat_list, session_info, date_of_run, source_code_info,
     keep_idx, entrez_id, gene_name,
     df_param, res_list,
     file = "../results/analysis_sensitivity.RData")
