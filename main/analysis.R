rm(list=ls()); set.seed(10)

library(networkSoSD); library(org.Mmu.eg.db); library(RSpectra); library(irlba)

load("../../../data/bakken_pnas/pnas.RData")
session_info <- devtools::session_info()
date_of_run <- Sys.time()
source_code_info <- readLines("../main/analysis.R")
run_suffix <- "_revision"

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

##########################

adj_list <- lapply(dat_list, function(mat){
  adj_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  adj_mat[which(abs(mat) >= 0.15)] <- 1
  adj_mat[keep_idx, keep_idx]
})

degree_list <- lapply(adj_list, function(adj_mat){
  colSums(adj_mat)
})

degree_mat <- do.call(rbind, degree_list)
degree_vec <- colSums(degree_mat)
# zero_vec <- apply(degree_mat, 2, function(x){length(which(x == 0))})
# quantile(degree_vec, probs = seq(0,1,length.out=11))
# table(zero_vec)
# idx <- which(zero_vec >= 5)
# quantile(degree_vec[idx], probs = seq(0,1,length.out=11))

keep_idx2 <- which(degree_vec >= 90) # above commented-out lines inform us of how 90 was chosen
# length(keep_idx2)
# length(keep_idx2)/nrow(adj_list[[1]])

adj_list <- lapply(adj_list, function(adj_mat){
  adj_mat[keep_idx2, keep_idx2]
})
entrez_id2 <- entrez_id[keep_idx2]
gene_name2 <- gene_name[keep_idx2]

# check power law distribution
power_law <- sapply(adj_list, function(adj_mat){
  deg_vec <- colSums(adj_mat)
  
  x <- table(deg_vec)+1
  y <- as.numeric(names(x))+1
  x <- as.numeric(x)
  
  abs(stats::cor(log(x),log(y)))
})

# lapply(adj_list, function(x){quantile(colSums(x))})

total_network <- networkSoSD::aggregate_networks(adj_list, method = "ss_debias", verbose = T)
svd_res <- RSpectra::svds(total_network, k = 30)
# abs(diff(svd_res$d))/svd_res$d[2:length(svd_res$d)]

set.seed(10)
K <- 8
clustering_res <- networkSoSD::spectral_clustering(total_network, K = K, weighted = F, row_normalize = F)

save_var <- c("adj_list", "entrez_id", "gene_name", "entrez_id2", "gene_name2", "power_law",
              "total_network", "svd_res", "clustering_res", "K", "session_info",
              "date_of_run", "source_code_info", "run_suffix")
var_list <- ls()
var_list <- var_list[!var_list %in% save_var]
rm(list = var_list)
save.image(paste0("../results/main_analysis", run_suffix, ".RData"))
