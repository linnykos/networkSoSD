rm(list=ls())
set.seed(10)
load("../../pnas/pnas.RData")

# manage gene names
gene_name <- read.csv("../../pnas/All_human_genes.txt", header = F)

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

adj_list <- lapply(dat_list, function(dat_list){
  adj_mat <- matrix(0, nrow = nrow(dat_list), ncol = ncol(dat_list))
  adj_mat[which(dat_list >= 0.15)] <- 1
  adj_mat[keep_idx, keep_idx]
})

degree_list <- lapply(adj_list, function(adj_mat){
  colSums(adj_mat)
})

# lapply(degree_list, table)
# idx_vec <- sort(unique(unlist(lapply(degree_list, function(vec){
#   which(vec <= 1)
# }))))
# length(idx_vec)

degree_mat <- do.call(rbind, degree_list)
degree_vec <- colSums(degree_mat)
zero_vec <- apply(degree_mat, 2, function(x){length(which(x == 0))})
quantile(degree_vec, probs = seq(0,1,length.out=11))
table(zero_vec)
# quantile(zero_vec[which(degree_vec <= 50)])
idx <- which(zero_vec >= 5)
quantile(degree_vec[idx], probs = seq(0,1,length.out=11))

keep_idx2 <- which(degree_vec >= 90)
length(keep_idx2)
length(keep_idx2)/nrow(adj_list[[1]])

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
  
  stats::cor(log(x),log(y))^2
})
quantile(power_law)

lapply(adj_list, function(x){quantile(colSums(x))})

total_network <- networkSoSD::aggregate_networks(adj_list, method = "ss_debias", verbose = T)
svd_res <- RSpectra::svds(total_network, k = 30)
abs(diff(svd_res$d))/svd_res$d[2:length(svd_res$d)]

set.seed(10)
K <- 8
clustering_res <- networkSoSD::spectral_clustering(total_network, K = K, weighted = F)

########################

# manually change the labeling of the clusters
desired_order <- c(6,1,4,3,8,7,5,2)
clustering_res2 <- rep(NA, length(clustering_res))
for(i in 1:length(desired_order)){
  clustering_res2[which(clustering_res == desired_order[i])] <- i
}
table(clustering_res2)
table(clustering_res2)/nrow(adj_list[[1]])

########################

ego_list <- vector("list", length = K)
for(i in 1:K){
  set.seed(10)
  ego_list[[i]] <- clusterProfiler::enrichGO(gene         = entrez_id2[which(clustering_res2 == i)],
                                             universe      = entrez_id,
                                             OrgDb         = org.Mmu.eg.db,
                                             ont           = "ALL",
                                             pAdjustMethod = "BH",
                                             pvalueCutoff  = 0.05,
                                             qvalueCutoff  = 0.05,
                                             readable      = TRUE)
}

ego_summary <- vector("list", length = K)
for(i in 1:K){
  tmp <- ego_list[[i]]@result
  if(length(tmp) > 0){
    tmp <- tmp[,c("Count", "Description", "pvalue")]
    ego_summary[[i]] <- tmp
  } else {
    ego_summary[[i]] <- NA
  }
}
ego_summary

# ggo_list <- vector("list", length = K)
# for(i in 1:K){
#   set.seed(10)
#   tmp <- clusterProfiler::groupGO(gene     = entrez_id2[which(clustering_res == i)],
#                                             OrgDb    = org.Mmu.eg.db,
#                                             ont      = "MF",
#                                             level    = 6,
#                                             readable = TRUE)
#   ggo_list[[i]] <- head(tmp[order(tmp$Count, decreasing = T),-5])
# }
# ggo_list
# # http://supfam.org/SUPERFAMILY/cgi-bin/go.cgi
# # https://www.biostars.org/p/237816/

##########################################

png(paste0("../figures/Writeup3_pnas_svd.png"), height = 1500, width = 2500, units = "px", res = 300)
par(mfrow = c(1,2))
plot(svd_res$d, pch = 16, xlab = "Index", ylab = "Eigenvalue magnitude")
lines(rep(8.5,2), c(-1e2,1e9), col = "red", lwd = 2, lty = 2)
plot(log(svd_res$d), ylim = c(min(log(svd_res$d)), 13.1), pch = 16, xlab = "Index", ylab = "Log of eigenvalue magnitude")
lines(rep(8.5,2), c(-1e2,1e9), col = "red", lwd = 2, lty = 2)
graphics.off()

file_vec <- list.files("/raid6/Fuchen/monkey_data/PretimeA/")
file_vec <- file_vec[grep("All.*new.*", file_vec)]
time_stamp <- as.vector(sapply(file_vec, function(x){
  tmp <- rev(strsplit(x, split = "_")[[1]])[1]
  strsplit(tmp, split = "\\.")[[1]][1]
}))

gene_idx <- unlist(lapply(1:K, function(x){which(clustering_res2 == x)}))
clockwise90 <- function(a) { t(a[nrow(a):1,]) } 
for(i in 1:length(dat_list)){
  print(i)
  png(paste0("../figures/Writeup3_pnas_adj", i, ".png"), height = 3000, width = 3000, units = "px", res = 300)
  image(clockwise90(adj_list[[i]][gene_idx, gene_idx]), main = time_stamp[i])
  
  for(j in 1:(K-1)){
    len <- length(which(clustering_res2 <= j))/length(clustering_res2)
    lines(c(-1e4,1e4), rep(1-len,2), lwd = 2, lty = 2)
    lines(rep(len,2), c(-1e4,1e4), lwd = 2, lty = 2)
  }
  
  graphics.off()
}

###########################################

# analyze the connectivity
connectivity_list <- lapply(1:length(adj_list), function(x){
  mat <- matrix(0, K, K)
  for(i in 1:K){
    for(j in 1:i){
      idx1 <- which(clustering_res == i)
      idx2 <- which(clustering_res == j)
      mat[i,j] <- mean(adj_list[[x]][idx1,idx2])
      mat[j,i] <- mat[i,j]
    }
  }
  
  mat
})

lapply(connectivity_list, function(x){round(x, 2)})

lapply(connectivity_list, function(x){
  eigen(x)$values
})

