rm(list=ls()); set.seed(10)
load("../results/main_analysis.RData")
library(clusterProfiler); library(igraph); library(org.Mmu.eg.db)

color_func <- function(alpha = 0.2){
  c(grDevices::rgb(240/255, 228/255, 66/255, alpha), #yellow
    grDevices::rgb(86/255, 180/255, 233/255, alpha), #skyblue
    grDevices::rgb(0/255, 158/255, 115/255, alpha), #bluish green
    grDevices::rgb(0/255, 114/255, 178/255,alpha), #blue
    grDevices::rgb(230/255, 159/255, 0/255,alpha), #orange
    grDevices::rgb(150/255, 150/255, 150/255, alpha), #gray
    grDevices::rgb(189/255, 57/255, 60/255, alpha), # red
    grDevices::rgb(245/255, 234/255, 204/255, alpha), #pale
    grDevices::rgb(204/255, 169/255, 221/255, alpha)) #lightpurple
}
color_name_vec <- c("yellow", "skyblue", "green", "blue", "orange", "gray", "red", "pale", "lightpurple")
color_vec <- color_func(1)

## first do some housekeeping
# manually change the labeling of the clusters
desired_order <- c(5,8,1,7,4,3,6,2)
clustering_res2 <- rep(NA, length(clustering_res))
for(i in 1:length(desired_order)){
  clustering_res2[which(clustering_res == desired_order[i])] <- i
}
table(clustering_res2)
table(clustering_res2)/nrow(adj_list[[1]])

## print out some stats
quantile(power_law)

## determine how many unique clusters there
component_num <- sapply(adj_list, function(adj_mat){
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  tmp <- igraph::components(g)
  idx <- which(tmp$csize <= 5)
  c(tmp$no, length(which(tmp$membership %in% idx))/nrow(adj_mat))
})

tmp <- networkSoSD::aggregate_networks(adj_list, method = "sum", verbose = T)
tmp[tmp != 0] <- 1
g <- igraph::graph_from_adjacency_matrix(tmp, mode = "undirected")
tmp <- igraph::components(g)
tmp$no
idx <- which(tmp$csize <= 5)
length(which(tmp$membership %in% idx))/nrow(adj_list[[1]])

##########

png(paste0("../figures/Writeup4_pnas_svd.png"), height = 1500, width = 2500, units = "px", res = 300)
par(mfrow = c(1,2))
plot(svd_res$d, pch = 16, xlab = "Index", ylab = "Singular value magnitude",
     main = "Singular values of aggregated matrix")
lines(rep(8.5,2), c(-1e2,1e9), col = color_vec[which(color_name_vec == "red")], lwd = 2, lty = 2)

vec <- -diff(svd_res$d)/(svd_res$d[-1])
plot(vec, xlim = c(1,length(vec)), ylim = range(vec), pch = 16, xlab = "Index", ylab = "Relative magnitude of singluar value diff.",
     main = "Singular value gap of\naggregated matrix")
lines(rep(8.5,2), c(-1e2,1e9), col = color_vec[which(color_name_vec == "red")], lwd = 2, lty = 2)
graphics.off()

#########

set.seed(10)
umap_embedding <- Seurat::RunUMAP(networkSoSD:::.mult_mat_vec(svd_res$u, svd_res$d), verbose = F)@cell.embeddings
col_umap <- color_vec[c(9,2:8)][clustering_res2]

png(paste0("../figures/Writeup4_umap.png"), height = 1500, width = 1500, units = "px", res = 300)
plot(NA, xlim = range(umap_embedding[,1]), ylim = range(umap_embedding[,2]),
     xlab = "UMAP dimension 1", ylab = "UMAP dimension 2",
     main = "UMAP of spectral embedding", asp = T)
idx <- which(clustering_res2 == 8)
points(umap_embedding[idx,1], umap_embedding[idx,2], pch = 16, col = col_umap[idx])
points(umap_embedding[-idx,1], umap_embedding[-idx,2], pch = 16, col = col_umap[-idx])
graphics.off()

###################

low_dim_mat <- do.call(cbind, lapply(1:length(adj_list), function(i){
  print(i)
  tmp <- networkSoSD:::.svd_truncated(adj_list[[i]], K = K, symmetric = T)
  networkSoSD:::.mult_mat_vec(tmp$u, tmp$d)
}))

set.seed(10)
umap_embedding2 <- Seurat::RunUMAP(low_dim_mat, verbose = F)@cell.embeddings
col_umap <- color_vec[c(9,2:8)][clustering_res2]

png(paste0("../figures/Writeup4_umap2.png"), height = 1500, width = 1500, units = "px", res = 300)
plot(NA, xlim = range(umap_embedding2[,1]), ylim = range(umap_embedding2[,2]),
     xlab = "UMAP dimension 1", ylab = "UMAP dimension 2",
     main = "UMAP of spectral embedding", asp = T)
idx <- which(clustering_res2 == 8)
points(umap_embedding2[idx,1], umap_embedding2[idx,2], pch = 16, col = col_umap[idx])
points(umap_embedding2[-idx,1], umap_embedding2[-idx,2], pch = 16, col = col_umap[-idx])
graphics.off()

####################

time_stamp <- c("0M", "12M", "3M", "48M", "E120", "E40", "E50", "E70", "E80", "E90")
time_order <- c(7,9,8,10,6,1,2,3,4,5)
# alternatively, run the following lines:
# file_vec <- list.files("/raid6/Fuchen/monkey_data/PretimeA/")
# file_vec <- file_vec[grep("All.*new.*", file_vec)]
# time_stamp <- as.vector(sapply(file_vec, function(x){
#   tmp <- rev(strsplit(x, split = "_")[[1]])[1]
#   strsplit(tmp, split = "\\.")[[1]][1]
# }))

idx_odd <- seq(1, length(clustering_res2), by = 10)
clustering_res2b <- clustering_res2[idx_odd]
gene_idx <- unlist(lapply(1:K, function(x){idx_odd[which(clustering_res2b == x)]}))
clockwise90 <- function(a) { t(a[nrow(a):1,]) } 
for(i in 1:length(adj_list)){
  print(i)
  tmp <- adj_list[[i]][gene_idx, gene_idx]
  diag(tmp) <- -1
  
  png(paste0("../figures/pnas_adj", time_order[i], ".png"), height = 1500, width = 1300, units = "px", res = 300)
  par(mar = c(2, 2, 4, 0.6))
  image(clockwise90(tmp), main = paste0("Network for time ", time_stamp[i]), 
        col = c("white", hcl.colors(2, "Cividis")),
        breaks = c(-1.5, -.5, .5, 1.5))
  
  for(j in 1:(K-1)){
    len <- length(which(clustering_res2b <= j))/length(clustering_res2b)
    lines(c(-1e4,1e4), rep(1-len,2), lwd = 2, lty = 2, col = 'white')
    lines(rep(len,2), c(-1e4,1e4), lwd = 2, lty = 2, col = 'white')
  }
  
  graphics.off()
}

# gene_idx <- unlist(lapply(1:K, function(x){which(clustering_res2 == x)}))
# clockwise90 <- function(a) { t(a[nrow(a):1,]) } 
# png(paste0("../figures/Writeup4_pnas_adj_all.png"), height = 3000, width = 3000, units = "px", res = 300)
# image(clockwise90(total_network[gene_idx, gene_idx]), main = "All", 
#       col = hcl.colors(12, "Cividis"))
# for(j in 1:(K-1)){
#   len <- length(which(clustering_res2 <= j))/length(clustering_res2)
#   lines(c(-1e4,1e4), rep(1-len,2), lwd = 2, lty = 2)
#   lines(rep(len,2), c(-1e4,1e4), lwd = 2, lty = 2)
# }
# graphics.off()

# tmp <- networkSoSD::aggregate_networks(adj_list, method = 'sum')
# gene_idx <- unlist(lapply(1:K, function(x){which(clustering_res2 == x)}))
# tmp <- tmp[gene_idx, gene_idx]
# diag(tmp) <- -1
# clockwise90 <- function(a) { t(a[nrow(a):1,]) } 
# png(paste0("../figures/pnas_adj_all.png"), height = 3000, width = 3000, units = "px", res = 300)
# image(clockwise90(tmp), main = "All", 
#       col = c("white", hcl.colors(2, "Cividis")),
#       breaks = c(-1.5, -.5, .5, 1.5))
# for(j in 1:(K-1)){
#   len <- length(which(clustering_res2 <= j))/length(clustering_res2)
#   lines(c(-1e4,1e4), rep(1-len,2), lwd = 2, lty = 2, col = 'white')
#   lines(rep(len,2), c(-1e4,1e4), lwd = 2, lty = 2, col = 'white')
# }
# graphics.off()
# 
# tmp_svd <- networkSoSD:::.svd_truncated(tmp, K = K)
# set.seed(10)
# umap_embedding3 <- Seurat::RunUMAP(networkSoSD:::.mult_mat_vec(tmp_svd$u, tmp_svd$d), verbose = F)@cell.embeddings
# col_umap <- color_vec[c(9,2:8)][clustering_res2]
# 
# png(paste0("../figures/Writeup4_umap3.png"), height = 1500, width = 1500, units = "px", res = 300)
# plot(NA, xlim = range(umap_embedding3[,1]), ylim = range(umap_embedding3[,2]),
#      xlab = "UMAP dimension 1", ylab = "UMAP dimension 2",
#      main = "UMAP of spectral embedding", asp = T)
# idx <- which(clustering_res2 == 8)
# points(umap_embedding3[idx,1], umap_embedding3[idx,2], pch = 16, col = col_umap[idx])
# points(umap_embedding3[-idx,1], umap_embedding3[-idx,2], pch = 16, col = col_umap[-idx])
# graphics.off()

###################

ego_list <- vector("list", length = K)
for(i in 1:K){
  set.seed(10)
  ego_list[[i]] <- clusterProfiler::enrichGO(gene          = entrez_id2[which(clustering_res2 == i)],
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

# see:
# # http://supfam.org/SUPERFAMILY/cgi-bin/go.cgi
# # https://www.biostars.org/p/237816/


###########################################

# # analyze the connectivity
# connectivity_list <- lapply(1:length(adj_list), function(x){
#   mat <- matrix(0, K, K)
#   for(i in 1:K){
#     for(j in 1:i){
#       idx1 <- which(clustering_res == i)
#       idx2 <- which(clustering_res == j)
#       mat[i,j] <- mean(adj_list[[x]][idx1,idx2])
#       mat[j,i] <- mat[i,j]
#     }
#   }
#   
#   mat
# })
# 
# lapply(connectivity_list, function(x){round(x, 2)})
# 
# lapply(connectivity_list, function(x){
#   eigen(x)$values
# })

