rm(list=ls())
load("../results/main_analysis_revision.RData")
names(adj_list) <- c("0M", "12M", "3M", "48M", "E120", "E40", "E50", "E70", "E80", "E90")

# first compute all the average matrices according to the clustering_res,
# so we know which genes fall into which clustering in the correct ordering
cor_list <- lapply(adj_list, function(adj_mat){
  K <- max(clustering_res)
  cor_mat <- matrix(0, K, K)
  for(i in 1:K){
    idx1 <- which(clustering_res == i)
    for(j in i:K){
      idx2 <- which(clustering_res == j)
      tmp <- adj_mat[idx1,idx2]
      cor_mat[i,j] <- mean(adj_mat[idx1,idx2])
      cor_mat[j,i] <- cor_mat[i,j]
    }
  }
  cor_mat
})
names(cor_list) <- names(adj_list)

for(i in 1:length(cor_list)){
  print(names(cor_list)[i])
  print(round(cor_list[[i]],2))
  print("===")
}

new_clustering_res <- clustering_res
order_vec <- c(3, 7, 1, 2, 6, 8, 4, 5)
for(i in 1:length(order_vec)){
  new_clustering_res[which(clustering_res == i)] <- order_vec[i]
}

gene_df <- data.frame(name = gene_name2, cluster = new_clustering_res)
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("gene_df")]
rm(list = ls_vec)
########################

load("../results/analysis_sensitivity.RData")

for(i in c(1:4,6:9)){
  name_intersect <- intersect(res_list[[i]]$gene_name2, res_list[[5]]$gene_name2)
  clust1 <- res_list[[i]]$clustering_res[which(res_list[[i]]$gene_name2 %in% name_intersect)]
  clust2 <- res_list[[5]]$clustering_res[which(res_list[[5]]$gene_name2 %in% name_intersect)]
  
  png(paste0("../figures/pnas_sensitivity", i, ".png"), 
      height = 1500, 
      width = 1500, 
      units = "px", 
      res = 300)
  par(mar = c(4,4,5,0.5))
  plot_table(clust1, clust2, 
             main = paste0("CorThres=", df_param[i,"cor_threshold"], " and K=", df_param[i,"num_clusters"], "\nvs. Baseline"),
             normalize_row = F,
             col_names = paste0("B", sort(unique(clust2))))
  graphics.off()
}

######################3

# reshuffle the labelings 
names(dat_list) <- c("0M", "12M", "3M", "48M", "E120", "E40", "E50", "E70", "E80", "E90")
x <- gene_df$cluster; y <- res_list[[5]]$clustering_res
table(x,y)
res_list[[5]]$clustering_res <- gene_df$cluster

i <- 3
name_intersect <- intersect(res_list[[i]]$gene_name2, res_list[[5]]$gene_name2)
clust1 <- res_list[[i]]$clustering_res[which(res_list[[i]]$gene_name2 %in% name_intersect)]
clust2 <- res_list[[5]]$clustering_res[which(res_list[[5]]$gene_name2 %in% name_intersect)]
order_vec <- c(2, 4, 5, 3, 7, 6, 1)
tmp <- clust1
for(kk in 1:max(clust1)){
  tmp[which(clust1 == kk)] <- order_vec[kk]
}
clust1 <- tmp

png(paste0("../figures/pnas_sensitivity", i, "_cleaned.png"), 
    height = 1500, 
    width = 1500, 
    units = "px", 
    res = 300)
par(mar = c(5,5,6,0.5))
plot_table(clust1, clust2, 
           main = paste0("Comparing results using correlation\nthreshold = ", df_param[i,"cor_threshold"], 
                         " and K = ", df_param[i,"num_clusters"], "\ncompared to baseline"),
           normalize_row = F,
           row_names = paste0("New cluster ",sort(unique(clust1))),
           col_names = paste0("Baseline cluster ", sort(unique(clust2))),
           col_xshift = 1/8+0.01,
           col_offset = -1/16+0.01)
graphics.off()

i <- 7
name_intersect <- intersect(res_list[[i]]$gene_name2, res_list[[5]]$gene_name2)
clust1 <- res_list[[i]]$clustering_res[which(res_list[[i]]$gene_name2 %in% name_intersect)]
clust2 <- res_list[[5]]$clustering_res[which(res_list[[5]]$gene_name2 %in% name_intersect)]
order_vec <- c(4, 6, 2, 8, 9, 3, 7, 5, 1)
tmp <- clust1
for(kk in 1:max(clust1)){
  tmp[which(clust1 == kk)] <- order_vec[kk]
}
clust1 <- tmp

png(paste0("../figures/pnas_sensitivity", i, "_cleaned.png"), 
    height = 1500, 
    width = 1500, 
    units = "px", 
    res = 300)
par(mar = c(4.3,3,4,0.5))
plot_table(clust1, clust2, 
           main = paste0("Comparing results using correlation\nthreshold = ", df_param[i,"cor_threshold"], 
                         " and K = ", df_param[i,"num_clusters"], "\ncompared to baseline"),
           normalize_row = F,
           row_names = paste0("New cluster ",sort(unique(clust1))),
           col_names = paste0("Baseline cluster ", sort(unique(clust2))),
           col_xshift = 1/8+0.03,
           col_offset = 1/8,
           row_offset = 0,
           cex_text = 0.9)
graphics.off()