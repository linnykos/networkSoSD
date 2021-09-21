rm(list=ls())
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
