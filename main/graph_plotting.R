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
desired_order <- c(6,1,4,3,8,7,5,2)
clustering_res2 <- rep(NA, length(clustering_res))
for(i in 1:length(desired_order)){
  clustering_res2[which(clustering_res == desired_order[i])] <- i
}

tmp <- networkSoSD::aggregate_networks(adj_list, method = 'sum')
set.seed(10)
tmp_svd <- networkSoSD:::.svd_truncated(tmp, K = K)
set.seed(10)
umap_embedding3 <- Seurat::RunUMAP(networkSoSD:::.mult_mat_vec(tmp_svd$u, tmp_svd$d), verbose = F)@cell.embeddings
col_umap <- color_vec[c(9,2:8)][clustering_res2]

#####################

for(i in 1:length(adj_list)){
  print(i)
  
  tmp <- adj_list[[i]]
  tmp <- tmp * matrix(stats::rbinom(prod(dim(tmp)), size = 1, prob = 0.01), nrow = nrow(tmp), ncol = ncol(tmp))
  
  g <- igraph::graph_from_adjacency_matrix(tmp, mode = 'undirected')
  igraph::V(g)$color <- col_umap
  vertex_size <- rep(2, nrow(adj_list[[1]]))
  vertex_size[clustering_res2 == 8] <- 1
  
  png(paste0("../figures/Writeup4_graph", i, ".png"), height = 1500, width = 1500, units = "px", res = 300)
  par(mar = rep(0.5, 4))
  graphics::plot(g, layout = umap_embedding3, vertex.label = NA, 
                 vertex.size = vertex_size,
                 edge.color = grDevices::rgb(0.5, 0.5, 0.5, 0.1), curved = T,
                 edge.width = 0.5)
  graphics.off()
}


