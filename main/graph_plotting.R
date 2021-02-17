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
gene_idx <- unlist(lapply(c(1,4,5,6), function(x){which(clustering_res2 == x)}))
clustering_res2 <- clustering_res2[gene_idx]

tmp <- networkSoSD::aggregate_networks(adj_list, method = 'sum')
tmp <- tmp[gene_idx, gene_idx]
set.seed(10)
tmp_svd <- networkSoSD:::.svd_truncated(tmp, K = K)
set.seed(10)
umap_embedding3 <- Seurat::RunUMAP(networkSoSD:::.mult_mat_vec(tmp_svd$u, tmp_svd$d), verbose = F)@cell.embeddings
col_umap <- color_vec[c(9,2:8)][clustering_res2]

#####################
time_stamp <- c("0M", "12M", "3M", "48M", "E120", "E40", "E50", "E70", "E80", "E90")
time_order <- c(7,9,8,10,6,1,2,3,4,5)

for(i in 1:length(adj_list)){
  print(i)
  
  set.seed(10)
  tmp <- adj_list[[i]][gene_idx, gene_idx]
  tmp <- tmp * matrix(stats::rbinom(prod(dim(tmp)), size = 1, prob = 0.1), nrow = nrow(tmp), ncol = ncol(tmp))
  
  g <- igraph::graph_from_adjacency_matrix(tmp, mode = 'undirected')
  igraph::V(g)$color <- color_vec[8]
  vertex_size <- rep(2, length(gene_idx))
  
  edge_mat <- igraph::get.edgelist(g)
  edge_width <- sapply(1:nrow(edge_mat), function(i){
    ifelse(clustering_res2[edge_mat[i,1]] == clustering_res2[edge_mat[i,2]], 1, 0.5)
  })
  edge_color <- sapply(1:nrow(edge_mat), function(i){
    if(clustering_res2[edge_mat[i,1]] == clustering_res2[edge_mat[i,2]]){
      grDevices::rgb(0.5, 0.5, 0.5, 0.5)
    } else {
      grDevices::rgb(0.5, 0.5, 0.5, 0.1)
    }
  })
  
  png(paste0("../figures/graph", time_order[i], ".png"), height = 1200, width = 1200, units = "px", res = 300)
  par(mar = rep(0.5, 4))
  graphics::plot(g, layout = umap_embedding3, vertex.label = NA, 
                 vertex.size = vertex_size,
                 edge.color = grDevices::rgb(0.5, 0.5, 0.5, 0.1), curved = T,
                 edge.width = edge_width)
  graphics.off()
}

idx_odd <- seq(1, length(gene_idx), by = 8)
clockwise90 <- function(a) { t(a[nrow(a):1,]) } 
for(i in 1:length(adj_list)){
  print(i)
  tmp <- adj_list[[i]][gene_idx[idx_odd], gene_idx[idx_odd]]
  diag(tmp) <- -1
  
  png(paste0("../figures/pnas_subset_adj", time_order[i], ".png"), height = 1500, width = 1300, units = "px", res = 300)
  par(mar = c(2, 2, 4, 0.6))
  image(clockwise90(tmp), main = paste0("Sub-network for time ", time_stamp[i]), 
        col = c("white", hcl.colors(2, "Cividis")),
        breaks = c(-1.5, -.5, .5, 1.5))
  
  graphics.off()
}


