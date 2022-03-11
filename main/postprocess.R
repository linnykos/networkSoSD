rm(list=ls()); set.seed(10)
load("../results/main_analysis_revision.RData")
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

## first do some rearranging
# manually change the labeling of the clusters

# first compute all the average matrices according to the clustering_res,
# so we know which genes fall into which clustering in the correct ordering
names(adj_list) <- c("0M", "12M", "3M", "48M", "E120", "E40", "E50", "E70", "E80", "E90")

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
  print(round(diag(cor_list[[i]]),2))
  print("===")
}  
print(table(clustering_res))

order_vec <- c(3, 7, 1, 2, 6, 8, 4, 5)
clustering_res2 <- rep(NA, length(clustering_res))
for(i in 1:length(order_vec)){
  clustering_res2[which(clustering_res == i)] <- order_vec[i]
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

# determine how many clusters the naively-summed graph has
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

####################

time_stamp <- c("0M", "12M", "3M", "48M", "E120", "E40", "E50", "E70", "E80", "E90")
time_order <- c(7,9,8,10,6,1,2,3,4,5) # HARD CODED

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


##################

cor_list <- lapply(adj_list, function(adj_mat){
  K <- max(clustering_res)
  cor_mat <- matrix(0, K, K)
  for(i in 1:K){
    idx1 <- which(clustering_res2 == i)
    for(j in i:K){
      idx2 <- which(clustering_res2 == j)
      tmp <- adj_mat[idx1,idx2]
      cor_mat[i,j] <- sum(adj_mat[idx1,idx2])/2
      cor_mat[j,i] <- cor_mat[i,j]
    }
  }
  cor_mat
})
names(cor_list) <- names(adj_list)

within_cluster_mat <- sapply(cor_list, function(cor_mat){
  num <- table(clustering_res2)
  diag(cor_mat)/(num * (num-1)/2)
})
colnames(within_cluster_mat) <- names(cor_list)
between_cluster_mat <- sapply(cor_list, function(cor_mat){
  offdiag_vec <- (colSums(cor_mat) - as.numeric(diag(cor_mat)))
  num <- table(clustering_res2)
  total_entries <- sapply(1:length(num), function(i){
    num[i]*sum(num[-i])
  })
  offdiag_vec/total_entries
})
colnames(between_cluster_mat) <- names(cor_list)

col_palette <- color_vec[c(7, 5, 8, 1, 3, 2, 4, 6)]
lty_vec <- c(1,5,3,1,5,3,1,5)
time_order <- c(6:10, 5, 1:4)
png(paste0("../figures/pnas_connectivity.png"), 
    height = 1500, width = 2500, units = "px", res = 300)
ylim <- c(0, max(within_cluster_mat))
plot(NA, xlim = c(1,ncol(within_cluster_mat)), 
     ylim = ylim,
     main = "Within-cluster connectivity over time",
     xlab = "Time",
     ylab = "Edge density",
     xaxt = "n")
for(i in 1:nrow(within_cluster_mat)){
  lines(x = 1:ncol(within_cluster_mat),
        y = within_cluster_mat[i,time_order],
        col = "white",
        lwd = 7)
  # lines(x = 1:ncol(within_cluster_mat),
  #       y = within_cluster_mat[i,time_order],
  #       col = "black",
  #       lwd = 2.1, lty = lty_vec[i])
  lines(x = 1:ncol(within_cluster_mat),
        y = within_cluster_mat[i,time_order],
        col = col_palette[i],
        lwd = 5, lty = lty_vec[i])
}
axis(1, at = 1:ncol(within_cluster_mat), 
     labels = colnames(within_cluster_mat)[time_order])
legend("bottomright", paste0("Cluster ", 1:8), 
       col=col_palette, 
       lty=c(1,2,3,1,2,3,1,2), lwd = 5)
graphics.off()
