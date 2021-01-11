
## first do some housekeeping
# manually change the labeling of the clusters
desired_order <- c(6,1,4,3,8,7,5,2)
clustering_res2 <- rep(NA, length(clustering_res))
for(i in 1:length(desired_order)){
  clustering_res2[which(clustering_res == desired_order[i])] <- i
}
table(clustering_res2)
table(clustering_res2)/nrow(adj_list[[1]])

## print out some stats
quantile(power_law)

##########


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

