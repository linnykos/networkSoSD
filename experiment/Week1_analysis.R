rm(list=ls())
load("pnas.RData")

for(i in 1:length(dat_list)){
  print(i)
  num_centers <- 500
  kmean_res <- stats::kmeans(dat_list[[i]], centers = num_centers)
  
  png(paste0("Writeup1_pnas_cor", i, "_magnitude.png"), height = 3000, width = 3000, units = "px", res = 300)
  plot(NA, main = i, xlim = c(0, nrow(dat_list[[i]])), ylim = c(0,1))
  for(j in 1:num_centers){
    lines(sort(abs(kmean_res$centers[j,])), col = rgb(0.5,0.5,0.5,0.1))
  }
  graphics.off()
}

adj_list <- lapply(dat_list, function(dat_list){
  adj_mat <- matrix(0, nrow = nrow(dat_list), ncol = ncol(dat_list))
  adj_mat[which(dat_list >= 0.15)] <- 1
  adj_mat
})
sapply(adj_list, function(adj_mat){sum(adj_mat)/(2*prod(dim(adj_mat)))})

adj_array <- array(0, dim = c(nrow(adj_list[[1]]), ncol(adj_list[[1]]), length(adj_list)))
for(i in 1:length(adj_list)){
  print(i)
  tmp <- adj_list[[i]]
  tmp <- tcrossprod(tmp) - diag(apply(tmp, 1, sum))
  adj_array[,,i] <- tmp
}
tot_adj_mat <- apply(adj_array, c(1,2), sum)
svd_res <- RSpectra::svds(tot_adj_mat, k = 50)

png(paste0("Writeup1_pnas_svd.png"), height = 1500, width = 2500, units = "px", res = 300)
par(mfrow = c(1,2))
plot(svd_res$d, pch = 16, xlab = "Index", ylab = "Eigenvalue magnitude")
plot(log(svd_res$d), pch = 16, xlab = "Index", ylab = "Log of eigenvalue magnitude")
graphics.off()

set.seed(10)
kmean_res <- stats::kmeans(svd_res$u[,1:20], centers = 20)

table(kmean_res$cluster)

gene_idx <- unlist(lapply(1:20, function(x){which(kmean_res$cluster == x)}))
clockwise90 = function(a) { t(a[nrow(a):1,]) } 
for(i in 1:length(dat_list)){
  print(i)
  png(paste0("Writeup1_pnas_adj", i, ".png"), height = 3000, width = 3000, units = "px", res = 300)
  image(clockwise90(adj_list[[i]][gene_idx, gene_idx]), main = i)
  
  for(j in 1:19){
    len <- length(which(kmean_res$cluster <= j))/length(kmean_res$cluster)
    lines(c(-1e4,1e4), rep(1-len,2), lwd = 2, lty = 2)
    lines(rep(len,2), c(-1e4,1e4), lwd = 2, lty = 2)
  }
  
  graphics.off()
}
