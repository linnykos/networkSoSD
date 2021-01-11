rm(list=ls())
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()

paramMat <- cbind(500, 100, seq(0.025, 0.2, length.out = 15), 0.4, 0.1, 0.5)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)
# round(B1, 2); round(B2, 2)
K <- 3

rm(list = c("vec1", "vec2", "vec3", "eigen_mat"))

rule <- function(vec){
  n <- vec["n"]; L <- vec["L"]; rho <- vec["rho"]
  mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list, prob_list = prob_list)
}

############

# first, measure the amount of bias caused by the diagonal entries
bias_vec <- sapply(1:nrow(paramMat), function(i){
  y <- 1; set.seed(y); dat <- rule(paramMat[i,])
  
  mat1 <- dat$prob_list[[1]]; mat2 <- dat$prob_list[[paramMat[i,"L"]/2+1]]
  sq_mat <- crossprod(mat1) + crossprod(mat2)
  diag_vec <- colSums(mat1 + mat2)
  
  sq_mat2 <- sq_mat + diag(diag_vec)
  
  # eigen1 <- eigen(sq_mat)$vectors; eigen2 <- eigen(sq_mat2)$vectors; par(mfrow = c(1,2)); image(t(eigen1[,1:3])); image(t(eigen2[,1:3]))
  
  # eigen1 <- eigen(sq_mat)$vectors[,1:3]
  # eigen2 <- eigen(sq_mat2)$vectors[,-c(1:3)]
  # svd(crossprod(eigen1, eigen2))$d[1]
  
  vec1 <- eigen(sq_mat)$values; vec2 <- eigen(sq_mat2)$values
  ratio1 <- (vec1[K]-vec1[K+1])/vec1[K+1]; ratio2 <- (vec2[K]-vec2[K+1])/vec2[K+1] 
  ratio2 # note: ratio1 is always high, since vec1[K+1] = 0
})

png("../figures/simulation_bias.png", height = 1200, width = 1400, units = "px", res = 300)
plot(NA, xlim = range(paramMat[,"rho"]), ylim = range(bias_vec), 
     xlab = "Sparisity (rho)", ylab = "Eigengap (ratio)",
     main = "Eigengap induced by diagonal bias")
for(x in paramMat[seq(1, nrow(paramMat), by = 2), "rho"]){
  lines(rep(x,2), c(-1e4,1e4), col = "gray", lwd = 0.5, lty = 2)
}
for(y in seq(0,0.1,length.out = 6)){
  lines(c(-1e4,1e4), rep(y,2), col = "gray", lwd = 0.5, lty = 2)
}
points(paramMat[,"rho"], bias_vec, pch = 16)
graphics.off()

#####################

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
color_vec2 <- color_func(0.1)

i <- 10; y <- 29; set.seed(y); dat <- rule(paramMat[i,])
agg_network1 <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
agg_network2 <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
eigen1 <- .svd_projection(agg_network1, K = 3, weighted = F)
eigen2 <- .svd_projection(agg_network2, K = 3, weighted = F)
est_clust1 <- stats::kmeans(eigen1, 3)$cluster
est_clust2 <- stats::kmeans(eigen2, 3)$cluster
# plot(eigen2[,1], eigen2[,2], asp = T, col = est_clust2)

agg_network <- networkSoSD::aggregate_networks(dat$prob_list, method = "ss")
eigen_true <- .svd_projection(agg_network1, K = 3, weighted = F)
# par(mfrow = c(1,2)); image(t(eigen1[,1:3])); image(t(eigen2[,1:3]))

#########

# reparameterize
svd_tmp <- svd(crossprod(eigen2, eigen1))
eigen2 <- eigen2 %*% svd_tmp$u %*% t(svd_tmp$v)

#########

vec <- paramMat[i,]
n <- vec["n"]; L <- vec["L"]; rho <- vec["rho"]
mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
idx <- c(1,which(diff(membership_vec) > 0)+1)
breaks <- seq(0.5, max(membership_vec)+1, by = 1)

png("../figures/simulation_3d.png", height = 2500, width = 2500, units = "px", res = 300)
par(mfrow = c(2,2), mar = c(1, 0.5, 1, 0.5))
max_diff <- max(apply(rbind(eigen1, eigen2), 2, function(x){diff(range(x))/2}))
max_diff <- max_diff/2.5
lim_list <- lapply(1:3, function(x){
  mean(c(eigen1[,x], eigen2[,x])) + c(-1,1)*max_diff
})
plot3D::scatter3D(x = eigen2[,1], y = eigen2[,2], z = eigen2[,3],
                  surface = FALSE, colvar = c(3,2,1)[est_clust2],
                  cex = 0.5,
                  breaks = breaks, col = color_vec[c(1,2,5)], colkey = F, pch = 16,
                  main = "SoS (Estimated clusters)", 
                  xlim = lim_list[[1]], ylim = lim_list[[2]], zlim = lim_list[[3]],
                  phi = 225, theta = 190, 
                  xlab = "Eigenvector 1", ylab = "Eigenvector 2", zlab = "Eigenvector 3")
plot3D::scatter3D(x = eigen2[,1], y = eigen2[,2], z = eigen2[,3],
                  surface = FALSE, colvar = membership_vec,
                  cex = 0.5,
                  breaks = breaks, col = color_vec[c(1,2,5)], colkey = F, pch = 16,
                  main = "SoS (True clusters)", 
                  xlim = lim_list[[1]], ylim = lim_list[[2]], zlim = lim_list[[3]],
                  phi = 225, theta = 190, 
                  xlab = "Eigenvector 1", ylab = "Eigenvector 2", zlab = "Eigenvector 3")

plot3D::scatter3D(x = eigen1[,1], y = eigen1[,2], z = eigen1[,3],
                  surface = FALSE, colvar = c(3,1,2)[est_clust1],
                  cex = 0.5,
                  breaks = breaks, col = color_vec[c(1,2,5)], colkey = F, pch = 16,
                  main = "SoS-Debias (Estimated clusters)", 
                  xlim = lim_list[[1]], ylim = lim_list[[2]], zlim = lim_list[[3]],
                  phi = 225, theta = 190, 
                  xlab = "Eigenvector 1", ylab = "Eigenvector 2", zlab = "Eigenvector 3")
plot3D::scatter3D(x = eigen1[,1], y = eigen1[,2], z = eigen1[,3],
                  surface = FALSE, colvar = membership_vec,
                  cex = 0.5,
                  breaks = breaks, col = color_vec[c(1,2,5)], colkey = F, pch = 16,
                  main = "SoS-Debias (True clusters)", 
                  xlim = lim_list[[1]], ylim = lim_list[[2]], zlim = lim_list[[3]],
                  phi = 225, theta = 190, 
                  xlab = "Eigenvector 1", ylab = "Eigenvector 2", zlab = "Eigenvector 3")
graphics.off()
