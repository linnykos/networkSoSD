rm(list=ls())
library(networkSoSD)

session_info <- devtools::session_info()
date_of_run <- Sys.time()

df_param <- cbind(3, 500, 100, seq(0.025, 0.2, length.out = 15), 0.4, 0.1, 0.5)
colnames(df_param) <- c("K", "n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")
df_param <- as.data.frame(df_param)

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
bias_vec <- sapply(1:nrow(df_param), function(i){
  y <- 1; set.seed(y); dat <- rule(df_param[i,])
  
  mat1 <- dat$prob_list[[1]]; mat2 <- dat$prob_list[[df_param[i,"L"]/2+1]]
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

set.seed(10)
row_setting <- 11; y <- 1; set.seed(y); dat <- rule(df_param[row_setting,])
agg_mat1 <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
eigen_val1 <- eigen(agg_mat1)$values

agg_mat2 <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss_debias")
eigen_val2 <- eigen(agg_mat2)$values

png("../figures/simulation_bias.png", height = 900, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(5,5,5,0.5))

hist(eigen_val1[-1], main = "Est. eigenvalues (excluding first)\nfor SoS aggregation", col = "gray",
     breaks = 50, xlab = "Estimate eigenvalues")
for(i in 2:3){
  lines(rep(eigen_val1[i], 2), c(0, 500), col = grDevices::rgb(189/255, 57/255, 60/255), lwd = 2, lty = 2)
}

hist(eigen_val2[-1], main = "Est. eigenvalues (excluding first)\nfor Bias-adjusted SoS aggregation", col = "gray",
     breaks = 50, xlab = "Estimated eigenvalues")
for(i in 2:3){
  lines(rep(eigen_val2[i], 2), c(0, 500), col = grDevices::rgb(189/255, 57/255, 60/255), lwd = 2, lty = 2)
}

plot(NA, xlim = range(df_param[,"rho"]), ylim = range(bias_vec), 
     xlab = expression(paste("Sparisity (", rho, ")")), ylab = "Population eigengap (ratio)",
     main = "Population eigengap induced by\ndiagonal bias")
for(x in df_param[seq(1, nrow(df_param), by = 2), "rho"]){
  lines(rep(x,2), c(-1e4,1e4), col = "gray", lwd = 0.5, lty = 2)
}
for(y in seq(0,0.1,length.out = 6)){
  lines(c(-1e4,1e4), rep(y,2), col = "gray", lwd = 0.5, lty = 2)
}
points(df_param[,"rho"], bias_vec, pch = 16)

graphics.off()

#####################################


