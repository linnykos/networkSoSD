rm(list=ls())
library(simulation)
library(networkSoSD)

session_info <- sessionInfo()
date_of_run <- Sys.time()
source_code_info <- readLines("../simulation/simulation_rho.R")
run_suffix <- "_writeup4"

paramMat <- cbind(500, 100, seq(0.025, 0.2, length.out = 15), 0.4, 0.1, 0.5)
colnames(paramMat) <- c("n", "L", "rho", "mem_prop1", "mem_prop2", "mem_prop3")

.l2norm <- function(x){sqrt(sum(x^2))}
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)
K <- 3

trials <- 50
ncores <- 10
doMC::registerDoMC(cores = ncores)
rm(list = c("vec1", "vec2", "vec3", "eigen_mat"))

###############################

rule <- function(vec){
  n <- vec["n"]; L <- vec["L"]; rho <- vec["rho"]
  mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
  
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- networkSoSD::compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- networkSoSD::compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  adj_list <- lapply(1:L, function(i){networkSoSD::generate_adjaceny_mat(prob_list[[i]])})
  
  list(adj_list = adj_list, prob_list = prob_list, membership_vec = membership_vec)
}

########################3

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
color_vec <- color_func(1)[c(9,5,4)]
color_vec2 <- color_func(0.5)[c(9,5,4)]

###########

row_setting <- 1; y <- 1; set.seed(y); dat <- rule(paramMat[row_setting,])
agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
deg_vec <- diag(agg_network)
deg_df <- data.frame(deg = jitter(deg_vec), mem = dat$membership_vec)
dens_list <- lapply(1:3, function(x){
  dens <- stats::density(deg_vec[which(dat$membership_vec == x)])
  dens$y <- dens$y[dens$x >= 0]; dens$x <- dens$x[dens$x >= 0]
  dens
})
xlim <- range(unlist(lapply(dens_list, function(dens){dens$x})))
ylim <- range(c(0,unlist(lapply(dens_list, function(dens){dens$y}))))
max_y <- max(ylim)
dens_list <- lapply(dens_list, function(dens){
  dens$y <- dens$y * 1.1/max_y
  dens
})

png("../figures/simulation_degree_rho_small.png", height = 1500, width = 1500, units = "px", res = 300)
plot(NA, xlim = xlim, ylim = c(0.8,max(dens_list[[1]]$y+3)), main = expression(paste("Degree across clusters for ", rho , "=0.025")),
     ylab = "Community", xlab = "Total degree across all layers", yaxt='n')
axis(2, at = 1:3, labels = c("3","2","1"))
for(i in 1:3){
  polygon(x = c(dens_list[[i]]$x[c(1,1:length(dens_list[[i]]$x),length(dens_list[[i]]$x))]),
          y = c(0, dens_list[[i]]$y, 0)+(4-i), col = color_vec2[i])
  
  tmp_df <- deg_df[deg_df$mem == i,]
  tmp <- boxplot(deg ~ mem, data = tmp_df,
          horizontal=TRUE, col = color_vec[i], range = ifelse(row_setting == 1, 1.5, 2.1), add = T, 
          pars = list(outpch = 16), at = 4-i)
}
graphics.off()

####################


row_setting <- 15; y <- 1; set.seed(y); dat <- rule(paramMat[row_setting,])
agg_network <- networkSoSD::aggregate_networks(dat$adj_list, method = "ss")
deg_vec <- diag(agg_network)
deg_df <- data.frame(deg = jitter(deg_vec), mem = dat$membership_vec)
dens_list <- lapply(1:3, function(x){
  dens <- stats::density(deg_vec[which(dat$membership_vec == x)])
  dens$y <- dens$y[dens$x >= 0]; dens$x <- dens$x[dens$x >= 0]
  dens
})
xlim <- range(unlist(lapply(dens_list, function(dens){dens$x})))
ylim <- range(c(0,unlist(lapply(dens_list, function(dens){dens$y}))))
max_y <- max(ylim)
dens_list <- lapply(dens_list, function(dens){
  dens$y <- dens$y * 1.1/max_y
  dens
})

png("../figures/simulation_degree_rho_large.png", height = 1500, width = 1500, units = "px", res = 300)
plot(NA, xlim = xlim, ylim = c(0.8,max(dens_list[[1]]$y+3)), main = expression(paste("Degree across clusters for ", rho , "=0.2")),
     ylab = "Community", xlab = "Total degree across all layers", yaxt='n')
axis(2, at = 1:3, labels = c("3","2","1"))
for(i in 1:3){
  polygon(x = c(dens_list[[i]]$x[c(1,1:length(dens_list[[i]]$x),length(dens_list[[i]]$x))]),
          y = c(0, dens_list[[i]]$y, 0)+(4-i), col = color_vec2[i])
  
  tmp_df <- deg_df[deg_df$mem == i,]
  tmp <- boxplot(deg ~ mem, data = tmp_df,
                 horizontal=TRUE, col = color_vec[i], range = ifelse(row_setting == 1, 1.5, 2.1), add = T, 
                 pars = list(outpch = 16), at = 4-i)
}
graphics.off()
