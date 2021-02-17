rm(list=ls())
load("../results/simulation_rho_simple.RData")
source("../simulation/simulation_key.R")

vec <- paramMat[1,]; n <- vec["n"]
mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]
membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n))
if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(2, n-length(membership_vec)))

res_mat <- sapply(1:length(res), function(i){
  tmp <- sapply(res[[i]], function(res_list){
    sapply(res_list, function(vec){
      vec <- align_two_membership_vectors(membership_vec, vec)
      tab <- table(membership_vec, vec)
      1-sum(diag(tab))/sum(tab)
    })
  })
  
  apply(tmp, 1, mean)
})

key_mat <- data.frame(method = rownames(res_mat), 
                      full_name = c("Bias-adjusted SoS", "Sum", "SoS", 
                                    "Tensor matricization"),
                      color = sapply(c("blue", "green", "orange", "lightpurple"), 
                                     function(x){ color_vec[which(color_name_vec == x)]}), 
                      lty = c(1,1,1,1), pch = c(16,15,17,18))
key_mat

png("../figures/simulation_rho_simple.png", height = 1200, width = 1500, units = "px", res = 300)
par(mar = c(5,5,5,0.5))
methods <- c( "res_sum_F", "res_ss_F", "res_ss_debias_F")
plot_func(methods, res_mat, key_mat, paramMat, 
          main = "Comparison against aggregation\nmethods (Two communities)",
          y_seq = seq(0, 0.4, length.out = 9), legend_cex = 0.7, legend_loc = "bottomleft")
graphics.off()
