rm(list=ls())
load("../results/simulation_rho_writeup4.RData")
source("../simulation/simulation_key.R")

vec <- paramMat[1,]; n <- vec["n"]
mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))

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
                      full_name = c("SoS-Debias", "SoS-Debias (Weighted)", 
                                    "Sum", "Sum (Weighted)", "SoS", "SoS (Weighted)",
                                    "Tensor matricization", "Tensor matricization (Weighted)",
                                    "Greedy clustering"),
                      color = sapply(c("blue", "blue", "green", "green", "orange", "orange", "lightpurple", "lightpurple", "gray"), 
                                     function(x){ color_vec[which(color_name_vec == x)]}), 
                      lty = c(1,3,1,3,1,3,1,3,1), pch = c(16,16, 15,15, 17,17, 18,18, 16))
key_mat

######################

png("../figures/simulation_rho.png", height = 1200, width = 1500, units = "px", res = 300)
par(mar = c(5,5,5,0.5))
methods <- c( "res_sum_F", "res_ss_F", "res_ss_debias_F", "res_flat_F")
plot_func(methods, res_mat, key_mat, paramMat, 
          main = "Comparison against aggregation\nmethods (Three communites)",
          y_seq = seq(0, 0.4, length.out = 9), legend_cex = 0.7, legend_loc = "bottomleft")
graphics.off()

png("../figures/simulation_rho_weighted.png", height = 1200, width = 1500, units = "px", res = 300)
par(mar = c(5,5,5,0.5))
methods <- c( "res_sum_T", "res_ss_T", "res_ss_debias_T", "res_flat_T")
plot_func(methods, res_mat, key_mat, paramMat, 
          main = "Comparison against aggregation\nmethods (Weighted, Three communites)",
          y_seq = seq(0, 0.4, length.out = 9), legend_cex = 0.55, legend_loc = "topright")
graphics.off()

# png("../figures/simulation_rho.png", height = 1000, width = 2500, units = "px", res = 300)
# par(mar = c(5,5,5,0.5), mfrow = c(1,3))
# 
# methods <- c("res_ss_debias_F", "res_sum_F", "res_ss_F")
# plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Unweighted")
# 
# methods <- c("res_ss_debias_F", "res_ss_debias_T", "res_sum_T", "res_ss_T")
# plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Weighted")
# 
# methods <- c("res_ss_debias_F", "res_flat_F", "res_flat_T", "res_greedy")
# plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against other\nmethods")
# graphics.off()
# 
# ######################################################
# ######################################################
# 
# rm(list=ls())
# load("../results/simulation_rho_same.RData")
# source("../simulation/simulation_key.R")
# 
# vec <- paramMat[1,]; n <- vec["n"]
# mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
# membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
# if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
# 
# res_mat <- sapply(1:length(res), function(i){
#   tmp <- sapply(res[[i]], function(res_list){
#     sapply(res_list, function(vec){
#       vec <- align_two_membership_vectors(membership_vec, vec, override = T)
#       tab <- table(membership_vec, vec)
#       1-sum(diag(tab))/sum(tab)
#     })
#   })
#   
#   apply(tmp, 1, mean)
# })
# 
# key_mat <- data.frame(method = rownames(res_mat)[c(1:8,11)], 
#                       full_name = c("Bias-adjusted SoS", "Bias-adjusted SoS (Weighted)", 
#                                     "Sum", "Sum (Weighted)", "SoS", "SoS (Weighted)",
#                                     "Tensor matricization", "Tensor matricization (Weighted)",
#                                     "Greedy clustering"),
#                       color = sapply(c("blue", "blue", "green", "green", "orange", "orange", "lightpurple", "lightpurple", "gray"), 
#                                      function(x){ color_vec[which(color_name_vec == x)]}), 
#                       lty = c(1,3,1,3,1,3,1,3,1))
# key_mat

######################

# png("../figures/simulation_rho_same.png", height = 1000, width = 2500, units = "px", res = 300)
# par(mar = c(5,5,5,0.5), mfrow = c(1,3))
# 
# methods <- c("res_ss_debias_F", "res_sum_F", "res_ss_F")
# plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Unweighted")
# 
# methods <- c("res_ss_debias_F", "res_ss_debias_T", "res_sum_T", "res_ss_T")
# plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Weighted")
# 
# methods <- c("res_ss_debias_F", "res_flat_F", "res_flat_T", "res_greedy")
# plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against other\nmethods")
# graphics.off()


