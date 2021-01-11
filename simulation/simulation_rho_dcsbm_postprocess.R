rm(list=ls())
load("../results/simulation_rho_dcsbm.RData")
source("../simulation/simulation_key.R")

vec <- paramMat[1,]; n <- vec["n"]
mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))

res_mat <- sapply(1:length(res), function(i){
  tmp <- sapply(res[[i]], function(res_list){
    sapply(res_list, function(vec){
      vec <- align_two_membership_vectors(membership_vec, vec)
      tab <- table(membership_vec, vec)
      sum(diag(tab))/sum(tab)
    })
  })
  
  apply(tmp, 1, mean)
})


key_mat <- data.frame(method = rownames(res_mat), 
                      full_name = c("SoS-Debias", "SoS-Debias (Weighted)", 
                                    "Sum", "Sum (Weighted)", "SoS", "SoS (Weighted)",
                                    "Tensor matricization", "Tensor matricization (Weighted)"),
                      color = sapply(c("blue", "blue", "green", "green", "orange", "orange", 
                                       "lightpurple", "lightpurple"), 
                                     function(x){ color_vec[which(color_name_vec == x)]}), 
                      lty = c(1,3,1,3,1,3,1,3))
key_mat

######################

png("../figures/simulation_rho_dcsbm.png", height = 1000, width = 2500, units = "px", res = 300)
par(mar = c(5,5,5,0.5), mfrow = c(1,3))

methods <- c("res_ss_debias_F", "res_sum_F", "res_ss_F")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Unweighted")

methods <- c("res_ss_debias_F", "res_ss_debias_T", "res_sum_T", "res_ss_T")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Weighted")

methods <- c("res_ss_debias_F", "res_flat_F", "res_flat_T")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against other\nmethods")
graphics.off()

###################################3

rm(list=ls())
load("../results/simulation_rho_dcsbm_equaldeg.RData")
source("../simulation/simulation_key.R")

vec <- paramMat[1,]; n <- vec["n"]
mem_prop1 <- vec["mem_prop1"]; mem_prop2 <- vec["mem_prop2"]; mem_prop3 <- vec["mem_prop3"]
membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))

res_mat <- sapply(1:length(res), function(i){
  tmp <- sapply(res[[i]], function(res_list){
    sapply(res_list, function(vec){
      vec <- align_two_membership_vectors(membership_vec, vec)
      tab <- table(membership_vec, vec)
      sum(diag(tab))/sum(tab)
    })
  })
  
  apply(tmp, 1, mean)
})

key_mat <- data.frame(method = rownames(res_mat)[-c(9:10)], 
                      full_name = c("SoS-Debias", "SoS-Debias (Weighted)", 
                                    "Sum", "Sum (Weighted)", "SoS", "SoS (Weighted)",
                                    "Tensor matricization", "Tensor matricization (Weighted)",
                                    "SoS-Debias SBM", "SoS SBM"),
                      color = sapply(c("blue", "blue", "green", "green", "orange", "orange", 
                                       "lightpurple", "lightpurple", "skyblue", "gray"), 
                                     function(x){ color_vec[which(color_name_vec == x)]}), 
                      lty = c(1,3,1,3,1,3,1,3,1,1))
key_mat

######################

png("../figures/simulation_rho_dcsbm_equaldeg.png", height = 1000, width = 2500, units = "px", res = 300)
par(mar = c(5,5,5,0.5), mfrow = c(1,3))

methods <- c("res_ss_debias_F", "res_sum_F", "res_ss_F", "res_ss_debias_T_sbm", "res_ss_T_sbm")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Unweighted")

methods <- c("res_ss_debias_F", "res_ss_debias_T", "res_sum_T", "res_ss_T")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Weighted")

methods <- c("res_ss_debias_F", "res_flat_F", "res_flat_T")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against other\nmethods")
graphics.off()
