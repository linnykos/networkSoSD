rm(list=ls())
load("../results/simulation_rho.RData")
source("../simulation/simulation_key.R")

vec <- df_param[1,]; n <- as.numeric(vec["n"])
mem_prop1 <- as.numeric(vec["mem_prop1"])
mem_prop2 <- as.numeric(vec["mem_prop2"])
mem_prop3 <- as.numeric(vec["mem_prop3"])
membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))

res_mat <- sapply(1:length(res), function(i){
  tmp <- sapply(1:length(res[[i]]), function(trial){
    res_list <- res[[i]][[trial]]
    tmp <- res_list$result
    tmp <- tmp[-grep("time", names(tmp))]
    
    sapply(tmp, function(vec){
      vec <- networkSoSD:::align_two_membership_vectors(membership_vec, vec)
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
                                    "LMF",
                                    "Co-reg"),
                      color = sapply(c("blue", "blue", 
                                       "green", "green", 
                                       "orange", "orange", 
                                       "lightpurple", "lightpurple", 
                                       "red", "gray"), 
                                     function(x){ color_vec[which(color_name_vec == x)]}), 
                      lty = c(1,3,1,3,1,3,1,3,1,1), 
                      pch = c(16,16, 15,15, 17,17, 18,18,16,15))
key_mat

######################

png("../figures/simulation_rho.png", height = 1200, width = 1500, units = "px", res = 300)
par(mar = c(5,5,5,0.5))
methods <- c( "res_sum_F", "res_ss_F", "res_ss_debias_F", "res_flat_F", 
              "chen_linked", "chen_coreg")
plot_func(methods, res_mat, key_mat, df_param, 
          main = "Comparison against aggregation\nmethods (Three communites)",
          y_seq = seq(0, 0.5, length.out = 11), 
          legend_cex = 0.55, 
          legend_loc = "bottomleft")
graphics.off()


png("../figures/simulation_rho_nolegend.png", height = 1200, width = 1500, units = "px", res = 300)
par(mar = c(5,5,5,0.5))
methods <- c( "res_sum_F", "res_ss_F", "res_ss_debias_F", "res_flat_F", 
              "chen_linked", "chen_coreg")
plot_func(methods, res_mat, key_mat, df_param, 
          main = "Comparison against aggregation\nmethods (Three communites)",
          y_seq = seq(0, 0.5, length.out = 11), 
          include_legend = F,
          legend_cex = NA, 
          legend_loc = NA)
graphics.off()


png("../figures/simulation_rho_weighted.png", height = 1200, width = 1500, units = "px", res = 300)
par(mar = c(5,5,5,0.5))
methods <- c( "res_sum_T", "res_ss_T", "res_ss_debias_T", "res_flat_T")
plot_func(methods, res_mat, key_mat, df_param, 
          main = "Comparison against aggregation\nmethods (Weighted, Three communites)",
          y_seq = seq(0, 0.4, length.out = 9), 
          legend_cex = 0.55, 
          legend_loc = "topright")
graphics.off()

##############

# now let's see what the simulation times are
time_list <- lapply(res, function(trial_list){
  sapply(trial_list, function(trial_res){
    name_vec <- names(trial_res$result)
    idx <- grep("time_*", name_vec)
    sapply(trial_res$result[idx], function(val){
      if(attributes(val)$units != "secs"){
        as.numeric(val)*60
      } else {
        as.numeric(val)
      }
    })
  })
})

time_mat <- do.call(cbind, time_list)
apply(time_mat, 1, mean)