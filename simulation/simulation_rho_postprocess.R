rm(list=ls())
load("../results/simulation_rho_writeup4.RData")

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

key_mat <- data.frame(method = rownames(res_mat), 
                      full_name = c("SoS-Debias", "SoS-Debias (Weighted)", 
                                    "Sum", "Sum (Weighted)", "SoS", "SoS (Weighted)",
                                    "Tensor matricization", "Tensor matricization (Weighted)",
                                    "Greedy clustering"),
                      color = sapply(c("blue", "blue", "green", "green", "orange", "orange", "lightpurple", "lightpurple", "gray"), 
                                     function(x){ color_vec[which(color_name_vec == x)]}), 
                      lty = c(1,3,1,3,1,3,1,3,1))
key_mat

######################

plot_func <- function(methods, res_mat, key_mat, paramMat, main){
  idx_vec <- sapply(methods, function(x){which(key_mat$method == x)})
  plot(NA, xlim = range(paramMat[,"rho"]), ylim = range(res_mat), xlab = "Sparisity (rho)",
       ylab = "Total clustering accuracy\n(Averaged over 100 trials)", main = main)
  # draw grid
  for(x in paramMat[seq(1, nrow(paramMat), length.out = 8), "rho"]){
    lines(rep(x,2), c(-1e4,1e4), col = "gray", lwd = 0.5, lty = 2)
  }
  for(y in seq(0.5,1,length.out = 11)){
    lines(c(-1e4,1e4), rep(y,2), col = "gray", lwd = 0.5, lty = 2)
  }
  for(i in idx_vec){
    graphics::lines(x = paramMat[,"rho"], y = res_mat[i,], col = key_mat$color[i], 
                    lwd = 2, lty = key_mat$lty[i])
    graphics::points(x = paramMat[,"rho"], y = res_mat[i,], col = key_mat$color[i],  
                     pch = 16, cex = 1)
  }
  legend("bottomright", key_mat$full_name[idx_vec], 
         lty = key_mat$lty[idx_vec], col = key_mat$color[idx_vec], lwd = 2, cex = 0.8)
  
}

######################

png("../figures/simulation_rho.png", height = 1000, width = 2500, units = "px", res = 300)
par(mar = c(5,5,5,0.5), mfrow = c(1,3))

methods <- c("res_ss_debias_F", "res_sum_F", "res_ss_F")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Unweighted")

methods <- c("res_ss_debias_F", "res_ss_debias_T", "res_sum_T", "res_ss_T")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against aggregation\nmethods: Weighted")

methods <- c("res_ss_debias_F", "res_flat_F", "res_flat_T", "res_greedy")
plot_func(methods, res_mat, key_mat, paramMat, main = "Comparison against other\nmethods")
graphics.off()
