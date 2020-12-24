rm(list=ls())
load("../results/simulation_rho.RData")

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

png("../figures/simulation_rho.png", height = 1500, width = 2000, units = "px", res = 300)
lty_vec <- c(1,1,1,2,1); lwd_vec <- c(2,2,2,1,2)
par(mar = c(5,5,5,1))
plot(NA, xlim = range(paramMat[,"rho"]), ylim = range(res_mat), xlab = "Sparisity (rho)",
     ylab = "Total clustering accuracy\n(Averaged over 100 trials)", main = "Clustering accuracy across different methods")
for(i in 1:nrow(res_mat)){
  graphics::lines(x = paramMat[,"rho"], y = res_mat[i,], col = i, 
                  lwd = lwd_vec[i], lty = lty_vec[i])
  graphics::points(x = paramMat[,"rho"], y = res_mat[i,], col = i, 
                  pch = 16, cex = lwd_vec[i])
}

legend("topleft", c("SoS-Debias", "Sum (Unweighted)", "SoS (Unweighted)",
                    "Sum (Weighted)", "SoS (Weighted)"),  
       bty="n", fill=1:5)
graphics.off()
