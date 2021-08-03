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


plot_func <- function(methods, res_mat, key_mat, df_param, main, y_seq = seq(0, 1,length.out = 11),
                      legend_cex = 1, legend_loc = 'bottomright'){
  idx_vec <- sapply(methods, function(x){which(key_mat$method == x)})
  idx_vec2 <- sapply(methods, function(x){which(rownames(res_mat) == x)}) 
  stopifnot(length(idx_vec) == length(idx_vec2))
  
  plot(NA, xlim = range(as.numeric(df_param[,"rho"])), ylim = range(res_mat), xlab = expression(paste("Sparisity (", rho, ")")),
       ylab = "% mis-clustered nodes\n(Averaged over 100 trials)", main = main)
  # draw grid
  for(x in as.numeric(df_param[seq(1, nrow(df_param), by = 2), "rho"])){
    lines(rep(x,2), c(-1e4,1e4), col = "gray", lwd = 0.5, lty = 2)
  }
  for(y in y_seq){
    lines(c(-1e4,1e4), rep(y,2), col = "gray", lwd = 0.5, lty = 2)
  }
  for(i in 1:length(idx_vec)){
    graphics::lines(x = as.numeric(df_param[,"rho"]), y = res_mat[idx_vec2[i],], col = key_mat$color[idx_vec[i]], 
                    lwd = 2, lty = key_mat$lty[idx_vec[i]])
    graphics::points(x = as.numeric(df_param[,"rho"]), y = res_mat[idx_vec2[i],], col = key_mat$color[idx_vec[i]],  
                     pch = key_mat$pch[idx_vec[i]], cex = 1)
  }
  
  legend(legend_loc, key_mat$full_name[idx_vec], 
         lty = key_mat$lty[idx_vec], col = key_mat$color[idx_vec], lwd = 2, 
         pch = key_mat$pch[idx_vec], cex = legend_cex, pt.cex = 0.8)
  
}
