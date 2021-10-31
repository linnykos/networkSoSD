plot_table <- function(cluster_vec1, cluster_vec2,
                       main, 
                       row_offset = 0,
                       col_offset = 0,
                       col_xshift = 0,
                       normalize_row = T,
                       row_names = paste0("C", sort(unique(cluster_vec1))),
                       col_names = paste0("C", sort(unique(cluster_vec2))),
                       asp = length(unique(cluster_vec1))/length(unique(cluster_vec2)),
                       base_color = grDevices::rgb(0.803, 0.156, 0.211),
                       sig_digs = 2,
                       cex_text = 1,
                       cex_label = 0.9,
                       ...){
  stopifnot(length(cluster_vec1) == length(cluster_vec2))
  
  tab <- table(cluster_vec1, cluster_vec2)
  p1 <- nrow(tab); p2 <- ncol(tab)
  if(normalize_row){
    for(i in 1:p1) tab[i,] <- tab[i,]/sum(tab[i,])
  } else {
    for(j in 1:p2) tab[,j] <- tab[,j]/sum(tab[,j])
  }
  col_ramp <- grDevices::colorRampPalette(c("white", base_color))(100)
  
  y_tic <- seq(-.5, p1-.5, by = 1)/(p1-1)
  x_tic <- seq(-.5, p2-.5, by = 1)/(p2-1)
  graphics::image(.rotate(tab), 
                  col = col_ramp, 
                  asp = asp, 
                  axes = F,
                  xlim = c(min(x_tic)-.05, max(x_tic)+.05),
                  ylim = c(min(y_tic)-.05, max(y_tic)+.05), 
                  main = main, ...)
  for(y in y_tic){
    graphics::lines(range(x_tic), rep(y, 2))
  }
  for(x in x_tic){
    graphics::lines(rep(x, 2), range(y_tic))
  }
  graphics::text(par("usr")[3] + row_offset, 
                 seq(0,1,length.out = p1), 
                 adj = 1, 
                 labels = rev(row_names), 
                 xpd = TRUE,
                 cex = cex_label)
  graphics::text(seq(0, 1, length.out = p2) + col_xshift, 
                 par("usr")[1] + col_offset, 
                 labels = col_names, 
                 srt = -45,
                 xpd = TRUE,
                 cex = cex_label)
  
  x_vec <- seq(0, 1, length.out = p2)
  y_vec <- seq(0, 1, length.out = p1)
  for(i in 1:length(x_vec)){
    for(j in 1:length(y_vec)){
      graphics::text(x_vec[i], 1-y_vec[j], 
                     label = round(tab[j,i], sig_digs),
                     cex = cex_text)
    }
  }
  
  invisible()
}

.rotate <- function(mat){t(mat)[,nrow(mat):1]}
