rm(list = ls()); dcsbm = T; K = 3
vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- round(eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat), 4)
B2 <- round(eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat), 4)
rho <- 1

n <- 200 # ifelse(dcsbm, 500, 100)
L <- 10
if(dcsbm){
  membership_vec <- c(rep(1, .4*n), rep(2, .3*n), rep(3, .3*n))
} else {
  membership_vec <- c(rep(1, .6*n), rep(2, .3*n), rep(3, .1*n))
}

prob_mat1 <- compute_prob_mat(rho*B1, membership_vec)
prob_mat2 <- compute_prob_mat(rho*B2, membership_vec)

if(dcsbm) {
  degree_vec <- c(seq(0.1, 1, length.out = 0.4*n), 
                  seq(0.1, 1, length.out = 0.3*n), 
                  seq(0.1, 1, length.out = 0.3*n))
  
  prob_mat1 <- .mult_mat_vec(.mult_vec_mat(degree_vec, prob_mat1), degree_vec)
  prob_mat2 <- .mult_mat_vec(.mult_vec_mat(degree_vec, prob_mat2), degree_vec)
}

prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})

####

# zz <- colSums(prob_list[[1]]); plot(zz)

agg_mat <- aggregate_networks(prob_list, method = "ss")
res <- spectral_clustering(agg_mat, K = 3, weighted = F, row_normalize = F)
res

svd_mat <- .svd_projection(agg_mat, K = K, weighted = F)
image(t(svd_mat))

#######

set.seed(10)
adj_list <- lapply(1:L, function(i){generate_adjaceny_mat(prob_list[[i]])})

n <- nrow(adj_list[[1]]); L <- length(adj_list)
agg_mat <- aggregate_networks(adj_list, method = "ss_debias")
res <- spectral_clustering(agg_mat, K = 3, weighted = F, row_normalize = T)
res

res <- spectral_clustering(agg_mat, K = 3, weighted = F, row_normalize = F)
res
