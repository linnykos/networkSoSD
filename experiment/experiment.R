rm(list=ls())
trials <- 50; K <- 3; n <- 100
mem_prop1 <- 0.4; mem_prop2 <- 0.1; mem_prop3 <- 0.5
membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
x = 4

generate_dataset <- function(rho = 0.5){
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B1 <- eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat)
  B2 <- eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat)
  K <- 3
  
  n <- 100; L <- 10
  mem_prop1 <- 0.4; mem_prop2 <- 0.1; mem_prop3 <- 0.5
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  if(length(membership_vec) < n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))
  
  prob_mat1 <- compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  lapply(1:L, function(i){generate_adjaceny_mat(prob_list[[i]])})
}

set.seed(x)
adj_list <- generate_dataset()

flat_mat <- networkSoSD::flatten(adj_list)
res1 <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = F)

# res2 <- greedy_refinement(adj_list, init_clust = res1, K = K)
# 
# vec <- align_two_membership_vectors(membership_vec, res1)
# tab1 <- table(membership_vec, vec)
# obj1 <- sum(diag(tab1))/sum(tab1)
# 
# vec <- align_two_membership_vectors(membership_vec, res2$cluster)
# tab2 <- table(membership_vec, vec)
# obj2 <- sum(diag(tab2))/sum(tab2)

##################################################33
init_clust <- res1
adj_array <- .form_array(adj_list)
idx <- init_clust; n <- length(init_clust)
# deal with cluster that are empty
idx <- .resolve_empty_clusters(idx, K)
# compute centers and distances 
centers <- .get_array_centers(adj_array, idx)
D <- .get_array_distances(adj_array, centers, idx)
image(t(D))
