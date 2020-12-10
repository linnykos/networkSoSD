rm(list=ls())
# vec1 <- c(1,2,2)/3
# vec2 <- c(-2,-1,2)/3
# vec3 <- c(2,-2,1)/3
# eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
# B1 <- eigen_mat %*% diag(c(0.8, 0.5, 0.3)) %*% t(eigen_mat)
# B2 <- eigen_mat %*% diag(c(0.9, 0.3, -0.3)) %*% t(eigen_mat)
# rho <- 0.1
# n <- 500
# L <- 100
# membership_vec <- c(rep(1, .2*n), rep(2, .2*n), rep(3, .6*n))

vec1 <- c(1,1,sqrt(2))
vec2 <- c(1,1,-sqrt(2))
vec3 <- c(-1,1,0)
eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
B1 <- eigen_mat %*% diag(c(0.9, 0.2, 0.4)) %*% t(eigen_mat)
B2 <-eigen_mat %*% diag(c(0.9, 0.2, -0.4)) %*% t(eigen_mat)
rho <- 0.1
n <- 500
L <- 100
membership_vec <- c(rep(1, .45*n), rep(2, .1*n), rep(3, 0.45*n))

if(length(membership_vec) != n) membership_vec <- c(membership_vec, rep(3, n-length(membership_vec)))

#################

prob_mat1 <- compute_prob_mat(rho*B1, membership_vec)
prob_mat2 <- compute_prob_mat(rho*B2, membership_vec)
prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})

set.seed(10)
adj_list <- lapply(1:L, function(i){generate_adjaceny_mat(prob_list[[i]])})

###############################

# # see if we can make ensure it'll fail
# tmp <- array(0, dim = c(nrow(prob_list[[1]]), ncol(prob_list[[1]]), L))
# for(i in 1:L){
#   tmp[,,i] <- prob_list[[i]]
# }
# agg_prob <- apply(tmp, c(1,2), sum)
# 
# tmp <- adj_list[[1]]
# eig_vec <- eigen(tmp)$vector[,1:3]
# 
# deg_vec <- rowSums(adj_list[[1]])
# deg_vec <- deg_vec/.l2norm(deg_vec)
# 
# proj_vec <- eig_vec %*% t(eig_vec) %*% deg_vec
# # resid_vec <- deg_vec - proj_vec
# # as.numeric(t(resid_vec) %*% proj_vec)
# # sqrt(sum(proj_vec^2))
# # sqrt(sum(resid_vec^2))
# acos(as.numeric(t(proj_vec) %*% deg_vec)/(.l2norm(proj_vec) * .l2norm(deg_vec))) * 180/pi

################################

# try the desired method
agg_network <- aggregate_networks(adj_list, method = "ss_debias")
res <- spectral_clustering(agg_network, K = 3, weighted = F)
table(res, membership_vec)

# try naive method of adding
agg_network <- aggregate_networks(adj_list, method = "sum")
res <- spectral_clustering(agg_network, K = 3, weighted = F)
table(res, membership_vec)

# try naive method of ss
agg_network <- aggregate_networks(adj_list, method = "ss")
res <- spectral_clustering(agg_network, K = 3, weighted = F)
table(res, membership_vec)

### now all the weighted versions
# try naive method of adding
agg_network <- aggregate_networks(adj_list, method = "sum")
res <- spectral_clustering(agg_network, K = 3, weighted = T)
table(res, membership_vec)

# try naive method of ss
agg_network <- aggregate_networks(adj_list, method = "ss")
res <- spectral_clustering(agg_network, K = 3, weighted = T)
table(res, membership_vec)

