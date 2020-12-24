context("Test aggregate")

generate_dataset <- function(){
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B1 <- round(eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat), 4)
  B2 <- round(eigen_mat %*% diag(c(1.5, 0.2, -0.4)) %*% t(eigen_mat), 4)
  rho <- 0.2
  
  n <- 100
  L <- 10
  membership_vec <- c(rep(1, .6*n), rep(2, .3*n), rep(3, .1*n))
  
  prob_mat1 <- compute_prob_mat(rho*B1, membership_vec)
  prob_mat2 <- compute_prob_mat(rho*B2, membership_vec)
  prob_list <- lapply(1:L, function(i){if(i <= L/2) prob_mat1 else prob_mat2})
  
  lapply(1:L, function(i){generate_adjaceny_mat(prob_list[[i]])})
}

## aggregate_networks is correct

test_that("aggregate_networks works for ss_debias", {
  set.seed(10)
  adj_list <- generate_dataset()

  res <- aggregate_networks(adj_list, method = "ss_debias")
  
  expect_true(nrow(res) == ncol(res))
  expect_true(sum(abs(res - t(res))) <= 1e-4)
})

test_that("aggregate_networks works for ss", {
  set.seed(10)
  adj_list <- generate_dataset()
  
  res <- aggregate_networks(adj_list, method = "ss")
  
  expect_true(nrow(res) == ncol(res))
  expect_true(sum(abs(res - t(res))) <= 1e-4)
})

test_that("aggregate_networks works for sum", {
  set.seed(10)
  adj_list <- generate_dataset()
  
  res <- aggregate_networks(adj_list, method = "sum")
  
  expect_true(nrow(res) == ncol(res))
  expect_true(sum(abs(res - t(res))) <= 1e-4)
})

###################################

## spectral_clustering is correct

test_that("spectral_clustering works", {
  set.seed(10)
  adj_list <- generate_dataset()
  
  agg_mat <- aggregate_networks(adj_list, method = "ss_debias")
  res <- spectral_clustering(agg_mat, K = 3, weighted = F)
  
  expect_true(length(res) == nrow(agg_mat))
  expect_true(all(res > 0))
  expect_true(all(res %% 1 == 0))
  expect_true(all(1:max(res) %in% res))
})