context("Test generate network")

## .generate_membership_matrix is correct

test_that(".generate_membership_matrix works", {
  res <- .generate_membership_matrix(sample(1:5, size = 100, replace = T))
  expect_true(all(dim(res) == c(100, 5)))
  expect_true(all(apply(res, 1, function(x){sum(x) == 1 & length(which(x==1)) == 1})))
})

## compute_prob_mat is correct

test_that("compute_prob_mat works", {
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B <- round(eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat), 4)
  rho <- 0.5
  
  n <- 500
  membership_vec <- c(rep(1, .6*n), rep(2, .3*n), rep(3, .1*n))
  
  prob_mat <- compute_prob_mat(rho*B, membership_vec)
  
  expect_true(nrow(prob_mat) == ncol(prob_mat))
  expect_true(sum(abs(prob_mat - t(prob_mat))) <= 1e-4)
})

## generate_adjaceny_mat is correct

test_that("generate_adjaceny_mat works", {
  vec1 <- c(1,1,sqrt(2))
  vec2 <- c(1,1,-sqrt(2))
  vec3 <- c(-1,1,0)
  eigen_mat <- cbind(vec1/.l2norm(vec1), vec2/.l2norm(vec2), vec3/.l2norm(vec3))
  B <- round(eigen_mat %*% diag(c(1.5, 0.2, 0.4)) %*% t(eigen_mat), 4)
  rho <- 0.5
  
  n <- 500
  membership_vec <- c(rep(1, .6*n), rep(2, .3*n), rep(3, .1*n))
  
  prob_mat <- compute_prob_mat(rho*B, membership_vec)
  adj_mat <- generate_adjaceny_mat(prob_mat)
  
  expect_true(nrow(adj_mat) == ncol(adj_mat))
  expect_true(sum(abs(adj_mat - t(adj_mat))) <= 1e-4)
  expect_true(all(adj_mat %in% c(0,1)))
})
