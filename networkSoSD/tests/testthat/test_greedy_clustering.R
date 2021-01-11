context("Test greedy clustering")

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

#####################

## .form_array is correct

test_that(".form_array works", {
  set.seed(10)
  adj_list <- generate_dataset()
  
  res <- .form_array(adj_list)
  
  expect_true(is.array(res))
  expect_true(all(dim(res) == c(10, 100, 100)))
})

#####################

## .get_array_centers is correct

test_that(".get_array_centers works", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  idx <- sample(1:2, size = 100, replace = T)
  adj_array <- .form_array(adj_list)
  
  res <- .get_array_centers(adj_array, idx, K)
  
  expect_true(is.array(res))
  expect_true(all(dim(res) == c(10, K, K)))
})

#######################

## .get_array_distances is correct

test_that(".get_array_distances works", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  idx <- sample(1:2, size = 100, replace = T)
  adj_array <- .form_array(adj_list)
  centers <- .get_array_centers(adj_array, idx, K)
  res <- .get_array_distances(adj_array, centers[,idx,])
  
  expect_true(is.array(res))
  expect_true(all(dim(res) == c(100, 2)))
})

##############################

## .assess_distances is correct

test_that(".assess_distances works", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  n <- 100
  idx <- sample(1:2, size = n, replace = T)
  adj_array <- .form_array(adj_list)
  centers <- .get_array_centers(adj_array, idx, K)
  D <- .get_array_distances(adj_array, centers[,idx,])
  
  res <- .assess_distances(D, 1:n, idx)
  
  expect_true(length(res) == n)
  expect_true(is.numeric(res))
})

test_that(".assess_distances is asssessing the correct values", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  n <- 100
  idx <- sample(1:2, size = n, replace = T)
  adj_array <- .form_array(adj_list)
  centers <- .get_array_centers(adj_array, idx, K)
  D <- .get_array_distances(adj_array, centers[,idx,])
  
  res <- .assess_distances(D, 1:n, idx)
  res2 <- numeric(n)
  for(i in 1:n){
    res2[i] <- D[i,idx[i]]
  }
  
  expect_true(sum(abs(res - res2)) <= 1e-4)
})

############################################

## .resolve_empty_clusters is correct

test_that(".resolve_empty_clusters is correct", {
  trials <- 100
  K <- 5
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    idx <- sample(1:K, size = 20, prob = c(1,100,100,100,1), replace = T)
    res <- .resolve_empty_clusters(idx, K)
    
    bool1 <- is.numeric(res) & length(res) == length(idx)
    if(length(unique(idx)) == K){
      bool2 <- all(idx == res); bool3 <- TRUE
    } else {
      non_empty_clusters <- unique(idx)
      empty_clusters <- setdiff(1:K, non_empty_clusters)
      bool2 <- all(sapply(empty_clusters, function(x){length(which(res == x)) == 1}))
      bool3 <- all(sapply(non_empty_clusters, function(x){
        abs(length(which(idx == x)) - length(which(res == x))) <= length(empty_clusters) 
      }))
    }
    
    bool1 & bool2 & bool3
  })
  
  expect_true(all(bool_vec))
})

test_that(".resolve_empty_clusters in extreme instances", {
  trials <- 100
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    res <- .resolve_empty_clusters(rep(3,10), K=3)
    length(res) == 10
  })
  
  expect_true(all(bool_vec))
})

###################################

## .move_indices is correct

test_that(".move_indices is asssessing the correct values", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  n <- 100
  idx <- sample(1:2, size = n, replace = T)
  adj_array <- .form_array(adj_list)
  centers <- .get_array_centers(adj_array, idx, K)
  D <- .get_array_distances(adj_array, centers[,idx,])
  
  res <- .move_indices(D, idx)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("idx", "none_moved"))))
})

test_that(".move_indices correctly moved indices", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  n <- 100
  idx <- sample(1:2, size = n, replace = T)
  adj_array <- .form_array(adj_list)
  centers <- .get_array_centers(adj_array, idx, K)
  D <- .get_array_distances(adj_array, centers[,idx,])
  vec1 <- .assess_distances(D, 1:n, idx)
  
  res <- .move_indices(D, idx)
  vec2 <- .assess_distances(D, 1:n, res$idx)
  
  expect_true(all(vec1 >= vec2))
})

#################################

## greedy_clustering is correct

test_that("greedy_clustering works ", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  
  res <- greedy_clustering(adj_list, K = K)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("cluster", "centers", "obj_vec"))))
})
