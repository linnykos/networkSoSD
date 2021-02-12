context("Test greedy clustering")

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

#####################

## .form_array is correct

test_that(".form_array works", {
  set.seed(10)
  adj_list <- generate_dataset()
  
  res <- .form_array(adj_list)
  
  expect_true(is.array(res))
  expect_true(all(dim(res) == c(100, 100, 10)))
})

test_that(".form_array gives the same values as aggregate_networks when summed", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    adj_list <- generate_dataset()
    
    val <- apply(.form_array(adj_list), c(1,2), sum)
    val2 <- aggregate_networks(adj_list, method = "sum", verbose = F)
    
    sum(abs(val - val2)) <= 1e-6
  })
  
  expect_true(all(bool_vec))
})

#####################

## .get_array_centers is correct

test_that(".get_array_centers works", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  idx <- sample(1:2, size = 100, replace = T)
  adj_array <- .form_array(adj_list)
  
  res <- .get_array_centers(adj_array, idx)
  
  expect_true(is.array(res))
  expect_true(all(dim(res) == c(K, K, 10)))
})

test_that(".get_array_centers is properly extracting", {
  vec1 <- c(1:10)
  vec2 <- (c(1:10)+1)*2
  vec3 <- (c(1:10)+2)*3
  len <- length(vec1)
  adj_list <- vector("list", length = len)
  
  for(i in 1:len){
    mat <- matrix(0, 6, 6)
    mat[1:5,1:5] <- vec1[i]
    mat[1:5,6] <- vec2[i]; mat[6,1:5] <- vec2[i]
    mat[6,6] <- vec3[i]
    adj_list[[i]] <- mat
  }
  
  adj_array <- .form_array(adj_list); idx <- c(rep(1,5), 2)
  res <- .get_array_centers(adj_array, idx)
  
  expect_true(sum(abs(res[1,1,] - vec1)) <= 1e-6)
  expect_true(sum(abs(res[1,2,] - vec2)) <= 1e-6)
  expect_true(sum(abs(res[2,1,] - vec2)) <= 1e-6)
  expect_true(sum(abs(res[2,2,] - vec3)) <= 1e-6)
})

#######################

## .get_array_distances is correct

test_that(".get_array_distances works", {
  set.seed(10)
  adj_list <- generate_dataset()
  n <- 100
  idx <- sample(1:2, size = n, replace = T)
  adj_array <- .form_array(adj_list)
  centers <- .get_array_centers(adj_array, idx)
  res <- .get_array_distances(adj_array, centers, idx)
  
  expect_true(is.array(res))
  expect_true(all(dim(res) == c(100, 2)))
})

test_that(".get_array_distances computes meaningful distances", {
  set.seed(10)
  adj_list <- generate_dataset(rho = 1)
  n <- 100
  mem_prop1 <- 0.4; mem_prop2 <- 0.1; mem_prop3 <- 0.5
  membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
  adj_array <- .form_array(adj_list)
  centers <- .get_array_centers(adj_array, membership_vec)
  res <- .get_array_distances(adj_array, centers, membership_vec)

  idx_min <- apply(res, 1, which.min)
  
  expect_true(all(idx_min == membership_vec))
})

test_that(".get_array_distances computes the correct values", {
  trials <- 50
  
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    adj_list <- generate_dataset(rho = 1)
    K <- 3; n <- 100
    mem_prop1 <- 0.4; mem_prop2 <- 0.1; mem_prop3 <- 0.5
    membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
    adj_array <- .form_array(adj_list)
    L <- length(adj_list)
    centers <- .get_array_centers(adj_array, membership_vec)
    res <- .get_array_distances(adj_array, centers, membership_vec)
    
    res2 <- matrix(0, nrow = n, ncol = K)
    for(i in 1:n){
      for(j in 1:K){
        for(l in 1:L){
          res2[i,j] <- res2[i,j] + .l2norm(adj_array[i,,l] - centers[j,membership_vec,l])^2
        }
      }
    }
    
    sum(abs(res - res2)) <= 1e-6
  })

  
  expect_true(all(bool_vec))
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
  centers <- .get_array_centers(adj_array, idx)
  D <- .get_array_distances(adj_array, centers, idx)
  
  res <- .assess_distances(D, idx)
  
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
  centers <- .get_array_centers(adj_array, idx)
  D <- .get_array_distances(adj_array, centers, idx)
  
  res <- .assess_distances(D, idx)
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

#################################

## greedy_refinement is correct

test_that("greedy_refinement works ", {
  set.seed(10)
  adj_list <- generate_dataset()
  K <- 2
  
  res <- greedy_refinement(adj_list, K = K)
  
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("cluster", "centers", "obj_vec"))))
})

# test_that("greedy_refinement improves upon tensoring ", {
#   trials <- 50; K <- 3; n <- 100
#   mem_prop1 <- 0.4; mem_prop2 <- 0.1; mem_prop3 <- 0.5
#   membership_vec <- c(rep(1, mem_prop1*n), rep(2, mem_prop2*n), rep(3, mem_prop3*n))
#   
#   bool_vec <- sapply(1:trials, function(x){
#     set.seed(x)
#     adj_list <- generate_dataset()
#     
#     flat_mat <- networkSoSD::flatten(adj_list)
#     res1 <- networkSoSD::spectral_clustering(flat_mat, K = K, weighted = F)
#     
#     res2 <- greedy_refinement(adj_list, init_clust = res1, K = K)
#     
#     vec <- align_two_membership_vectors(membership_vec, res1)
#     tab1 <- table(membership_vec, vec)
#     obj1 <- sum(diag(tab1))/sum(tab1)
#     
#     vec <- align_two_membership_vectors(membership_vec, res2$cluster)
#     tab2 <- table(membership_vec, vec)
#     obj2 <- sum(diag(tab2))/sum(tab2)
#     
#     obj1 <= obj2
#   })
#   
#   expect_true(all(bool_vec))
# })
