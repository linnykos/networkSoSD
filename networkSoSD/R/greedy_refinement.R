# from https://www.dropbox.com/s/o4vbqzd98i41gbp/code.zip?dl=0
# For Lei, Chen, Lynch 2019 Consistent community detection in multi-layer network data.

#' Greedy clustering
#'
#' @param adj_list list of adjacency matrices
#' @param K positive integer
#' @param init_clust initial membership vector. If \code{NA}, we matrixify the tensor and perform K-means
#' @param iter_max positive integer
#' @param verbose boolean
#'
#' @return membership vector
#' @export
greedy_refinement <- function(adj_list, K, init_clust = NA,
                              iter_max = 100, verbose = F){
  stopifnot(length(adj_list) > 1)
  stopifnot(length(unique(as.numeric(sapply(adj_list, dim)))) == 1)
  
  if(all(is.na(init_clust))){
    adj_mat <- flatten(adj_list)
    init_clust <- spectral_clustering(adj_mat, K = K)
  }

  stopifnot(length(init_clust) == nrow(adj_list[[1]]))
  adj_array <- .form_array(adj_list)
  idx <- init_clust; n <- length(init_clust)
  iter <- 0; totsum_DPrev <- Inf; obj_vec <- numeric(0)
  
  # greedily reassign nodes
  while (iter < iter_max){
    if (verbose) cat('step = ', iter, '\n' )
    iter <- iter + 1
    
    # deal with cluster that are empty
    idx <- .resolve_empty_clusters(idx, K)
    
    # compute centers and distances 
    centers <- .get_array_centers(adj_array, idx)
    D <- .get_array_distances(adj_array, centers, idx)
    totsum_D <- sum(.assess_distances(D, idx))
    obj_vec <- c(obj_vec, totsum_D)
    
    # determine which indices need to move, and to where
    idx_prev <- idx 
    idx <- apply(D, 1, which.min)
    
    if(all(idx == idx_prev)) break 
  }
  
  list(cluster = idx, centers = centers, obj_vec = obj_vec)
}

##############################3

.form_array <- function(adj_list){
  stopifnot(length(unique(as.numeric(sapply(adj_list, dim)))) == 1)
  
  len <- length(adj_list)
  tmp <- array(0, dim = c(nrow(adj_list[[1]]), ncol(adj_list[[1]]), len))
  
  for(i in 1:len){
    tmp[,,i] <- adj_list[[i]]
  }
  
  tmp
}

# compute p different kxk matrices that average the SBMs
.get_array_centers <-function(adj_array, idx){
  stopifnot(is.array(adj_array), length(idx) == dim(adj_array)[1],
            all(idx > 0), length(unique(idx)) == max(idx), all(idx %% 1 == 0))
  
  n <- dim(adj_array)[1]
  len <- dim(adj_array)[3]
  K <- max(idx)
  centers <- array(0, dim = c(K, K, len))
  
  for (i in 1:K){
    for (j in i:K){
      centers[i,j,] <- apply(adj_array[idx == i, idx == j, ,drop = F], 3, mean)
      centers[j,i,] <- centers[i,j,]
    }
  }
  
  centers
}

# A: n*n*p
# centers: k*k*p
.get_array_distances <- function(adj_array, centers, idx){
  stopifnot(is.array(adj_array), is.array(centers),
            dim(adj_array)[3] == dim(centers)[3], 
            dim(adj_array)[1] == dim(adj_array)[2],
            dim(centers)[1] == max(idx), dim(centers)[2] == max(idx),
            all(idx > 0), length(unique(idx)) == max(idx), all(idx %% 1 == 0))
  
  n <- dim(adj_array)[1]; K <- max(idx)
  centers_expand <- centers[,idx,,drop = F]
  D <- matrix(0, nrow = n, ncol = K)
  stopifnot(dim(adj_array)[-1] == dim(centers_expand)[-1])
  
  for (j in 1:n){
    for (i in 1:K){
      D[j,i] <- .l2norm(adj_array[j,,] - centers_expand[i,,])^2
    }
  }
  
  D
}

.assess_distances <- function(D, idx){
  stopifnot(length(idx) == nrow(D), max(idx) <= ncol(D))
  
  n <- nrow(D); D[(idx[1:n]-1)*n + 1:n]
}

# assign a randomly-chosen node into empty clusters
.resolve_empty_clusters <- function(idx, K){
  counts <- summary(factor(idx, 1:K))
  empties <- which(counts == 0)
  for (i in empties){
    from <- sample(which(counts > 1), 1)
    if(length(from) == 1){
      lonely <- which(idx == from)[1]
    } else {
      lonely <- sample(which(idx == from), 1)
    }
    idx[lonely] <- i
    counts[i] <- 1
    counts[from] <- counts[from] - 1
  }   
  
  idx
}