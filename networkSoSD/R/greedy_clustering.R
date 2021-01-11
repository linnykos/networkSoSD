# from https://www.dropbox.com/s/o4vbqzd98i41gbp/code.zip?dl=0
# For Lei, Chen, Lynch 2019 Consistent community detection in multi-layer network data.

#' Greedy clustering
#'
#' @param adj_list list of adjacency matrices
#' @param K positive integer
#' @param nstart positive integer
#' @param iter_max positive integer
#' @param verbose boolean
#'
#' @return membership vector
#' @export
greedy_clustering <- function(adj_list, K, nstart = 10, iter_max = 100, verbose = F){
  stopifnot(length(adj_list) > 1)
  stopifnot(length(unique(as.numeric(sapply(adj_list, dim)))) == 1)
  adj_array <- .form_array(adj_list)
  
  totsum_DBest <- Inf; idx_Best <- NA; center_Best <- NA; objvec_Best <- NA
  dim_vec <- dim(adj_array); n <- dim_vec[3]
  adj_vec <- t(matrix(adj_array, dim_vec[1]*dim_vec[2], dim_vec[3]))
  
  # if there are too few unique rows
  if (dim(unique(adj_vec))[1] < K) {
    return(list(idx = idx_Best, centers = center_Best))
  }
  
  # do nstart number of restarts
  for (rep in 1:nstart){
    idx <- stats::kmeans(adj_vec, K)$cluster 
    iter <- 0; totsum_DPrev <- Inf; objvec <- numeric(0)
    
    # greedily reassign nodes
    while (iter < iter_max){
      if (verbose) cat('step = ', iter, '\n' )
      iter <- iter + 1
      
      # deal with cluster that are empty
      idx <- .resolve_empty_clusters(idx, K)
      
      # compute centers and distances 
      centers <- .get_array_centers(adj_array, idx, K)
      D <- .get_array_distances(adj_array, centers[,idx,,drop = F])
      totsum_D <- sum(.assess_distances(D, 1:n, idx))
      if (totsum_DPrev <  totsum_D | totsum_DPrev == totsum_D) break  
      totsum_DPrev <- totsum_D; objvec <- c(objvec, totsum_D)
      
      # determine which indices need to move, and to where
      res <- .move_indices(D, idx)
      if(res$none_moved) break else idx <- res$idx
    }
    
    if (totsum_D < totsum_DBest){
      totsum_DBest <- totsum_D; idx_Best <- idx; center_Best <- centers; objvec_Best <- objvec
    }
  }
  
  list(cluster = idx_Best, centers = center_Best, obj_vec = objvec_Best)
}



##############################3

.form_array <- function(adj_list){
  stopifnot(length(unique(as.numeric(sapply(adj_list, dim)))) == 1)
  
  len <- length(adj_list)
  tmp <- array(0, dim = c(len, nrow(adj_list[[1]]), ncol(adj_list[[1]])))
  
  for(i in 1:len){
    tmp[i,,] <- adj_list[[i]]
  }
  
  tmp
}

# compute p different kxk matrices that average the SBMs
.get_array_centers <-function(adj_array, idx, K){
  stopifnot(is.array(adj_array))
  
  p <- dim(adj_array)[1]
  n <- dim(adj_array)[2]
  idx_n <- idx[1:n]
  centers <- array(0, dim = c(p,K,K))
  
  for (i in 1:K){
    for (j in i:K){
      # only one element in the cluster
      if ((sum(idx_n == i) == 1) & (sum(idx == j ) == 1) ) {
        centers[,i,j] <- 0.5
        centers[,j,i] <- centers[,i,j]
        
      } else{
        centers[,i,j] <- apply(adj_array[, idx_n == i, idx == j], 1, mean)
        centers[,j,i] <- centers[,i,j]
      }
    }
  }
  
  centers
}

# A: p*n*n
# centers: p*n*k
.get_array_distances <- function(adj_array, centers){
  stopifnot(is.array(adj_array), is.array(centers),
            dim(adj_array)[1] == dim(centers)[1], 
            dim(adj_array)[2] == dim(centers)[2],
            dim(adj_array)[2] == dim(adj_array)[3])
  
  k <- dim(centers)[3]
  n <- dim(adj_array)[3]
  D <- matrix(0, nrow = n, ncol = k)
  
  for (j in 1:n){
    for (i in 1:k){
      D[j,i] <- .l2norm(adj_array[,,j] - centers[,,i])^2
    }
  }
  
  D
}

.assess_distances <- function(D, n_idx, cluster_assignment){
  stopifnot(length(n_idx) <= nrow(D), length(cluster_assignment) == nrow(D),
            max(cluster_assignment) <= ncol(D))
  
  n <- nrow(D); D[(cluster_assignment[n_idx]-1)*n + n_idx]
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

.move_indices <- function(D, idx){
  n <- nrow(D)
  idx_min <- apply(D, 1, which.min)
  moved <- which(idx_min != idx)
  
  # resolve tie in favor of not move
  if (length(moved) > 0) {
    moved <- moved[.assess_distances(D, moved, idx) > apply(D[moved,,drop = F], 1, min)] 
    idx[moved] <- idx_min[moved]
  } 
  
  list(idx = idx, none_moved = (length(moved) == 0))
}

