#' Aggregate networks
#'
#' @param adj_list list of adjacency matrices
#' @param method method, either \code{"ss_debias"}, \code{"ss"} or \code{"sum"}
#' @param verbose boolean
#'
#' @return a matrix
#' @export
aggregate_networks <- function(adj_list, method = "ss_debias", 
                               verbose = F){
  stopifnot(length(unique(as.numeric(sapply(adj_list, dim)))) == 1)
  n <- nrow(adj_list[[1]]); len <- length(adj_list)
  
  tmp <- array(0, dim = c(n, n, len))
  for(i in 1:len){
    if(verbose && i %% floor(len/10) == 0) cat('*')
    if(method == "sum"){
      tmp[,,i] <- adj_list[[i]]
    } else if(method == "ss"){
      tmp[,,i] <- crossprod(adj_list[[i]])
    } else if(method == "ss_debias"){
      tmp[,,i] <- crossprod(adj_list[[i]]); diag(tmp[,,i]) <- 0
      # equivalent to: tmp[,,i] <- crossprod(adj_list[[i]]) - diag(colSums(adj_list[[i]])) 
    } else if(method == "ss_debias2"){
      tmp[,,i] <- crossprod(adj_list[[i]]); diag(tmp[,,i]) <- 0; diag(tmp[,,i]) <- colSums(tmp[,,i])/n
    } else {
      stop("method not found")
    }
  }
  
  apply(tmp, c(1,2), sum)
}

#' Spectral clustering
#' 
#' When the model is DC-SBM, both \code{weighted} and \code{row_normalize} should be \code{TRUE}.
#'
#' @param mat matrix
#' @param K positive integer
#' @param weighted boolean
#' @param row_normalize boolean.
#'
#' @return membership vector
#' @export
spectral_clustering <- function(mat, K, weighted = F, row_normalize = F){
  svd_res <- .svd_truncated(mat, K = K, symmetric = T)
  if(weighted){
    svd_mat <- .mult_mat_vec(svd_res$u, svd_res$d)
  } else {
    svd_mat <- svd_res$u * sqrt(nrow(mat)) # multiplying for numerical stability
  }
  
  if(row_normalize){
    svd_mat <- t(apply(svd_mat, 1, function(x){x/.l2norm(x)}))
  }
  
  .safe_kmeans(svd_mat, K)
}

#' Flatten collection of adjacency matrices
#'
#' @param adj_list list of adjacency matrices
#'
#' @return a matrix
#' @export
flatten <- function(adj_list){
  stopifnot(length(unique(as.numeric(sapply(adj_list, dim)))) == 1)
  n <- nrow(adj_list[[1]]); L <- length(adj_list)
  
  mat <- matrix(NA, nrow = n, ncol = n*L)
  for(l in 1:L){
    mat[,((l-1)*n+1):(l*n)] <- adj_list[[l]]
  }
  
  mat
}

##################

.safe_kmeans <- function(mat, K){
  tryCatch({
    stats::kmeans(mat, centers = K, nstart = 10)$cluster
  }, error = function(e){
    range_vec <- apply(mat, 2, function(x){abs(diff(range(x)))})
    stopifnot(!all(range_vec == 0))
    idx <- which.min(range_vec)
    
    mat[,idx] <- mat[,idx] + stats::rnorm(nrow(mat), sd = min(range_vec[range_vec > 0]))
    
    stats::kmeans(mat, centers = K, nstart = 10)$cluster
  })
}
