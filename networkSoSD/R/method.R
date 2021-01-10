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
  len <- length(adj_list)
  
  tmp <- array(0, dim = c(nrow(adj_list[[1]]), ncol(adj_list[[1]]), len))
  for(i in 1:len){
    if(verbose && i %% floor(len/10) == 0) cat('*')
    if(method == "sum"){
      tmp[,,i] <- adj_list[[i]]
    } else if(method == "ss"){
      tmp[,,i] <- crossprod(adj_list[[i]])
    } else if(method == "ss_debias"){
      tmp[,,i] <- crossprod(adj_list[[i]]) - diag(colSums(adj_list[[i]]))
    } else {
      stop("method not found")
    }
  }
  
  apply(tmp, c(1,2), sum)
}

#' Spectral clustering
#'
#' @param mat matrix
#' @param K positive integer
#' @param weighted boolean
#' @param row_normalize boolean
#'
#' @return membership vector
#' @export
spectral_clustering <- function(mat, K, weighted = F, row_normalize = F){
  svd_mat <- .svd_projection(mat, K = K, weighted = weighted)
  
  if(row_normalize){
    svd_mat <- t(apply(svd_mat, 1, function(x){x/.l2norm(x)}))
  }
  
  stats::kmeans(svd_mat, centers = K, nstart = 20)$cluster
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

#' Do an SVD projection
#'
#' Uses \code{RSpectra::svds} to compute the \code{k} leading singular vectors, but
#' sometimes there are numerical instability issues. In case of crashes, the code
#' then uses the default \code{svd} function.
#'
#' @param mat numeric matrix with \code{n} rows and \code{n} columns
#' @param K positive integer less than \code{n}
#' @param weighted boolean
#'
#' @return numeric matrix
.svd_projection <- function(mat, K, weighted = F){
  stopifnot(nrow(mat) >= K, ncol(mat) >= nrow(mat))
  
  if(min(dim(mat)) > K+2){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      RSpectra::svds(mat, k = K + 2)
    }, error = function(e){
      svd(mat)
    })
  } else {
    res <- svd(mat)
  }
  
  if(weighted){
    .mult_mat_vec(res$u[,1:K, drop = F], sqrt(abs(res$d[1:K])))
  } else {
    res$u[,1:K,drop = F]
  }
}