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
#' @param sym_mat symmetric matrix
#' @param K positive integer
#' @param weighted boolean
#'
#' @return membership vector
#' @export
spectral_clustering <- function(sym_mat, K, weighted = F){
  svd_mat <- .svd_projection(sym_mat, K = K, weighted = weighted)
  
  stats::kmeans(svd_mat, centers = K, nstart = 20)$cluster
}

##################

#' Do an SVD projection
#'
#' Uses \code{RSpectra::svds} to compute the \code{k} leading singular vectors, but
#' sometimes there are numerical instability issues. In case of crashes, the code
#' then uses the default \code{svd} function.
#'
#' @param adj_mat numeric matrix with \code{n} rows and \code{n} columns
#' @param K positive integer less than \code{n}
#' @param weighted boolean
#'
#' @return numeric matrix
.svd_projection <- function(adj_mat, K, weighted = F){
  stopifnot(nrow(adj_mat) >= K, nrow(adj_mat) == ncol(adj_mat))
  
  if(nrow(adj_mat) > K+2){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      RSpectra::svds(adj_mat, k = K + 2)
    }, error = function(e){
      svd(adj_mat)
    })
  } else {
    res <- svd(adj_mat)
  }
  
  if(weighted){
    diag_mat <- .diag_matrix(sqrt(res$d[1:K]))
    res$u[,1:K, drop = F] %*% diag_mat
  } else {
    res$u[,1:K,drop = F]
  }
}

.diag_matrix <- function(vec){
  K <- length(vec)
  if(K == 1) {
    matrix(vec, 1, 1)
  } else {
    diag(vec)
  }
}
