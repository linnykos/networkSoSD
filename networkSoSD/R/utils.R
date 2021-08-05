# requires both vec1 and vec2 to have the same number of unique clusters
align_two_membership_vectors <- function(vec1, vec2, override = F){
  if(override & length(unique(vec1)) > length(unique(vec2))){
    len <- length(unique(vec1))
    vec2[1:len] <- 1:len
  }
  stopifnot(length(unique(vec1)) == length(unique(vec2)))
  
  K <- length(unique(vec1))
  tab <- table(vec1, vec2)
  
  permn_list <- .permn(K)
  similarity_vec <- sapply(permn_list, function(x){
    tab2 <- tab
    tab2 <- tab2[,x]
    sum(diag(tab2))/sum(tab2)
  })
  
  idx <- which.max(similarity_vec)
  
  vec2_new <- vec2
  for(i in 1:K){
    vec2_new[which(vec2 == permn_list[[idx]][i])] <- i
  }
  
  vec2_new
}


## https://github.com/cran/combinat/blob/master/R/permn.R
.permn <- function(x, fun = NULL, ...) {
  if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x) x <- seq(x)
  n <- length(x)
  nofun <- is.null(fun)
  out <- vector("list", gamma(n + 1))
  p <- ip <- seqn <- 1:n
  d <- rep(-1, n)
  d[1] <- 0
  m <- n + 1
  p <- c(m, p, m)
  i <- 1
  use <-  - c(1, n + 2)
  
  while(m != 1) {
    out[[i]] <- if(nofun) x[p[use]] else fun(x[p[use]], ...)
    i <- i + 1
    m <- n
    chk <- (p[ip + d + 1] > seqn)
    m <- max(seqn[!chk])
    if(m < n)
      d[(m + 1):n] <-  - d[(m + 1):n]
    index1 <- ip[m] + 1
    index2 <- p[index1] <- p[index1 + d[m]]
    p[index1 + d[m]] <- m
    tmp <- ip[index2]
    ip[index2] <- ip[m]
    ip[m] <- tmp
  }
  
  out
}

######################

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == nrow(mat))
  vec * mat
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(is.matrix(mat), !is.matrix(vec), length(vec) == ncol(mat))
  mat * rep(vec, rep(nrow(mat), length(vec)))
}


.svd_truncated <- function(mat, K, symmetric){
  if(is.na(K)) K <- min(dim(mat))
  stopifnot(min(dim(mat)) >= K)
  
  if(min(dim(mat)) > 2*(K+2)){
    res <- tryCatch({
      # ask for more singular values than needed to ensure stability
      if(symmetric){
        tmp <- irlba::partial_eigen(mat, n = K)
        list(u = tmp$vectors, d = tmp$values, v = tmp$vectors)
      } else {
        irlba::irlba(mat, nv = K)
      }
    }, warning = function(e){
      RSpectra::svds(mat, k = K)
    }, error = function(e){
      RSpectra::svds(mat, k = K)
    })
  } else {
    res <- svd(mat)
  }
  
  res$u <- res$u[,1:K, drop = F]; res$v <- res$v[,1:K, drop = F]; res$d <- res$d[1:K]
  
  res
}
