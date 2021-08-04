# All functions from https://u.osu.edu/subhadeep/codes/

laplacian <- function(x, weights = rep(1, length(x))) {
  M = length(x)
  laplist <- lapply(1:M, function(m) {
    d = rowSums(x[[m]])
    d = d + mean(d)
    deg <- diag(1 / sqrt(d))
    lap = deg %*% x[[m]] %*% deg
    return(weights[M] * lap)
  })
  return(laplist)
}

tr <- function(mat){
  sum(diag(mat))
}

###### An implementation of Co-regularized Spectral Clustering from Kumar et al. NeurIPS 2011.######

coregfunction <- function(laplist, snrlist, ulist, ustar) {
  M = length(laplist)
  objlist <- lapply(1:M, function(m) {
    specobj <- tr(t(ulist[[m]]) %*% laplist[[m]] %*% ulist[[m]])
    regterm <- (Matrix::norm(t(ulist[[m]]) %*% ustar, type = "F")) ^ 2
    return(specobj + snrlist[m] * regterm)
  })
  obj <- Reduce("+", objlist)
  return(obj)
}


#' Co-regularized  spectral  clustering
#'
#' The main Coreg function. Inputs are the list of adjacency matrices x, number of nodes n, nunber of communities
#' k and regularization parameter beta.
#' @param x list of adjacency matrices
#' @param k number of clusters
#' @param beta positive numeric
#' @param max_iter positive integer
#' @param verbose boolean
#'
#' @return list
#' @export
coreg <- function(x, k, beta, max_iter = 100, verbose = F)
{
  n = nrow(x[[1]])
  M = length(x)
  laplist <- laplacian(x)
  snrlist <- rep(beta, M)
  ulist <- lapply(1:M, function(m) {
    spectra <- RSpectra::eigs_sym(laplist[[m]], k = k)
    specmat <- spectra$vectors[, 1:k]
    return(specmat)
  })
  usumlist <- lapply(1:M, function(m) {
    return(snrlist[m] * ulist[[m]] %*% t(ulist[[m]]))
  })
  usum <- Reduce('+', usumlist)
  ustareigen <- RSpectra::eigs_sym(usum, k = k)
  ustar <- ustareigen$vectors[, 1:k]
  value <- coregfunction(laplist, snrlist, ulist, ustar)
  if(verbose) print(value)
  t = 0
  valdiff = 1
  while ((valdiff > 1e-04) & (t < max_iter))
  {
    valuelast = value
    ulist <- lapply(1:M, function(m) {
      um <- RSpectra::eigs_sym(laplist[[m]] + snrlist[m] * ustar %*% t(ustar), k = k)
      specmat <- um$vectors[, 1:k]
      return(specmat)
    })
    
    usumlist <- lapply(1:M, function(m) {
      return(snrlist[m] * ulist[[m]] %*% t(ulist[[m]]))
    })
    usum <- Reduce('+', usumlist)
    ustareigen <- RSpectra::eigs_sym(usum, k = k)
    ustar <- ustareigen$vectors[, 1:k]
    
    value = coregfunction(laplist, snrlist, ulist, ustar)
    valdiff = (value - valuelast)
    t = t + 1
    if(verbose) print(value)
    if(verbose) print(t)
  }
  specclus <- lapply(ulist, function(r) {
    specu <- stats::kmeans(r, k)
    return(specu$cluster)
  })
  specstar <- stats::kmeans(ustar, k)
  specclus[[M + 1]] <- specstar$cluster
  specclus[[M + 2]] <- Reduce("+", lapply(1:M, function(m) {
    y = clue::solve_LSAP(table(specclus[[m]], specclus[[M + 1]]), maximum = TRUE)
    for (i in 1:n) {
      specclus[[m]][i] <- y[specclus[[m]][i]]
    }
    return(table(specclus[[M + 1]], specclus[[m]]))
  }))
  specclus[[M + 2]] <- specclus[[M + 2]] / rowSums(specclus[[M + 2]])
  return(specclus)
}


################# Implementation of orthogonal LMF used in Paul and Chen, Annals of Statistics 2020.  ##################


## The objective function

lmffunctiono <- function(param, laplist, n, k) {
  M = length(laplist)
  ustar <- matrix(param[1:(n * k)], n, k)
  lambda <-
    lapply(1:M, function(m) {
      return(matrix(param[(n * k + (m - 1) * k ^ 2 + 1):(n * k + m * k ^ 2)], k, k))
    })
  objloop <- sum(unlist(lapply(1:M, function(m) {
    specobj <-
      norm(laplist[[m]] - ustar %*% lambda[[m]] %*% t(ustar), type = "F") ^ 2
    return(specobj)
  })))
  obj = objloop
  return(obj)
}


##  The gradients

lmfdero <- function(param, laplist, n, k) {
  M = length(laplist)
  ustar <- matrix(param[1:(n * k)], n, k)
  lambda <-
    lapply(1:M, function(m) {
      return(matrix(param[(n * k + (m - 1) * k ^ 2 + 1):(n * k + m * k ^ 2)], k, k))
    })
  derlist1 <- lapply(1:M, function(m) {
    specobj = -(diag(n) - ustar %*% t(ustar)) %*% laplist[[m]] %*% ustar %*% lambda[[m]]
    return(specobj)
  })
  derlist2 <- lapply(1:M, function(m) {
    specobj = -t(ustar) %*% (laplist[[m]] - ustar %*% lambda[[m]] %*% t(ustar)) %*%
      ustar
    return(specobj)
  })
  der1 <- Reduce("+", derlist1)
  der2 <- unlist(derlist2)
  return(c(as.vector(der1), as.vector(der2)))
}



#' Linked matrix factorization
#' 
#' The main function with BFGS optimization
#'
#' @param x list of adjacency matrices
#' @param k number of clusters
#'
#' @return list
#' @export
lmfo <- function(x, k) {
  n = nrow(x[[1]])
  M = length(x)
  laplist <- laplacian(x)
  # Initialize with mean laplacian
  lapmean <- Reduce("+", laplist)
  spectra <- RSpectra::eigs_sym(lapmean, k = k)
  ustar <- spectra$vectors[, 1:k]
  lambda <- lapply(1:M, function(m) {
    return(diag(spectra$values[1:k]))
  })
  param <- c(as.vector(ustar), as.vector(unlist(lambda)))
  optimized <-
    stats::optim(
      par = param,
      fn = lmffunctiono,
      gr = lmfdero,
      method = "BFGS",
      control = list(reltol = 0.0001, maxit = 200),
      laplist = laplist,
      n = n,
      k = k
    )
  param <- optimized$par
  
  ustar <- matrix(param[1:(n * k)], n, k)
  lambda <-
    lapply(1:M, function(m) {
      return(matrix(param[(n * k + (m - 1) * k ^ 2 + 1):(n * k + m * k ^ 2)], k, k))
    })
  
  specstar <- stats::kmeans(ustar, k)
  specclus <- specstar$cluster
  return(specclus)
  
}
