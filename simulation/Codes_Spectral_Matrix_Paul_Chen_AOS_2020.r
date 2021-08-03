library(kernlab)
library(clue)
library(psych)
library(igraph)


######## Codes to implement methods and reproduce the simulations in Paul and Chen, Annals of Statistics, 2020.####
### The paper:
#### Paul, S., & Chen, Y. (2020). Spectral and matrix factorization methods for consistent community detection in multi-layer networks. The Annals of Statistics, 48(1), 230-250.


##### Note many of the implementations here use the list of Laplacian matrices instead of Adjacency matrices
##### If desired, replace the Laplacian matrices with Adjacency matrices by changing "laplist" to "x" .

mycheck2 <- function(x)
{
  if (x < 10 ^ (-6) || is.infinite(x) || is.nan(x))
  {
    x = 10 ^ (-6)
  }
  return(x)
}


## Calculate NMI between two community assignments
nmi <- function(clus, comm)
{
  n = length(clus)
  entclus = sum((table(clus) / n) * log(table(clus) / n))
  entcomm = sum((table(comm) / n) * log(table(comm) / n))
  mi = sum((table(clus, comm) / n) * log(n * table(clus, comm) / (table(clus) %*%
                                                                    t(table(
                                                                      comm
                                                                    )))), na.rm = TRUE)
  nmi = -mi / ((entclus + entcomm) / 2)
  return(nmi)
}



### MCrate calculates misclassification rate between two community assignments
pMatrix.min <- function(A, B) {
  n <- nrow(A)
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      D[j, i] <- (sum((B[j,] - A[i,]) ^ 2))
    }
  }
  vec <- c(clue::solve_LSAP(D))
  list(A = A[vec, ], pvec = vec)
}

# estimate, true
MCrate <- function (E, C) {
  E = factor(E)
  C = factor(C)
  N = length(E)
  k = length(levels(C))
  E = factor(x = E, levels = 1:k)
  A = table(E, C)
  X <- pMatrix.min(A, diag(1, k))
  A = X$A
  N_err = sum(A) - sum(diag(A))
  P_err = N_err / N
  return(1 - P_err)
}




### A generic Spectral Clustering of Laplacian matrix with row normalization and regularization
#### This function will be used to initialize some optimization methods
dcspectral <- function(x, n, k)
{
  d = rowSums(x)
  d = d + mean(d)
  deg <- diag(1 / sqrt(d))
  l = deg %*% x %*% deg
  spectra <- RSpectra::eigs_sym(l, k = k)
  specmat <- spectra$vectors[, 1:k]
  rownorm <- apply(specmat, 1, function(a) {
    (sum(a ^ 2)) ^ 0.5
  })
  rownorm <- ifelse(rownorm < 10 ^ (-06), 10 ^ (-06), rownorm)
  specnorm <- specmat / rownorm
  speck <- kmeans(specnorm, k)
  return(speck$cluster)
}






########## From here onwards the input x is a list of adjacency matrices, n is the number of nodes and k is the
#########  number of communities.

#### Create the list of Laplacian matrices

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


# The main Coreg function. Inputs are the list of adjacency matrices x, number of nodes n, nunber of communities
# k and regularization parameter beta.
coreg <- function(x, n, k, beta, max_iter = 100, verbose = F)
{
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
  print(value)
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




### Mean Laplacian matrix community detection
meancluster <- function(x, n, k)
{
  laplist <- laplacian(x)
  lapmean <- Reduce("+", laplist)
  spectra <- RSpectra::eigs_sym(lapmean, k = k)
  umean <- spectra$vectors[, 1:k]
  specmean <- stats::kmeans(umean, k)
  return(specmean$cluster)
}


### Laplacian Aggregate Spectral Kernel
speck <- function(x, n, k) {
  M = length(x)
  laplist <- laplacian(x)
  ulist <- lapply(1:M, function(m) {
    spectra <- RSpectra::eigs_sym(laplist[[m]], k = k)
    specmat <- spectra$vectors[, 1:k]
    return(specmat)
  })
  usumlist <- lapply(1:M, function(m) {
    return(ulist[[m]] %*% t(ulist[[m]]))
  })
  usum <- Reduce('+', usumlist)
  ustareigen <- RSpectra::eigs_sym(usum, k = k)
  ustar <- ustareigen$vectors[, 1:k]
  specmean <- kmeans(ustar, k)
  return(specmean$cluster)
}






##### The mean Adjacency matrix

meanadj <- function(x, n, k)
{
  adjsum <- Reduce("+", x)
  d = rowSums(adjsum)
  d = d + mean(d)
  deg <- diag(1 / sqrt(d))
  lap = deg %*% adjsum %*% deg
  spectra <- RSpectra::eigs_sym(adjsum, k = k)
  umean <- spectra$vectors[, 1:k]
  specmean <- kmeans(umean, k)
  return(specmean$cluster)
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


## The main function with BFGS optimization

lmfo <- function(x, n, k) {
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



##### Implementation of Co-OSNTF method from Paul and Chen, Annals of Applied Statistics, 2020  #####




conmffunction <- function(laplist, ulist, ustar, slist, alpha) {
  M = length(laplist)
  objlist <- lapply(1:M, function(m) {
    specobj <-
      (norm(laplist[[m]] - ulist[[m]] %*% slist[[m]] %*% t(ulist[[m]]), type =
              "F")) ^ 2
    regterm <- (norm(t(ulist[[m]]) %*% ustar, type = "F")) ^ 2
    return(specobj - alpha * regterm)
  })
  obj <- Reduce("+", objlist)
  return(obj)
}


# The main Co-OSNTF function

conmf <- function(x, n, k, alpha)
{
  M = length(x)
  laplist <- laplacian(x)
  uspectra <- lapply(1:M, function(m) {
    spectra <- eigen(laplist[[m]])
    return(spectra)
  })
  #unmf<-lapply(1:M,function(m){
  #  osntf<-osntfcal(laplist[[m]],n,k)
  #  return(osntf)
  #  })
  ulist <- lapply(1:M, function(m) {
    tau <- matrix(0.01, n, k)
    specmat <- uspectra[[m]]$vectors[, 1:k]
    rownorm <- apply(specmat, 1, function(a) {
      (sum(a ^ 2)) ^ 0.5
    })
    rownorm <- ifelse(rownorm < 10 ^ (-06), 10 ^ (-06), rownorm)
    specnorm <- specmat / rownorm
    specm <- kmeans(specnorm, k, nstart = 5)$cluster
    for (i in 1:k) {
      tau[, i] <- ifelse(specm == i, (1 - 0.01 * k), 0.01)
    }
    colnorm <- apply(tau, 2, function(a) {
      (sum(a ^ 2)) ^ 0.5
    })
    colnorm <- ifelse(colnorm < 10 ^ (-06), 10 ^ (-06), colnorm)
    tau <- t(t(tau) / colnorm)
    return(tau)
  })
  slist <- lapply(1:M, function(m) {
    s <- diag(abs(uspectra[[m]]$values[1:k]))
    return(s)
  })
  
  
  usumlist <- lapply(1:M, function(m) {
    return(ulist[[m]] %*% t(ulist[[m]]))
  })
  usum <- Reduce('+', usumlist) / M
  
  ustar <- matrix(0.01, n, k)
  specmat <- eigen(usum)$vectors[, 1:k]
  rownorm <- apply(specmat, 1, function(a) {
    (sum(a ^ 2)) ^ 0.5
  })
  rownorm <- ifelse(rownorm < 10 ^ (-06), 10 ^ (-06), rownorm)
  specnorm <- specmat / rownorm
  specm <- kmeans(specnorm, k, nstart = 5)$cluster
  for (i in 1:k) {
    ustar[, i] <- ifelse(specm == i, (1 - 0.01 * k), 0.01)
  }
  colnorm <- apply(ustar, 2, function(a) {
    (sum(a ^ 2)) ^ 0.5
  })
  colnorm <- ifelse(colnorm < 10 ^ (-06), 10 ^ (-06), colnorm)
  ustar <- t(t(ustar) / colnorm)
  #ustareigen <- eigen(usum)
  #ustar<- ustareigen$vectors[,1:k]
  ustarlast <- ustar
  t = 0
  valdiff = 1
  value <- conmffunction(laplist, ulist, ustar, slist, alpha)
  #print(value)
  while ((valdiff > 1e-04) & (t < 100))
  {
    valuelast = value
    #  ustarlast<-ustar
    ulistnew <- lapply(1:M, function(m) {
      x <-
        ulist[[m]] * sqrt((
          laplist[[m]] %*% ulist[[m]] %*% slist[[m]] + alpha * ustar
          %*% t(ustar) %*% ulist[[m]]
        ) / (
          ulist[[m]] %*% t(ulist[[m]]) %*% (
            laplist[[m]] %*% ulist[[m]] %*% slist[[m]] + alpha *
              ustar %*% t(ustar) %*% ulist[[m]]
          )
        ))
      x <- apply(x, c(1, 2), mycheck2)
      return(x)
    })
    ulist <- ulistnew
    ustarsum1 <- Reduce("+", lapply(1:M, function(m) {
      return(alpha * ulist[[m]] %*% t(ulist[[m]]) %*% ustar)
    }))
    ustarsum2 <- Reduce("+", lapply(1:M, function(m) {
      return(alpha * ustar
             %*% t(ustar) %*% ulist[[m]] %*% t(ulist[[m]]) %*% ustar)
    }))
    ustar <- ustar * sqrt(ustarsum1 / ustarsum2)
    ustar <- apply(ustar, c(1, 2), mycheck2)
    slistnew <- lapply(1:M, function(m) {
      x <-
        slist[[m]] * sqrt((t(ulist[[m]]) %*% laplist[[m]] %*% ulist[[m]]) / (t(ulist[[m]]) %*%
                                                                               ulist[[m]] %*% slist[[m]] %*%
                                                                               t(ulist[[m]]) %*% ulist[[m]]))
      x <- apply(x, c(1, 2), mycheck2)
      return(x)
    })
    slist <- slistnew
    
    value <- conmffunction(laplist, ulist, ustar, slist, alpha)
    
    valdiff = valuelast - value
    
    
    t = t + 1
    
  }
  nmfclus <- lapply(ulist, function(x) {
    nmfcomm <- apply(x, 1, which.max)
    return(nmfcomm)
  })
  nmfstar <- apply(ustar, 1, which.max)
  nmfclus[[M + 1]] <- nmfstar
  nmfclus[[M + 2]] <- Reduce("+", lapply(1:M, function(m) {
    y = clue::solve_LSAP(table(nmfclus[[m]], nmfclus[[M + 1]]), maximum = TRUE)
    for (i in 1:n) {
      nmfclus[[m]][i] <- y[nmfclus[[m]][i]]
    }
    return(table(nmfclus[[M + 1]], nmfclus[[m]]))
  }))
  nmfclus[[M + 2]] <- nmfclus[[M + 2]] / rowSums(nmfclus[[M + 2]])
  nmfclus[[M + 3]] <- list()
  for (m in 1:M) {
    clusters = nmfclus[[m]]
    A = x[[m]]
    ones = matrix(1, n, n)
    diag(ones) = 0
    phat = matrix(0, k, k)
    for (i in 1:k)
    {
      for (j in i:k)
      {
        phat[i, j] = sum(A[clusters == i, clusters == j]) / length(ones[clusters ==
                                                                          i, clusters == j])
        phat[j, i] = phat[i, j]
      }
    }
    nmfclus[[M + 3]][[m]] = phat
  }
  return(nmfclus)
}
