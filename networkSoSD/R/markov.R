#' Compute the appropriate transition matrix to form desired stationary vector
#'
#' Authors note: This does NOT always work. To fix this, I suspect
#' I would need to reparameterize the problem such that the D matrix
#' is positive definite (by having only 3 unknowns instead of 6).
#'
#' @param vec desired stationary vector of length 3 
#' @param delta number between 0 and 1 (inclusive)
#' @param sig_dig number of significant digits
#'
#' @return a row-stochastic matrix
#' @export
compute_markov_transition <- function(vec, delta, sig_dig = 3){
  stopifnot(length(vec) == 3, sum(vec) == 1, all(vec > 0),
            delta >= 0, delta <= 1)
  if(delta == 0){
    return(diag(3))
  }
  
  val1 <- vec[1]; val2 <- vec[2]; val3 <- vec[3]
  n <- 6 # number of unknowns
  
  D <- matrix(0, n, n)
  diag(D) <- rep(c(val1^2, val2^2, val3^2), each = 2)
  D[1,6] <- D[6,1] <- val1*val3
  D[2,4] <- D[4,2] <- val1*val2
  D[3,5] <- D[5,3] <- val2*val3
  D <- 2*D
  
  tmp <- c(val1*val2, val1*val3, 
           val1*val2, val2*val3, 
           val1*val3, val2*val3)
  d <- 2*tmp - 2*tmp*(1-delta)
  
  A <- matrix(0, nrow = 3, ncol = n)
  A[1,1:2] <- 1
  A[2,3:4] <- 1
  A[3,5:6] <- 1
  
  b0 <- rep(delta, 3)
 
  res <- kernlab::ipop(c = matrix(-d, ncol = 1), 
                       H = D, A = A, 
                       b = matrix(b0, ncol = 1), 
                       l = matrix(rep(0,6), ncol = 1), 
                       u = matrix(rep(1,6), ncol = 1), 
                       r = matrix(rep(0, 3), ncol = 1),
                       maxiter = 1000)
  
  P <- matrix(0, 3, 3)
  diag(P) <- 1-delta
  P[1,2:3] <- res@primal[1:2]
  P[2,c(1,3)] <- res@primal[3:4]
  P[3,1:2] <- res@primal[5:6]
  
  round(P, sig_dig)
}