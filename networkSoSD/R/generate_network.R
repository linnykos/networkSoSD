.generate_membership_matrix <- function(membership_vector){
  K <- max(membership_vector)
  
  sapply(1:K, function(x){
    as.numeric(membership_vector == x)
  })
}

#' Compute probability matrix
#'
#' @param B_mat symmetric connectivity matrix
#' @param membership_vec vector containing values \code{1} through \code{ncol(B_mat)}
#'
#' @return symmetric matrix of dimension \code{length(membership_vec)}
#' @export
compute_prob_mat <- function(B_mat, membership_vec){
  membership_mat <- .generate_membership_matrix(membership_vec)
  prob_mat <- membership_mat %*% B_mat %*% t(membership_mat)
  
  prob_mat
}

#' Simulate adjacency matrix
#'
#' @param prob_mat symmetric probability matrix
#'
#' @return symmetric adjacency matrix
#' @export
generate_adjaceny_mat <- function(prob_mat){
  upper_tri_idx <- upper.tri(prob_mat)
  prob_upper <- prob_mat[upper_tri_idx]
  
  n <- nrow(prob_mat)
  adj_upper <- stats::rbinom(n*(n-1)/2, 1, prob_upper)
  adj_mat <- matrix(0, ncol = n, nrow = n)
  adj_mat[upper_tri_idx] <- adj_upper
  adj_mat <- adj_mat + t(adj_mat)
  
  adj_mat
}

.l2norm <- function(x){sqrt(sum(x^2))}
