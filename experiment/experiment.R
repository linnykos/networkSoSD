D <- matrix(0, nrow = m, ncol = k)
n <- nrow(D)
nidx <- apply(D, 1, which.min)
moved <- which(nidx != idx_prev)

# resolve tie in favor of not move
if (length(moved) > 0) {
  moved <- moved[D[(idx_prev[moved]-1)*n+moved] > min(D[moved,])] 
  idx[moved] <- nidx[moved]
} 

list(idx = idx, any_moved = (length(moved) == 0))


######

[(idx_prev[moved]-1)*n+moved]  
# vector of length length(moved)
# (idx_prev[moved]-1) pulls the values in idx_prev
# D is a matrix of n x k
# idx_prev is length n, but values go from 1 through k

D[(idx_prev[moved]-1)*n+moved]


min(D[moved,])
# min D[moved,] takes the minimum among all the rows in moved = 
### all the samples where the minimizing cluster isn't what was in the previous iteration