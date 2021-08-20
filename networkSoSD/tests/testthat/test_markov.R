context("Test markov")

## compute_markov_transition is correct

test_that("compute_markov_transition works", {
  vec <- c(0.5, 0.4, 0.1)
  res <- compute_markov_transition(vec, delta = 0.1)
  expect_true(all(res >= 0) & all(dim(res) == 3) & all(abs(rowSums(res)-1) <= 1e-3))
  expect_true(sum(abs(vec %*% res - vec)) <= 1e-3)
  
  res <- compute_markov_transition(vec, delta = 0)
  expect_true(all(res >= 0) & all(dim(res) == 3) & all(abs(rowSums(res)-1) <= 1e-3))
  expect_true(sum(abs(vec %*% res - vec)) <= 1e-3)
  
  res <- compute_markov_transition(vec, delta = 1)
  expect_true(all(res >= 0) & all(dim(res) == 3) & all(abs(rowSums(res)-1) <= 1e-3))
  expect_true(sum(abs(vec %*% res - vec)) <= 1e-3)
})

# test_that("compute_markov_transition works for random inputs", {
#   trials <- 100
#   
#   bool_vec <- sapply(1:trials, function(x){
#     set.seed(x)
#     vec <- runif(3); vec <- vec/sum(vec)
#     delta <- runif(1)
#     
#     res <- compute_markov_transition(vec, delta)
#     
#     sum(abs(vec %*% res - vec)) <= 1e-3
#   })
#   
#   expect_true(all(bool_vec))
# })


