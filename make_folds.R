make_folds <- function(n, k = 10, seed = 123) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}