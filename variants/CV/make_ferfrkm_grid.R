make_ferfrkm_grid <- function(
  gamma_grid = default_gamma_grid,
  lambda_grid = default_lambda_grid
) {
  expand.grid(
    gamma = gamma_grid,
    lambda = lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}