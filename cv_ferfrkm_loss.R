cv_ferfrkm_loss <- function(
  Xtr,
  G,
  Q,
  K = kspline(seq_len(ncol(Xtr)))$K,
  Pk = NULL,
  Lk = NULL,
  folds = 10,
  fold_ids = NULL, # vector of label assignment to folds
  gamma_grid = default_gamma_grid,
  lambda_grid = default_lambda_grid,
  max_iter = Inf,
  tol = 1e-8,
  nstart_kmeans = 10,
  parallel = FALSE,
  ncores = max(1L, parallel::detectCores() - 1L),
  seed = 123
) {
  stopifnot(is.matrix(Xtr), G >= 2, Q >= 1, Q <= min(G, ncol(Xtr)))
  #
  if (is.null(Pk) || is.null(Lk)) {
    ks <- kspline(seq_len(ncol(Xtr)))
    Pk <- ks$Pk
    Lk <- ks$Lk
  }
  #
  if (is.null(fold_ids)) {
    fold_ids <- make_folds(nrow(Xtr), k = folds, seed = seed)
  }
  folds <- sort(unique(fold_ids)) # fold labels
  #
  grid <- make_ferfrkm_grid(gamma_grid, lambda_grid)
  grid$id <- seq_len(nrow(grid))
  #
  fit_one_combo <- function(grid_row) {
    gamma <- grid_row$gamma
    lambda <- grid_row$lambda
    fold_scores <- numeric(length(folds))
#
    for (fold_idx in seq_along(folds)) {
      fold <- folds[fold_idx]
      train_idx <- fold_ids != fold
      valid_idx <- !train_idx
      X_train <- Xtr[train_idx, , drop = FALSE]
      X_valid <- Xtr[valid_idx, , drop = FALSE]
#
      init <- init_ferfrkm(
          X = X_train,
          G = G,
          Q = Q,
          seed = seed,
          nstart_kmeans = nstart_kmeans
        )
#
        fit <- FERFRKM(
          C = X_train,
          K = K,
          Pk = Pk,
          Lk = Lk,
          U = init$U,
          A = init$A,
          B = init$B,
          lambda = lambda,
          gamma = gamma,
          max_iter = max_iter,
          tol = tol
        )
#
      fold_scores[fold_idx] <- fit$loss_function_unpen
    }
#
    data.frame(
      gamma = gamma,
      lambda = lambda,
      mean_cv_score = mean(fold_scores),
      sd_cv_score = sd(fold_scores),
      stringsAsFactors = FALSE
    )
  }
#
  if (parallel) {
    cl <- parallel::makeCluster(ncores) # Creates a cluster with ncores
    on.exit(parallel::stopCluster(cl), add = TRUE) # shut down automatically when the function exits
    doParallel::registerDoParallel(cl) # %dopar% sends iterations to those workers.
#
    results <- foreach(
      grid_idx = seq_len(nrow(grid)),
      .combine = rbind,
      .packages = c("stats"), # Packages for each cluster
      .export = c( # Variables and function to consider
        "FERFRKM",
        "init_ferfrkm",
        "loss_function"
      )
    ) %dopar% {
      fit_one_combo(grid[grid_idx, , drop = FALSE])
    }
  } else {
    results <- do.call(
      rbind,
      lapply(seq_len(nrow(grid)), function(grid_idx) {
        fit_one_combo(grid[grid_idx, , drop = FALSE])
      })
    )
  }
  # sorts first by mean_cv_score, and uses sd_cv_score to break ties.
  results <- results[order(results$mean_cv_score, results$sd_cv_score), ]
  rownames(results) <- NULL
#
  list(
    results = results,
    best = results[1, , drop = FALSE]
  )
}