cv_ferfrkm_XB <- function(
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
  seed = 123,
  randomstarts = 5
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
      score_fin <- Inf
      for(start in seq_len(randomstarts)){
        if(start == 1){
          init <- init_ferfrkm(
          X = X_train,
          G = G,
          Q = Q,
          seed = seed,
          nstart_kmeans = nstart_kmeans
          )
        }else{
          U_init <- randgenuc(nrow(X_train), G)
          A_init <- rand_orthogonal(G, Q)
          B_init <- t(t(A_init) %*% solve(t(U_init) %*% U_init) %*% t(U_init) %*% X_train)
          init <- list(U = U_init, A = A_init, B = B_init)
        }
#
        fit <- tryCatch(
          FERFRKM(
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
          ),
          error = function(e) NULL
        )
      if (is.null(fit)) {
          next
        }
      # Xie-Beni index on the left out fold if successful fit
      centers <- fit$A %*% t(fit$B)   # G x J
      cnorm2_valid <- rowSums(X_valid^2)     # length n_valid
      vnorm2 <- rowSums(centers^2)           # length G
      dist2_valid <- outer(cnorm2_valid, vnorm2, "+") - 2 * (X_valid %*% t(centers))
      U_valid <- exp(-dist2_valid / gamma)
      U_valid <- U_valid / rowSums(U_valid)
      U_valid[U_valid < 1e-12] <- 1e-12
      sep2 <- as.matrix(dist(centers))^2
      sep2[sep2 == 0] <- Inf
      min_sep2 <- min(sep2)
      cur_score <- sum(U_valid * dist2_valid) / (nrow(X_valid) * min_sep2)
        if(cur_score<score_fin){
          score_fin <- cur_score
        }
      }
      fold_scores[fold_idx] <- score_fin
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
        "loss_function",
        "randgenuc",
        "rand_orthogonal"
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
