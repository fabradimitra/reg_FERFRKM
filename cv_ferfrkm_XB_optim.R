cv_ferfrkm_XB_optim <- function(
  Xtr,
  G,
  Q,
  K = kspline(seq_len(ncol(Xtr)))$K,
  Pk = NULL,
  Lk = NULL,
  folds = 10,
  fold_ids = NULL, 
  max_iter = Inf,
  tol = 1e-8,
  nstart_kmeans = 10,
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
  xie_beni_optim <- function(vars, A_init, B_init, U_init, X_train, X_valid, K, Pk, Lk, max_iter, tol){
    fit <- FERFRKM(
          C = X_train,
          K = K,
          Pk = Pk,
          Lk = Lk,
          U = U_init,
          A = A_init,
          B = B_init,
          lambda = vars[1],
          gamma = vars[2],
          max_iter = max_iter,
          tol = tol
          )
    centers <- fit$A %*% t(fit$B)   # G x J
    cnorm2_valid <- rowSums(X_valid^2)     # length n_valid
    vnorm2 <- rowSums(centers^2)           # length G
    dist2_valid <- outer(cnorm2_valid, vnorm2, "+") - 2 * (X_valid %*% t(centers))
    U_valid <- exp(-dist2_valid / vars[2])
    U_valid <- U_valid / rowSums(U_valid)
    U_valid[U_valid < 1e-12] <- 1e-12
    sep2 <- as.matrix(dist(centers))^2
    sep2[sep2 == 0] <- Inf
    min_sep2 <- min(sep2)
    sum(U_valid * dist2_valid) / (nrow(X_valid) * min_sep2)
  }
  fold_scores <- numeric(length(folds))
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
      res_optim <- tryCatch(
        optim(
          par = c(1,1),
          fn = xie_beni_optim,
          method = "Nelder-Mead",
          A_init = init$A,
          B_init = init$B,
          U_init = init$U,
          X_train = X_train,
          X_valid = X_valid,
          K = K,
          Pk = Pk,
          Lk = Lk,
          max_iter = max_iter,
          tol = tol
        ),
        error = function(e) NULL
      )
      if (is.null(res_optim)) {
          next
      }
      cur_score <- res_optim$value
      if(cur_score<score_fin){
          lambda <- res_optim$par[1]
          gamma <- res_optim$par[2]
          score_fin <- cur_score
      }
    }
      fold_scores[fold_idx] <- score_fin
  }
  list(
    gamma = gamma,
    lambda = lambda,
    mean_cv_score = mean(fold_scores),
    sd_cv_score = sd(fold_scores),
    stringsAsFactors = FALSE
  )
}
