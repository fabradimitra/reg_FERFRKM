# Function Reference

Short reference for the functions defined in the project root.

## `FERFRKM(C, K, Pk, Lk, U, A, B, lambda, gamma, max_iter = Inf, tol = 1e-6)`
Runs the FERFRKM clustering algorithm and alternates updates for memberships, centroid coefficients, and the reduced basis.

- `C`: data matrix of spline coefficients, one row per sample.
- `K`: roughness penalty matrix for the spline basis.
- `Pk`: eigenvector matrix from the spline decomposition. Present in the API but not used inside the current implementation.
- `Lk`: diagonal matrix of spline eigenvalues. Present in the API but not used inside the current implementation.
- `U`: initial fuzzy membership matrix.
- `A`: initial coefficients for cluster centroids in the reduced subspace.
- `B`: initial coefficients of the reduced basis in the spline basis.
- `lambda`: smoothness penalty strength for the centroids.
- `gamma`: fuzziness parameter controlling how soft the memberships are.
- `max_iter`: maximum number of update iterations.
- `tol`: convergence tolerance on the penalized loss decrease.

## `CV_FERFRKM(Xtr, G, Q, K = kspline(seq_len(ncol(Xtr)))$K, Pk = NULL, Lk = NULL, lambda_init = 1, gamma_init = 1, folds = 10, fold_ids = NULL, max_iter = Inf, tol = 1e-8, nstart_kmeans = 10, seed = 123, randomstarts = 5)`
Tunes `lambda` and `gamma` by cross-validation using a Xie-Beni style score.

- `Xtr`: training data matrix.
- `G`: number of clusters.
- `Q`: reduced-rank dimension for the latent subspace.
- `K`: spline roughness penalty matrix.
- `Pk`: spline eigenvectors. If `NULL`, they are computed from `kspline()`.
- `Lk`: spline eigenvalues. If `NULL`, they are computed from `kspline()`.
- `lambda_init`: starting value for the optimizer over `lambda`.
- `gamma_init`: starting value for the optimizer over `gamma`.
- `folds`: number of cross-validation folds.
- `fold_ids`: optional fold assignment vector; if `NULL`, it is generated automatically.
- `max_iter`: maximum iterations for each FERFRKM fit inside CV.
- `tol`: convergence tolerance for each FERFRKM fit inside CV.
- `nstart_kmeans`: number of random starts for the initial k-means step.
- `seed`: random seed used for fold creation and initialization.
- `randomstarts`: number of FERFRKM restarts evaluated per fold.

### Internal helper: `cv_xie_beni_optim(vars, folds, fold_ids, G, Q, X, K, Pk, Lk, randomstarts, nstart_kmeans, seed, max_iter, tol)`
Evaluates the mean cross-validated Xie-Beni score for a candidate `(lambda, gamma)` pair. Function to be passed to optim.

- `vars`: numeric vector with `lambda` in `vars[1]` and `gamma` in `vars[2]`.
- `folds`: fold labels to evaluate.
- `fold_ids`: fold assignment for each sample.
- `G`: number of clusters.
- `Q`: reduced-rank dimension.
- `X`: data matrix used in cross-validation.
- `K`: spline roughness penalty matrix.
- `Pk`: spline eigenvectors.
- `Lk`: spline eigenvalues.
- `randomstarts`: number of initialization restarts per fold.
- `nstart_kmeans`: number of k-means starts for the first initialization.
- `seed`: random seed for initialization.
- `max_iter`: maximum number of FERFRKM iterations.
- `tol`: convergence tolerance for FERFRKM.

## `init_FERFRKM(X, G, Q, seed = NULL, nstart_kmeans = 10)`
Builds a deterministic starting point for FERFRKM using k-means memberships and an SVD-based reduced basis.

- `X`: data matrix used for initialization.
- `G`: number of clusters.
- `Q`: reduced-rank dimension.
- `seed`: optional random seed for reproducible initialization.
- `nstart_kmeans`: number of random starts passed to `kmeans()`.

## `kspline(t)`
Constructs the natural cubic spline penalty matrix and its spectral decomposition for a grid `t`.

- `t`: ordered vector of grid points where the spline basis is defined.

## `loss_function(U, C, Cbar, D, A, B, K, lambda, gamma)`
Computes the unpenalized and penalized FERFRKM loss function.

- `U`: fuzzy membership matrix.
- `C`: data matrix of spline coefficients.
- `Cbar`: cluster-centroid coefficient matrix in the spline basis.
- `D`: diagonal matrix built from the cluster membership weights.
- `A`: centroid coefficients in the reduced subspace.
- `B`: reduced basis coefficients in the spline basis.
- `K`: roughness penalty matrix.
- `lambda`: smoothness penalty strength.
- `gamma`: fuzziness parameter.

## `make_folds(n, k = 10, seed = 123)`
Creates balanced fold labels for cross-validation by randomly assigning observations to `k` folds.

- `n`: number of observations.
- `k`: number of folds.
- `seed`: random seed used before sampling fold labels.

## `perm_hungarian_fast(M_true, M_est, J)`
Finds the permutation of estimated cluster centroids that best matches the true centroids using the Hungarian algorithm.

- `M_true`: matrix of true centroids.
- `M_est`: matrix of estimated centroids.
- `J`: number of grid points or columns used to normalize the squared-distance cost.

## `randgenuc(n, k)`
Generates a random hard membership matrix with one active cluster per observation.

- `n`: number of observations.
- `k`: number of clusters.

## `randgenuf(n, k)`
Generates a random fuzzy membership matrix with positive entries normalized by column.

- `n`: number of observations.
- `k`: number of clusters.

## `rand_orthogonal(G, Q)`
Generates a random `G x Q` orthonormal matrix using QR decomposition.

- `G`: number of rows.
- `Q`: number of orthogonal columns to generate.

## Script-local helpers in `sim12_kfoldCV.R`
Replication of the simulation in the parent paper of 2012 by Gattone and Rocci. 
