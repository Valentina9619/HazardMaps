# @noRd
.hm_matern <- function(D, theta) {
  range      <- theta[1]
  smoothness <- theta[2]
  nugget     <- theta[3]
  h          <- D / range
  corr       <- matrix(0, nrow = nrow(D), ncol = ncol(D))
  zero       <- h == 0
  corr[zero] <- 1
  pos        <- !zero & is.finite(h)
  if (any(pos)) {
    hv         <- h[pos]
    corr[pos]  <- (2^(1 - smoothness) / gamma(smoothness)) *
                  (sqrt(2 * smoothness) * hv)^smoothness *
                  besselK(sqrt(2 * smoothness) * hv, smoothness)
  }
  corr       <- (1 - nugget) * corr
  diag(corr) <- 1
  corr
}

# @noRd
.hm_exponential <- function(D, theta) {
  corr       <- exp(-D / theta[1])
  corr       <- (1 - theta[2]) * corr
  diag(corr) <- 1
  corr
}

# @noRd
.hm_gaussian_cov <- function(D, theta) {
  corr       <- exp(-(D / theta[1])^2)
  corr       <- (1 - theta[2]) * corr
  diag(corr) <- 1
  corr
}

# @noRd
.hm_spherical <- function(D, theta) {
  h          <- D / theta[1]
  corr       <- ifelse(h <= 1, 1 - 1.5 * h + 0.5 * h^3, 0)
  corr       <- (1 - theta[2]) * corr
  diag(corr) <- 1
  corr
}


# Internal: weighted least-squares fit for one covariance family.
# Returns a list with params, rmse, and converged.
.hm_fit_kernel <- function(fam, corg, range_init) {

  w       <- corg$n_pairs / sum(corg$n_pairs)
  d_row   <- matrix(corg$bin_center, nrow = 1)

  kernel  <- switch(fam,
    matern      = .hm_matern,
    exponential = .hm_exponential,
    gaussian    = .hm_gaussian_cov,
    spherical   = .hm_spherical
  )

  is_matern <- fam == "matern"

  obj <- function(p) {
    valid <- if (is_matern) p[1] > 0 && p[2] > 0 && p[3] >= 0 && p[3] < 1
             else           p[1] > 0 && p[2] >= 0 && p[2] < 1
    if (!valid) return(1e10)
    pred <- kernel(d_row, p)[1, ]
    sum(w * (corg$mean_tau - pred)^2)
  }

  init <- if (is_matern) c(range_init, 0.5, 0.1) else c(range_init, 0.1)

  tryCatch({
    opt    <- stats::optim(init, obj, method = "Nelder-Mead",
                           control = list(maxit = 2000))
    if (opt$convergence != 0) stop("did not converge")
    pred   <- kernel(d_row, opt$par)[1, ]
    params <- if (is_matern) {
      c(range = opt$par[1], smoothness = opt$par[2], nugget = opt$par[3])
    } else {
      c(range = opt$par[1], nugget = opt$par[2])
    }
    rmse   <- sqrt(sum(w * (corg$mean_tau - pred)^2))
    list(family = fam, params = params, rmse = rmse, converged = TRUE)
  }, error = function(e) {
    list(family = fam, params = NA, rmse = Inf, converged = FALSE)
  })
}



#' Fit a parametric spatial correlation model to the empirical correlogram
#'
#' Fits Matérn, Exponential, Gaussian, and/or Spherical covariance families
#' to the correlogram from [hm_empirical_correlogram()] by weighted least
#' squares (Eq. 6--7 of Clavijo Mesa et al., 2026, SF paper). Selects the
#' best family by RMSE, matching Table 1 of that paper.
#'
#' @param corg data.frame from [hm_empirical_correlogram()].
#' @param families Covariance families to fit. Any subset of
#'   \code{c("matern", "exponential", "gaussian", "spherical")}.
#'   Defaults to all four.
#' @param max_dist Maximum distance in km to include. Defaults to the
#'   full correlogram range.
#'
#' @return An \code{hm_corrmodel} list with \code{best_family},
#'   \code{best_params}, \code{best_rmse}, \code{all_fits}, and
#'   \code{corg}.
#'
#' @references
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @seealso [hm_empirical_correlogram()], [hm_simulate_sf_kle()]
#'
#' @examples
#' \dontrun{
#' corg  <- hm_empirical_correlogram(tau, D)
#' model <- hm_fit_correlation_model(corg)
#' model$best_family
#' model$all_fits
#' }
#'
#' @export
hm_fit_correlation_model <- function(corg,
                                      families = c("matern", "exponential",
                                                   "gaussian", "spherical"),
                                      max_dist = NULL) {

  if (!is.data.frame(corg))
    stop("`corg` must be a data.frame from `hm_empirical_correlogram()`.")
  if (!all(c("bin_center", "mean_tau", "n_pairs") %in% names(corg)))
    stop("`corg` must contain `bin_center`, `mean_tau`, and `n_pairs`.")

  families <- match.arg(families,
                        c("matern", "exponential", "gaussian", "spherical"),
                        several.ok = TRUE)

  if (!is.null(max_dist)) {
    corg <- corg[corg$bin_center <= max_dist, , drop = FALSE]
    if (nrow(corg) == 0L)
      stop("No correlogram bins within `max_dist` = ", max_dist, " km.")
  }

  range_init <- corg$bin_center[which(corg$mean_tau < 0.5)[1L]]
  if (is.na(range_init))
    range_init <- max(corg$bin_center) / 2

  results    <- lapply(families, .hm_fit_kernel, corg = corg,
                       range_init = range_init)
  names(results) <- families

  rmse_vals <- sapply(results, function(r) r$rmse)
  best_fam  <- names(which.min(rmse_vals))
  best_fit  <- results[[best_fam]]

  all_fits <- do.call(rbind, lapply(results, function(r) {
    data.frame(family    = r$family,
               rmse      = round(r$rmse, 6),
               converged = r$converged,
               stringsAsFactors = FALSE)
  }))
  rownames(all_fits) <- NULL

  structure(
    list(best_family = best_fam,
         best_params = best_fit$params,
         best_rmse   = best_fit$rmse,
         all_fits    = all_fits,
         corg        = corg),
    class = c("hm_corrmodel", "list")
  )
}



#' Simulate hazard scenarios via Karhunen-Loève Expansion
#'
#' Generates \code{n_scenarios} synthetic hazard fields from the KLE of
#' the fitted covariance model, following Eq. 8--11 of Clavijo Mesa et al.
#' (2026, SF paper): build \eqn{\bar{C}}, eigen-decompose and truncate to
#' \code{var_retain} variance, sample KLE coefficients, back-transform
#' via \eqn{\Phi} and the inverse marginal CDFs.
#'
#' @param corrmodel An \code{hm_corrmodel} from [hm_fit_correlation_model()].
#' @param D Symmetric distance matrix in km from [hm_distance_matrix()].
#'   Must have \code{nrow(fits)} rows.
#' @param fits data.frame from [hm_fit_marginals()] for the same land cells.
#' @param n_scenarios Number of synthetic scenarios. Defaults to 1000.
#' @param var_retain Variance fraction retained in KLE truncation.
#'   Defaults to 0.95.
#' @param seed Optional integer random seed.
#'
#' @return A numeric matrix (\eqn{S \times N}) of simulated hazard
#'   intensities. Attributes \code{family}, \code{k_modes},
#'   \code{var_retain}, and \code{var_explained} are attached.
#'
#' @references
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @seealso [hm_fit_correlation_model()], [hm_fit_marginals()]
#'
#' @examples
#' \dontrun{
#' model     <- hm_fit_correlation_model(corg)
#' scenarios <- hm_simulate_sf_kle(model, D, fits,
#'                                  n_scenarios = 1000, seed = 42)
#' attr(scenarios, "k_modes")
#' attr(scenarios, "var_explained")
#' }
#'
#' @export
hm_simulate_sf_kle <- function(corrmodel, D, fits,
                                n_scenarios = 1000L,
                                var_retain  = 0.95,
                                seed        = NULL) {

  if (!inherits(corrmodel, "hm_corrmodel"))
    stop("`corrmodel` must be from `hm_fit_correlation_model()`.")
  if (!is.matrix(D) || !is.numeric(D))
    stop("`D` must be a numeric matrix from `hm_distance_matrix()`.")
  if (!isSymmetric(D))
    stop("`D` must be symmetric.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame from `hm_fit_marginals()`.")
  if (nrow(fits) != nrow(D))
    stop("`fits` and `D` must have the same number of rows.")
  if (!is.numeric(var_retain) || length(var_retain) != 1L ||
      var_retain <= 0 || var_retain > 1)
    stop("`var_retain` must be in (0, 1].")
  if (!is.numeric(n_scenarios) || n_scenarios < 1L)
    stop("`n_scenarios` must be a positive integer.")

  n_scenarios <- as.integer(n_scenarios)
  N           <- nrow(D)

  if (!is.null(seed))
    set.seed(seed)

  fam    <- corrmodel$best_family
  params <- corrmodel$best_params

  C_bar <- switch(fam,
    matern      = .hm_matern(D, params),
    exponential = .hm_exponential(D, params),
    gaussian    = .hm_gaussian_cov(D, params),
    spherical   = .hm_spherical(D, params),
    stop("Unknown covariance family: ", fam)
  )

  eig     <- eigen(C_bar, symmetric = TRUE)
  lam     <- pmax(eig$values, 0)
  cum_var <- cumsum(lam) / sum(lam)
  k       <- which(cum_var >= var_retain)[1L]
  if (is.na(k)) k <- N

  var_explained <- cum_var[k]

  message("KLE: retaining ", k, " of ", N, " eigenmodes (",
          round(100 * var_explained, 1), "% variance explained).")

  V_k   <- eig$vectors[, seq_len(k), drop = FALSE]
  lam_k <- lam[seq_len(k)]

  xi    <- matrix(stats::rnorm(k * n_scenarios), nrow = k, ncol = n_scenarios)
  Z_sim <- t(V_k %*% (sqrt(lam_k) * xi))
  U_sim <- stats::pnorm(Z_sim)

  H_sim <- matrix(NA_real_, nrow = n_scenarios, ncol = N)

  for (j in seq_len(N)) {
    dist     <- fits$best_dist[j]
    params_j <- fits$params[[j]]
    if (is.na(dist) || is.null(params_j)) next
    H_sim[, j] <- tryCatch(
      .hm_inv_marginal(U_sim[, j], dist, params_j),
      error = function(e) rep(NA_real_, n_scenarios)
    )
  }

  attr(H_sim, "family")        <- fam
  attr(H_sim, "k_modes")       <- k
  attr(H_sim, "var_retain")    <- var_retain
  attr(H_sim, "var_explained") <- var_explained

  H_sim
}
