# @noRd
.hm_inv_marginal <- function(u, dist, params) {
  switch(dist,
    weibull  = stats::qweibull(u, params$shape, params$scale),
    gamma    = stats::qgamma(u, params$shape, rate = params$rate),
    lnorm    = stats::qlnorm(u, params$meanlog, params$sdlog),
    gumbel   = .qgumbel(u, params$mu, params$beta),
    tnorm    = .qtnorm(u, params$mean, params$sd),
    tlaplace = .qtlaplace(u, params$mu, params$b),
    stop("Unknown distribution: ", dist)
  )
}



#' Fit a spatial copula model to the uniform PIT data
#'
#' Fits a C-vine, D-vine, or Gaussian copula to the uniform PIT matrix
#' from [hm_pit_transform()], following Eq. 6 of Clavijo Mesa et al.
#' (2026). The three structures correspond to the comparison in Table IV
#' of that paper.
#'
#' Vine copulas are fitted via \code{rvinecopulib::vinecop()}. The
#' Gaussian copula uses the sin-transform
#' \eqn{\rho = \sin(\pi \hat{\tau} / 2)} to convert Kendall \eqn{\tau}
#' to a Pearson correlation matrix. Requires \code{rvinecopulib}
#' (listed in \code{Suggests}).
#'
#' @param U Numeric matrix (\eqn{T \times N}) of uniform PIT values from
#'   [hm_pit_transform()].
#' @param type Copula structure: \code{"cvine"} (default), \code{"dvine"},
#'   or \code{"gaussian"}.
#' @param ... Additional arguments passed to \code{rvinecopulib::vinecop()}.
#'
#' @return An \code{hm_copula} list with \code{type}, \code{model},
#'   \code{n_cells}, and \code{n_obs}.
#'
#' @references
#' Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026).
#' \emph{Reliability Engineering & System Safety}, 275, 112830.
#'
#' @seealso [hm_pit_transform()], [hm_simulate_copula()]
#'
#' @examples
#' \dontrun{
#' X_land <- hm_filter_land_cells(X)
#' fits   <- hm_fit_marginals(X_land)
#' pit    <- hm_pit_transform(X_land, fits)
#'
#' cop_c <- hm_fit_copula(pit$U, type = "cvine")
#' cop_d <- hm_fit_copula(pit$U, type = "dvine")
#' cop_g <- hm_fit_copula(pit$U, type = "gaussian")
#' }
#'
#' @export
hm_fit_copula <- function(U, type = "cvine", ...) {

  type <- match.arg(type, c("cvine", "dvine", "gaussian"))

  if (!is.matrix(U) || !is.numeric(U))
    stop("`U` must be a numeric matrix of uniform PIT values.")
  if (nrow(U) < 3L)
    stop("`U` must have at least 3 rows (extreme-event days).")
  if (ncol(U) < 2L)
    stop("`U` must have at least 2 columns (land cells).")
  if (any(U < 0 | U > 1, na.rm = TRUE))
    stop("`U` values must lie in [0, 1]. Run `hm_pit_transform()` first.")

  n_obs   <- nrow(U)
  n_cells <- ncol(U)

  if (type %in% c("cvine", "dvine")) {

    if (!requireNamespace("rvinecopulib", quietly = TRUE))
      stop("Package `rvinecopulib` is required for vine copulas.")

    vine_structure <- if (type == "cvine") {
      rvinecopulib::cvine_structure(seq_len(n_cells))
    } else {
      rvinecopulib::dvine_structure(seq_len(n_cells))
    }

    model <- rvinecopulib::vinecop(data = U, structure = vine_structure, ...)

  } else {

    Z       <- stats::qnorm(U)
    tau_mat <- stats::cor(Z, method = "kendall", use = "pairwise.complete.obs")
    R       <- .hm_nearPD(sin(pi / 2 * tau_mat))
    model   <- list(R = R, n_cells = n_cells)
  }

  structure(
    list(type = type, model = model, n_cells = n_cells, n_obs = n_obs),
    class = c("hm_copula", "list")
  )
}


# @noRd
.hm_nearPD <- function(R) {
  eig  <- eigen(R, symmetric = TRUE)
  lam  <- pmax(eig$values, 1e-8)
  R_pd <- eig$vectors %*% diag(lam) %*% t(eig$vectors)
  d    <- sqrt(diag(R_pd))
  R_pd <- R_pd / outer(d, d)
  diag(R_pd) <- 1
  R_pd
}



#' Simulate spatially coherent hazard scenarios from a fitted copula
#'
#' Draws \code{n_scenarios} realisations from a fitted copula and
#' back-transforms to physical units via the inverse marginal CDFs.
#' Implements Eq. 7--8 of Clavijo Mesa et al. (2026).
#'
#' Vine copulas are sampled via \code{rvinecopulib::rvinecop()}. The
#' Gaussian copula uses \code{MASS::mvrnorm()} followed by the standard
#' normal CDF. Both packages are listed in \code{Suggests}.
#'
#' @param copula An \code{hm_copula} object from [hm_fit_copula()].
#' @param fits data.frame from [hm_fit_marginals()] for the same land
#'   cells used to fit \code{copula}.
#' @param n_scenarios Number of scenarios. Defaults to 1000.
#' @param seed Optional integer random seed.
#'
#' @return A numeric matrix (\eqn{S \times N}) of simulated hazard
#'   intensities. Attributes \code{type}, \code{n_scenarios}, and
#'   \code{n_cells} are attached.
#'
#' @references
#' Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026).
#' \emph{Reliability Engineering & System Safety}, 275, 112830.
#'
#' @seealso [hm_fit_copula()], [hm_fit_marginals()]
#'
#' @examples
#' \dontrun{
#' cop       <- hm_fit_copula(pit$U, type = "cvine")
#' scenarios <- hm_simulate_copula(cop, fits, n_scenarios = 1000)
#' dim(scenarios)
#' range(scenarios, na.rm = TRUE)
#' }
#'
#' @export
hm_simulate_copula <- function(copula, fits, n_scenarios = 1000L,
                                seed = NULL) {

  if (!inherits(copula, "hm_copula"))
    stop("`copula` must be an object returned by `hm_fit_copula()`.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame returned by `hm_fit_marginals()`.")
  if (nrow(fits) != copula$n_cells)
    stop("`fits` must have ", copula$n_cells, " rows; got ", nrow(fits), ".")
  if (!is.numeric(n_scenarios) || length(n_scenarios) != 1L ||
      n_scenarios < 1L)
    stop("`n_scenarios` must be a positive integer.")

  n_scenarios <- as.integer(n_scenarios)
  N           <- copula$n_cells

  if (!is.null(seed))
    set.seed(seed)

  if (copula$type %in% c("cvine", "dvine")) {

    if (!requireNamespace("rvinecopulib", quietly = TRUE))
      stop("Package `rvinecopulib` is required for vine copula simulation.")

    U_sim <- as.matrix(rvinecopulib::rvinecop(n_scenarios, copula$model))

  } else {

    if (!requireNamespace("MASS", quietly = TRUE))
      stop("Package `MASS` is required for Gaussian copula simulation.")

    Z_sim <- MASS::mvrnorm(n_scenarios, mu = rep(0, N), Sigma = copula$model$R)
    U_sim <- stats::pnorm(Z_sim)
  }

  H_sim <- matrix(NA_real_, nrow = n_scenarios, ncol = N)

  for (j in seq_len(N)) {
    dist   <- fits$best_dist[j]
    params <- fits$params[[j]]
    if (is.na(dist) || is.null(params)) next
    H_sim[, j] <- tryCatch(
      .hm_inv_marginal(U_sim[, j], dist, params),
      error = function(e) rep(NA_real_, n_scenarios)
    )
  }

  attr(H_sim, "type")        <- copula$type
  attr(H_sim, "n_scenarios") <- n_scenarios
  attr(H_sim, "n_cells")     <- N

  H_sim
}
