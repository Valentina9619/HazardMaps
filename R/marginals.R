# Supported distribution families and their fitdistrplus names.
# Zero-truncated variants are handled internally via truncated MLE.
.hm_supported_dists <- c("weibull", "gamma", "lnorm", "gumbel",
                          "tnorm", "tlaplace")

# Internal: display-friendly names for messages and output columns.
.hm_dist_labels <- c(
  weibull  = "Weibull",
  gamma    = "Gamma",
  lnorm    = "Lognormal",
  gumbel   = "Gumbel",
  tnorm    = "Zero-trunc. Gaussian",
  tlaplace = "Zero-trunc. Laplace"
)



# ------------------------------------------------------------------------------
# Internal distribution helpers
# ------------------------------------------------------------------------------

# Density, CDF and random generation for the Gumbel distribution.
# parametrized by location (mu) and scale (beta > 0).
.dgumbel <- function(x, mu, beta, log = FALSE) {
  if (beta <= 0) return(rep(if (log) -Inf else 0, length(x)))
  z   <- (x - mu) / beta
  val <- -log(beta) - z - exp(-z)
  if (log) val else exp(val)
}

.pgumbel <- function(q, mu, beta, lower.tail = TRUE) {
  p <- exp(-exp(-(q - mu) / beta))
  if (lower.tail) p else 1 - p
}

.qgumbel <- function(p, mu, beta) {
  mu - beta * log(-log(p))
}


# Density, CDF and quantile for the zero-truncated Gaussian (x > 0).
.dtnorm <- function(x, mean, sd, log = FALSE) {
  if (sd <= 0) return(rep(if (log) -Inf else 0, length(x)))
  phi0 <- stats::pnorm(0, mean, sd, lower.tail = FALSE)
  val  <- stats::dnorm(x, mean, sd, log = FALSE) / phi0
  val[x <= 0] <- 0
  if (log) log(val) else val
}

.ptnorm <- function(q, mean, sd, lower.tail = TRUE) {
  phi0 <- stats::pnorm(0, mean, sd, lower.tail = FALSE)
  p    <- pmax(0, (stats::pnorm(q, mean, sd) -
                   stats::pnorm(0, mean, sd)) / phi0)
  p[q <= 0] <- 0
  if (lower.tail) p else 1 - p
}

.qtnorm <- function(p, mean, sd) {
  phi0 <- stats::pnorm(0, mean, sd, lower.tail = FALSE)
  stats::qnorm(stats::pnorm(0, mean, sd) + p * phi0, mean, sd)
}


# Density, CDF and quantile for the zero-truncated Laplace (x > 0).
# parametrized by location (mu) and scale (b > 0).
.dtlaplace <- function(x, mu, b, log = FALSE) {
  if (b <= 0) return(rep(if (log) -Inf else 0, length(x)))
  # Laplace CDF at 0: normalisation constant
  F0  <- if (mu >= 0) 0.5 * exp(-mu / b) else 1 - 0.5 * exp(mu / b)
  val <- (1 / (2 * b)) * exp(-abs(x - mu) / b) / (1 - F0)
  val[x <= 0] <- 0
  if (log) log(val) else val
}

.ptlaplace <- function(q, mu, b, lower.tail = TRUE) {
  .plaplace <- function(x, mu, b) {
    ifelse(x < mu,
           0.5 * exp((x - mu) / b),
           1 - 0.5 * exp(-(x - mu) / b))
  }
  F0 <- .plaplace(0, mu, b)
  p  <- pmax(0, (.plaplace(q, mu, b) - F0) / (1 - F0))
  p[q <= 0] <- 0
  if (lower.tail) p else 1 - p
}

.qtlaplace <- function(p, mu, b) {
  .qlaplace <- function(p, mu, b) {
    ifelse(p < 0.5,
           mu + b * log(2 * p),
           mu - b * log(2 * (1 - p)))
  }
  .plaplace <- function(x, mu, b) {
    ifelse(x < mu,
           0.5 * exp((x - mu) / b),
           1 - 0.5 * exp(-(x - mu) / b))
  }
  F0 <- .plaplace(0, mu, b)
  .qlaplace(F0 + p * (1 - F0), mu, b)
}



# ------------------------------------------------------------------------------
# Internal: fit one distribution to one cell's data
# ------------------------------------------------------------------------------

# Returns a named list with: dist, params, loglik, aic, bic, converged.
# Returns NULL on failure so the caller can skip gracefully.
.hm_fit_one <- function(x, dist) {

  n <- length(x)

  tryCatch({

    if (dist == "weibull") {
      fit    <- fitdistrplus::fitdist(x, "weibull")
      params <- as.list(fit$estimate)

    } else if (dist == "gamma") {
      fit    <- fitdistrplus::fitdist(x, "gamma")
      params <- as.list(fit$estimate)

    } else if (dist == "lnorm") {
      fit    <- fitdistrplus::fitdist(x, "lnorm")
      params <- as.list(fit$estimate)

    } else if (dist == "gumbel") {
      # MLE via numerical optimisation
      nll <- function(par) {
        mu <- par[1]; beta <- par[2]
        if (beta <= 0) return(1e15)
        -sum(.dgumbel(x, mu, beta, log = TRUE))
      }
      opt    <- stats::optim(c(mean(x), stats::sd(x)), nll,
                             method = "Nelder-Mead")
      if (opt$convergence != 0) return(NULL)
      params <- list(mu = opt$par[1], beta = opt$par[2])

    } else if (dist == "tnorm") {
      nll <- function(par) {
        m <- par[1]; s <- par[2]
        if (s <= 0) return(1e15)
        -sum(.dtnorm(x, m, s, log = TRUE))
      }
      opt    <- stats::optim(c(mean(x), stats::sd(x)), nll,
                             method = "Nelder-Mead")
      if (opt$convergence != 0) return(NULL)
      params <- list(mean = opt$par[1], sd = opt$par[2])

    } else if (dist == "tlaplace") {
      nll <- function(par) {
        mu <- par[1]; b <- par[2]
        if (b <= 0) return(1e15)
        -sum(.dtlaplace(x, mu, b, log = TRUE))
      }
      opt    <- stats::optim(c(mean(x), stats::mad(x) + 1e-6), nll,
                             method = "Nelder-Mead")
      if (opt$convergence != 0) return(NULL)
      params <- list(mu = opt$par[1], b = opt$par[2])

    } else {
      return(NULL)
    }

    # Log-likelihood at the MLE
    loglik <- if (dist %in% c("weibull", "gamma", "lnorm")) {
      fit$loglik
    } else {
      -opt$value
    }

    k   <- length(params)
    aic <- -2 * loglik + 2 * k
    bic <- -2 * loglik + log(n) * k

    list(dist      = dist,
         params    = params,
         loglik    = loglik,
         aic       = aic,
         bic       = bic,
         converged = TRUE)

  }, error = function(e) NULL)
}



# ------------------------------------------------------------------------------
# Public functions
# ------------------------------------------------------------------------------

#' Fit marginal distributions to each grid cell
#'
#' For each grid cell (column) in the hazard matrix, fits a set of candidate
#' parametric distributions by maximum likelihood and selects the best-fitting
#' one using AIC or BIC.
#'
#' @details
#' Fitting is performed on the event-focused sample produced by
#' [hm_select_extreme_events()], where each row is an extreme-event day and
#' each column is a grid cell. This matches the methodology in Clavijo Mesa
#' et al. (2026), where marginal distributions are calibrated independently
#' at each spatial location to capture site-specific hazard behaviour.
#'
#' Six candidate families are available:
#' \itemize{
#'   \item \code{"weibull"}: two-parameter Weibull (scale, shape). Dominant
#'     in wind-gust fields.
#'   \item \code{"gamma"}: two-parameter Gamma (shape, rate).
#'   \item \code{"lnorm"}: Lognormal (meanlog, sdlog). Common in pre-Alpine
#'     locations.
#'   \item \code{"gumbel"}: Gumbel extreme-value (location mu, scale beta).
#'   \item \code{"tnorm"}: zero-truncated Gaussian (mean, sd). Enforces
#'     non-negative support.
#'   \item \code{"tlaplace"}: zero-truncated Laplace (location mu, scale b).
#' }
#'
#' Parameters are estimated by maximum likelihood. \code{"weibull"},
#' \code{"gamma"} and \code{"lnorm"} use [fitdistrplus::fitdist()];
#' \code{"gumbel"}, \code{"tnorm"} and \code{"tlaplace"} use
#' [stats::optim()] with Nelder-Mead.
#'
#' The best-fitting distribution at each cell is selected by minimising either
#' AIC or BIC (set with \code{criterion}). Both criteria are always computed
#' and stored in the output. Cells for which all candidate fits fail are
#' recorded with \code{NA} in the \code{best_dist} column.
#'
#' Requires the \code{fitdistrplus} package (listed in \code{Suggests}).
#'
#' @param X A numeric matrix of dimension \eqn{T \times N} (time steps by
#'   grid cells), as returned by [hm_to_matrix()] applied to the output of
#'   [hm_select_extreme_events()]. All values must be strictly positive.
#' @param dists Character vector of candidate distribution families to fit.
#'   Any subset of \code{c("weibull", "gamma", "lnorm", "gumbel", "tnorm",
#'   "tlaplace")}. Defaults to all six.
#' @param criterion Model-selection criterion. Either \code{"bic"}
#'   (default, used in Clavijo Mesa et al. 2026 SF paper) or \code{"aic"}.
#'
#' @return A data.frame with one row per grid cell and columns:
#'   \describe{
#'     \item{\code{cell_id}}{Integer cell index (column number of \code{X}).}
#'     \item{\code{best_dist}}{Name of the selected distribution family.}
#'     \item{\code{best_dist_label}}{Human-readable distribution name.}
#'     \item{\code{params}}{Named list of MLE parameters for the best fit.}
#'     \item{\code{loglik}}{Log-likelihood of the best fit.}
#'     \item{\code{aic}}{AIC of the best fit.}
#'     \item{\code{bic}}{BIC of the best fit.}
#'     \item{\code{n_obs}}{Number of non-\code{NA} observations in the cell.}
#'     \item{\code{n_failed}}{Number of candidate fits that failed to converge.}
#'   }
#'   The returned data.frame also carries an attribute \code{"fits_full"}: a
#'   list of length \eqn{N}, where each element is a list of raw fit results
#'   (one per candidate) for that cell. Use it with [hm_marginal_gof()].
#'
#' @seealso [hm_marginal_gof()] for goodness-of-fit diagnostics,
#'   [hm_pit_transform()] to apply the fitted CDFs.
#'
#' @references
#' Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). A modeling framework
#' for the inoperability assessment of interdependent critical infrastructures
#' exposed to spatially distributed natural hazards.
#' \emph{Reliability Engineering & System Safety}, 275, 112830.
#'
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' Inoperability assessment of interdependent critical infrastructures
#' exposed to natural hazards considering climate change.
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @examples
#' \dontrun{
#' f <- system.file("extdata", "Copernicus_dailymax_1994_2021.nc",
#'                  package = "HazardMaps")
#' x <- hm_read_netcdf(f, var = "i10fg")
#' x <- hm_standardize_coords(x)
#' x <- hm_decode_time(x)
#' x <- hm_select_extreme_events(x, threshold = 25)
#' X <- hm_to_matrix(x)
#'
#' fits <- hm_fit_marginals(X)
#' head(fits[, c("cell_id", "best_dist", "aic", "bic")])
#' table(fits$best_dist)
#' }
#'
#' @export
hm_fit_marginals <- function(X,
                              dists     = c("weibull", "gamma", "lnorm",
                                            "gumbel", "tnorm", "tlaplace"),
                              criterion = "bic") {

  if (!requireNamespace("fitdistrplus", quietly = TRUE))
    stop("Package `fitdistrplus` is required. ",
         "Install it with install.packages(\"fitdistrplus\").")

  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (nrow(X) < 2L)
    stop("`X` must have at least 2 rows (time steps).")

  criterion <- match.arg(criterion, c("bic", "aic"))

  bad_dists <- setdiff(dists, .hm_supported_dists)
  if (length(bad_dists) > 0L)
    stop("Unknown distribution(s): ", paste(bad_dists, collapse = ", "),
         ". Supported: ", paste(.hm_supported_dists, collapse = ", "))

  N <- ncol(X)

  # Pre-allocate output vectors
  best_dist  <- character(N)
  best_label <- character(N)
  params_out <- vector("list", N)
  loglik_out <- numeric(N)
  aic_out    <- numeric(N)
  bic_out    <- numeric(N)
  n_obs      <- integer(N)
  n_failed   <- integer(N)
  fits_full  <- vector("list", N)

  for (j in seq_len(N)) {

    xj <- X[, j]
    xj <- xj[is.finite(xj) & xj > 0]   # drop NA and non-positive values
    n_obs[j] <- length(xj)

    if (n_obs[j] < 5L) {
      # Not enough data to fit: leave as NA
      best_dist[j]  <- NA_character_
      best_label[j] <- NA_character_
      params_out[j] <- list(NULL)
      loglik_out[j] <- NA_real_
      aic_out[j]    <- NA_real_
      bic_out[j]    <- NA_real_
      n_failed[j]   <- length(dists)
      fits_full[[j]] <- list()
      next
    }

    # Fit all candidate distributions
    cell_fits  <- lapply(dists, .hm_fit_one, x = xj)
    names(cell_fits) <- dists
    fits_full[[j]]   <- cell_fits

    # Keep only the ones that converged
    ok <- Filter(Negate(is.null), cell_fits)
    n_failed[j] <- length(dists) - length(ok)

    if (length(ok) == 0L) {
      best_dist[j]  <- NA_character_
      best_label[j] <- NA_character_
      params_out[j] <- list(NULL)
      loglik_out[j] <- NA_real_
      aic_out[j]    <- NA_real_
      bic_out[j]    <- NA_real_
      next
    }

    # Select best by chosen criterion
    scores <- sapply(ok, function(f) f[[criterion]])
    best   <- ok[[which.min(scores)]]

    best_dist[j]  <- best$dist
    best_label[j] <- .hm_dist_labels[best$dist]
    params_out[j] <- list(best$params)
    loglik_out[j] <- best$loglik
    aic_out[j]    <- best$aic
    bic_out[j]    <- best$bic
  }

  out <- data.frame(
    cell_id        = seq_len(N),
    best_dist      = best_dist,
    best_dist_label = best_label,
    loglik         = loglik_out,
    aic            = aic_out,
    bic            = bic_out,
    n_obs          = n_obs,
    n_failed       = n_failed,
    stringsAsFactors = FALSE
  )

  # Store MLE params as a list-column
  out$params <- params_out

  # Attach full per-cell fits for use in hm_marginal_gof()
  attr(out, "fits_full") <- fits_full
  attr(out, "criterion") <- criterion

  out
}



#' Goodness-of-fit diagnostics for a single grid cell
#'
#' Computes the Kolmogorov-Smirnov statistic and returns the data needed to
#' draw a QQ plot for the best-fitting marginal distribution at a chosen cell.
#'
#' @details
#' The KS statistic is provided as a **complementary diagnostic**, not as a
#' model-selection criterion. As discussed in Clavijo Mesa et al. (2026),
#' the best-fitting distribution is selected by BIC or AIC; the KS test
#' and QQ plot are used afterwards to confirm that the selected model
#' provides an adequate representation of the data (see Section 3.2 of the
#' SF paper and Appendix B of the Copernicus paper).
#'
#' Note that the KS p-value reported here is computed with fitted (not
#' fixed) parameters, so it should be interpreted with caution: the test
#' is anti-conservative under MLE-fitted parameters, meaning p-values tend
#' to be larger than their nominal level. The statistic itself remains a
#' useful measure of the maximum deviation between the empirical and fitted
#' CDFs.
#'
#' @param X A numeric matrix (same one passed to [hm_fit_marginals()]).
#' @param fits The data.frame returned by [hm_fit_marginals()].
#' @param cell_id Integer. The cell to diagnose (column index of \code{X}).
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{cell_id}}{The cell examined.}
#'     \item{\code{dist}}{Selected distribution name.}
#'     \item{\code{dist_label}}{Human-readable distribution name.}
#'     \item{\code{params}}{MLE parameters of the selected fit.}
#'     \item{\code{ks_statistic}}{KS test statistic D.}
#'     \item{\code{ks_pvalue}}{KS p-value (interpret with caution; see Details).}
#'     \item{\code{qq}}{data.frame with \code{theoretical} and
#'       \code{empirical} quantiles for QQ plotting.}
#'     \item{\code{n_obs}}{Number of observations used.}
#'   }
#'
#' @seealso [hm_fit_marginals()], [hm_pit_transform()]
#'
#' @examples
#' \dontrun{
#' fits <- hm_fit_marginals(X)
#' gof  <- hm_marginal_gof(X, fits, cell_id = 1)
#'
#' # QQ plot
#' plot(gof$qq$theoretical, gof$qq$empirical,
#'      xlab = "Theoretical quantiles", ylab = "Empirical quantiles",
#'      main = paste("QQ plot —", gof$dist_label, "| cell", gof$cell_id))
#' abline(0, 1, col = "firebrick")
#' }
#'
#' @export
hm_marginal_gof <- function(X, fits, cell_id) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame returned by `hm_fit_marginals()`.")
  if (!is.numeric(cell_id) || length(cell_id) != 1L ||
      cell_id < 1L || cell_id > ncol(X))
    stop("`cell_id` must be a single integer between 1 and ncol(X).")

  cell_id <- as.integer(cell_id)
  row     <- fits[fits$cell_id == cell_id, , drop = FALSE]

  if (nrow(row) == 0L)
    stop("No entry found in `fits` for cell_id = ", cell_id)

  dist   <- row$best_dist
  params <- row$params[[1L]]

  if (is.na(dist) || is.null(params))
    stop("No successful fit available for cell_id = ", cell_id,
         ". Check `fits$n_failed`.")

  # Observed values for this cell
  xj <- X[, cell_id]
  xj <- sort(xj[is.finite(xj) & xj > 0])
  n  <- length(xj)

  # CDF and quantile functions for the selected distribution
  pfun <- switch(dist,
    weibull  = function(q) stats::pweibull(q, params$shape, params$scale),
    gamma    = function(q) stats::pgamma(q, params$shape, rate = params$rate),
    lnorm    = function(q) stats::plnorm(q, params$meanlog, params$sdlog),
    gumbel   = function(q) .pgumbel(q, params$mu, params$beta),
    tnorm    = function(q) .ptnorm(q, params$mean, params$sd),
    tlaplace = function(q) .ptlaplace(q, params$mu, params$b),
    stop("Unknown distribution: ", dist)
  )

  qfun <- switch(dist,
    weibull  = function(p) stats::qweibull(p, params$shape, params$scale),
    gamma    = function(p) stats::qgamma(p, params$shape, rate = params$rate),
    lnorm    = function(p) stats::qlnorm(p, params$meanlog, params$sdlog),
    gumbel   = function(p) .qgumbel(p, params$mu, params$beta),
    tnorm    = function(p) .qtnorm(p, params$mean, params$sd),
    tlaplace = function(p) .qtlaplace(p, params$mu, params$b),
    stop("Unknown distribution: ", dist)
  )

  # KS statistic: max |F_n(x) - F(x)|
  ecdf_vals <- seq_len(n) / n
  fitted_cdf <- pfun(xj)
  ks_stat    <- max(abs(ecdf_vals - fitted_cdf))

  # Approximate p-value via ks.test.
  # The "ties" warning is suppressed because wind-gust data from gridded
  # climate datasets often contains repeated values. The p-value is
  # slightly anti-conservative in that case, which is already documented
  # in the @details section. The KS statistic itself is unaffected.
  ks_pval <- tryCatch({
    suppressWarnings(stats::ks.test(xj, pfun)$p.value)
  }, error = function(e) NA_real_)

  # QQ data: theoretical vs empirical quantiles
  probs <- (seq_len(n) - 0.5) / n
  qq    <- data.frame(
    theoretical = qfun(probs),
    empirical   = xj
  )

  list(
    cell_id      = cell_id,
    dist         = dist,
    dist_label   = .hm_dist_labels[dist],
    params       = params,
    ks_statistic = ks_stat,
    ks_pvalue    = ks_pval,
    qq           = qq,
    n_obs        = n
  )
}



#' Apply the Probability Integral Transform to all grid cells
#'
#' Transforms the hazard matrix to uniform margins by applying each cell's
#' fitted CDF, then optionally applies Gaussian anamorphosis (probit
#' transform) to obtain standardised Gaussian scores.
#'
#' @details
#' This function implements Equations (3) and (4) from Clavijo Mesa et al.
#' (2026, SF paper):
#'
#' \deqn{u_l(t) = F_l(H^*(d_l, t)) \in [0, 1]}
#'
#' \deqn{Z_l(t) = \Phi^{-1}(u_l(t))}
#'
#' The uniform matrix \eqn{U} is the input to copula fitting
#' (\code{fit_copula.R}). The Gaussian scores matrix \eqn{Z} is the input
#' to spatial dependence estimation and KLE-based stochastic field generation
#' (\code{fit_SF.R}).
#'
#' Uniform values are clipped to \eqn{[10^{-6}, 1 - 10^{-6}]} before the
#' probit step to avoid infinite Gaussian scores at the boundaries.
#'
#' @param X A numeric matrix of dimension \eqn{T \times N}, as returned by
#'   [hm_to_matrix()] applied to the output of [hm_select_extreme_events()].
#' @param fits The data.frame returned by [hm_fit_marginals()].
#' @param gaussian Logical. If \code{TRUE} (default), also returns the
#'   Gaussian scores matrix \eqn{Z = \Phi^{-1}(U)}, centred by subtracting
#'   the empirical column means (as in Eq. 4 of the SF paper).
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{U}}{Numeric matrix \eqn{T \times N} of uniform PIT values
#'       in \eqn{[0, 1]}.}
#'     \item{\code{Z}}{Numeric matrix \eqn{T \times N} of centred Gaussian
#'       scores. \code{NULL} when \code{gaussian = FALSE}.}
#'     \item{\code{n_cells_ok}}{Number of cells successfully transformed.}
#'     \item{\code{n_cells_skipped}}{Number of cells skipped due to missing
#'       fit (recorded as \code{NA} columns in \code{U} and \code{Z}).}
#'   }
#'
#' @seealso [hm_fit_marginals()], [hm_marginal_gof()]
#'
#' @examples
#' \dontrun{
#' fits <- hm_fit_marginals(X)
#' pit  <- hm_pit_transform(X, fits)
#'
#' # Uniform matrix for copula fitting
#' dim(pit$U)
#'
#' # Gaussian scores matrix for stochastic field
#' dim(pit$Z)
#' }
#'
#' @export
hm_pit_transform <- function(X, fits, gaussian = TRUE) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame returned by `hm_fit_marginals()`.")
  if (ncol(X) != nrow(fits))
    stop("`ncol(X)` must equal `nrow(fits)`. ",
         "Make sure `X` and `fits` correspond to the same dataset.")

  N <- ncol(X)
  T_steps <- nrow(X)

  U <- matrix(NA_real_, nrow = T_steps, ncol = N)

  n_ok      <- 0L
  n_skipped <- 0L

  for (j in seq_len(N)) {

    row    <- fits[j, , drop = FALSE]
    dist   <- row$best_dist
    params <- row$params[[1L]]

    # Skip cells with no successful fit
    if (is.na(dist) || is.null(params)) {
      n_skipped <- n_skipped + 1L
      next
    }

    # Build the CDF function for the selected distribution.
    # Note: switch() with next inside a for loop does not behave as expected
    # in R, so unknown distributions are caught with an explicit check instead.
    if (!dist %in% .hm_supported_dists) {
      n_skipped <- n_skipped + 1L
      next
    }

    xj <- X[, j]

    pfun <- switch(dist,
      weibull  = function(q) stats::pweibull(q, params$shape, params$scale),
      gamma    = function(q) stats::pgamma(q, params$shape, rate = params$rate),
      lnorm    = function(q) stats::plnorm(q, params$meanlog, params$sdlog),
      gumbel   = function(q) .pgumbel(q, params$mu, params$beta),
      tnorm    = function(q) .ptnorm(q, params$mean, params$sd),
      tlaplace = function(q) .ptlaplace(q, params$mu, params$b)
    )

    U[, j] <- pfun(xj)
    n_ok   <- n_ok + 1L
  }

  # Clip to avoid infinite values in the probit step.
  # pmax/pmin drop the matrix class when NA values are present,
  # so dimensions are restored explicitly afterwards.
  eps <- 1e-6
  U[] <- pmax(eps, pmin(1 - eps, U))

  # Gaussian scores: probit transform + centering (Eq. 4 in SF paper)
  Z <- NULL
  if (isTRUE(gaussian)) {
    Z <- stats::qnorm(U)
    # Centre each column by subtracting its empirical mean.
    # apply() is used instead of colMeans() to safely handle NA columns
    # without dropping matrix dimensions.
    col_means <- apply(Z, 2, mean, na.rm = TRUE)
    Z <- sweep(Z, 2, col_means, FUN = "-")
  }

  list(
    U               = U,
    Z               = Z,
    n_cells_ok      = n_ok,
    n_cells_skipped = n_skipped
  )
}
