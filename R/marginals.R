.hm_supported_dists <- c("weibull", "gamma", "lnorm", "gumbel",
                          "tnorm", "tlaplace")

.hm_dist_labels <- c(
  weibull  = "Weibull",
  gamma    = "Gamma",
  lnorm    = "Lognormal",
  gumbel   = "Gumbel",
  tnorm    = "Zero-trunc. Gaussian",
  tlaplace = "Zero-trunc. Laplace"
)


# --- Gumbel distribution -----------------------------------------------------

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

.qgumbel <- function(p, mu, beta) mu - beta * log(-log(p))


# --- Zero-truncated Gaussian (x > 0) -----------------------------------------

.dtnorm <- function(x, mean, sd, log = FALSE) {
  if (sd <= 0) return(rep(if (log) -Inf else 0, length(x)))
  phi0        <- stats::pnorm(0, mean, sd, lower.tail = FALSE)
  val         <- stats::dnorm(x, mean, sd, log = FALSE) / phi0
  val[x <= 0] <- 0
  if (log) log(val) else val
}

.ptnorm <- function(q, mean, sd, lower.tail = TRUE) {
  phi0      <- stats::pnorm(0, mean, sd, lower.tail = FALSE)
  p         <- pmax(0, (stats::pnorm(q, mean, sd) -
                        stats::pnorm(0, mean, sd)) / phi0)
  p[q <= 0] <- 0
  if (lower.tail) p else 1 - p
}

.qtnorm <- function(p, mean, sd) {
  phi0 <- stats::pnorm(0, mean, sd, lower.tail = FALSE)
  stats::qnorm(stats::pnorm(0, mean, sd) + p * phi0, mean, sd)
}


# --- Zero-truncated Laplace (x > 0) ------------------------------------------

.plaplace <- function(x, mu, b) {
  ifelse(x < mu,
         0.5 * exp((x - mu) / b),
         1 - 0.5 * exp(-(x - mu) / b))
}

.dtlaplace <- function(x, mu, b, log = FALSE) {
  if (b <= 0) return(rep(if (log) -Inf else 0, length(x)))
  F0          <- .plaplace(0, mu, b)
  val         <- (1 / (2 * b)) * exp(-abs(x - mu) / b) / (1 - F0)
  val[x <= 0] <- 0
  if (log) log(val) else val
}

.ptlaplace <- function(q, mu, b, lower.tail = TRUE) {
  F0        <- .plaplace(0, mu, b)
  p         <- pmax(0, (.plaplace(q, mu, b) - F0) / (1 - F0))
  p[q <= 0] <- 0
  if (lower.tail) p else 1 - p
}

.qtlaplace <- function(p, mu, b) {
  .qlaplace <- function(p, mu, b) {
    ifelse(p < 0.5,
           mu + b * log(2 * p),
           mu - b * log(2 * (1 - p)))
  }
  F0 <- .plaplace(0, mu, b)
  .qlaplace(F0 + p * (1 - F0), mu, b)
}


# --- CDF / quantile dispatch --------------------------------------------------

.hm_pfun <- function(dist, params) {
  switch(dist,
    weibull  = function(q) stats::pweibull(q, params$shape, params$scale),
    gamma    = function(q) stats::pgamma(q, params$shape, rate = params$rate),
    lnorm    = function(q) stats::plnorm(q, params$meanlog, params$sdlog),
    gumbel   = function(q) .pgumbel(q, params$mu, params$beta),
    tnorm    = function(q) .ptnorm(q, params$mean, params$sd),
    tlaplace = function(q) .ptlaplace(q, params$mu, params$b),
    stop("Unknown distribution: ", dist)
  )
}

.hm_qfun <- function(dist, params) {
  switch(dist,
    weibull  = function(p) stats::qweibull(p, params$shape, params$scale),
    gamma    = function(p) stats::qgamma(p, params$shape, rate = params$rate),
    lnorm    = function(p) stats::qlnorm(p, params$meanlog, params$sdlog),
    gumbel   = function(p) .qgumbel(p, params$mu, params$beta),
    tnorm    = function(p) .qtnorm(p, params$mean, params$sd),
    tlaplace = function(p) .qtlaplace(p, params$mu, params$b),
    stop("Unknown distribution: ", dist)
  )
}


# --- MLE fitting for one cell ------------------------------------------------

.hm_fit_one <- function(x, dist) {

  n <- length(x)

  tryCatch({

    if (dist %in% c("weibull", "gamma", "lnorm")) {
      fit    <- fitdistrplus::fitdist(x, dist)
      params <- as.list(fit$estimate)
      loglik <- fit$loglik

    } else {
      nll <- switch(dist,
        gumbel   = function(par) {
          if (par[2] <= 0) return(1e15)
          -sum(.dgumbel(x, par[1], par[2], log = TRUE))
        },
        tnorm    = function(par) {
          if (par[2] <= 0) return(1e15)
          -sum(.dtnorm(x, par[1], par[2], log = TRUE))
        },
        tlaplace = function(par) {
          if (par[2] <= 0) return(1e15)
          -sum(.dtlaplace(x, par[1], par[2], log = TRUE))
        }
      )

      init <- c(mean(x), stats::sd(x) + 1e-6)
      opt  <- stats::optim(init, nll, method = "Nelder-Mead")
      if (opt$convergence != 0) return(NULL)

      params <- switch(dist,
        gumbel   = list(mu = opt$par[1], beta = opt$par[2]),
        tnorm    = list(mean = opt$par[1], sd = opt$par[2]),
        tlaplace = list(mu = opt$par[1], b = opt$par[2])
      )
      loglik <- -opt$value
    }

    k <- length(params)
    list(
      dist      = dist,
      params    = params,
      loglik    = loglik,
      aic       = -2 * loglik + 2 * k,
      bic       = -2 * loglik + log(n) * k,
      converged = TRUE
    )

  }, error = function(e) NULL)
}



#' Fit marginal distributions to each grid cell
#'
#' Fits candidate parametric distributions to each column of \code{X} by
#' maximum likelihood and selects the best fit using AIC or BIC, following
#' the methodology in Clavijo Mesa et al. (2026).
#'
#' @param X A numeric matrix (\eqn{T \times N}) from [hm_to_matrix()].
#'   All values must be strictly positive.
#' @param dists Candidate families. Any subset of
#'   \code{c("weibull", "gamma", "lnorm", "gumbel", "tnorm", "tlaplace")}.
#' @param criterion Selection criterion: \code{"bic"} (default) or
#'   \code{"aic"}.
#'
#' @return A data.frame with one row per cell and columns \code{cell_id},
#'   \code{best_dist}, \code{best_dist_label}, \code{params}, \code{loglik},
#'   \code{aic}, \code{bic}, \code{n_obs}, \code{n_failed}. The attribute
#'   \code{"fits_full"} stores raw per-cell fit results for use in
#'   [hm_marginal_gof()].
#'
#' @seealso [hm_marginal_gof()], [hm_pit_transform()]
#' @export
hm_fit_marginals <- function(X,
                              dists     = .hm_supported_dists,
                              criterion = "bic") {

  if (!requireNamespace("fitdistrplus", quietly = TRUE))
    stop("Package `fitdistrplus` is required. ",
         "Install with install.packages(\"fitdistrplus\").")
  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (nrow(X) < 2L)
    stop("`X` must have at least 2 rows.")

  criterion <- match.arg(criterion, c("bic", "aic"))

  bad <- setdiff(dists, .hm_supported_dists)
  if (length(bad) > 0L)
    stop("Unknown distribution(s): ", paste(bad, collapse = ", "))

  N <- ncol(X)

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

    xj       <- X[, j]
    xj       <- xj[is.finite(xj) & xj > 0]
    n_obs[j] <- length(xj)

    if (n_obs[j] < 5L) {
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

    cell_fits        <- lapply(dists, .hm_fit_one, x = xj)
    names(cell_fits) <- dists
    fits_full[[j]]   <- cell_fits

    ok       <- Filter(Negate(is.null), cell_fits)
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
    cell_id         = seq_len(N),
    best_dist       = best_dist,
    best_dist_label = best_label,
    loglik          = loglik_out,
    aic             = aic_out,
    bic             = bic_out,
    n_obs           = n_obs,
    n_failed        = n_failed,
    stringsAsFactors = FALSE
  )

  out$params            <- params_out
  attr(out, "fits_full") <- fits_full
  attr(out, "criterion") <- criterion
  out
}



#' Goodness-of-fit diagnostics for a single grid cell
#'
#' Returns the KS statistic and QQ data for the best-fitting distribution at
#' a chosen cell. The KS test is a complementary diagnostic — distribution
#' selection is done by AIC/BIC in [hm_fit_marginals()].
#'
#' @param X Numeric matrix passed to [hm_fit_marginals()].
#' @param fits data.frame from [hm_fit_marginals()].
#' @param cell_id Column index to diagnose.
#'
#' @return A list with \code{cell_id}, \code{dist}, \code{dist_label},
#'   \code{params}, \code{ks_statistic}, \code{ks_pvalue}, \code{qq},
#'   and \code{n_obs}.
#'
#' @seealso [hm_fit_marginals()], [hm_pit_transform()]
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
    stop("No entry in `fits` for cell_id = ", cell_id)

  dist   <- row$best_dist
  params <- row$params[[1L]]

  if (is.na(dist) || is.null(params))
    stop("No successful fit for cell_id = ", cell_id,
         ". Check `fits$n_failed`.")

  xj <- sort(X[, cell_id][is.finite(X[, cell_id]) & X[, cell_id] > 0])
  n  <- length(xj)

  pfun <- .hm_pfun(dist, params)
  qfun <- .hm_qfun(dist, params)

  ecdf_vals  <- seq_len(n) / n
  ks_stat    <- max(abs(ecdf_vals - pfun(xj)))
  ks_pval    <- tryCatch(
    suppressWarnings(stats::ks.test(xj, pfun)$p.value),
    error = function(e) NA_real_
  )

  probs <- (seq_len(n) - 0.5) / n
  qq    <- data.frame(theoretical = qfun(probs), empirical = xj)

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



#' Apply the Probability Integral Transform
#'
#' Transforms \code{X} to uniform margins using each cell's fitted CDF, then
#' optionally applies the probit transform to obtain centred Gaussian scores.
#' Implements Equations (3)--(4) from Clavijo Mesa et al. (2026, SF paper).
#'
#' @param X Numeric matrix (\eqn{T \times N}) from [hm_to_matrix()].
#' @param fits data.frame from [hm_fit_marginals()].
#' @param gaussian If \code{TRUE} (default), also returns centred Gaussian
#'   scores \eqn{Z = \Phi^{-1}(U)}.
#'
#' @return A list with \code{U} (uniform matrix), \code{Z} (Gaussian scores
#'   or \code{NULL}), \code{n_cells_ok}, and \code{n_cells_skipped}.
#'
#' @seealso [hm_fit_marginals()], [hm_marginal_gof()]
#' @export
hm_pit_transform <- function(X, fits, gaussian = TRUE) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame returned by `hm_fit_marginals()`.")
  if (ncol(X) != nrow(fits))
    stop("`ncol(X)` (", ncol(X), ") must equal `nrow(fits)` (", nrow(fits), ").")

  N       <- ncol(X)
  U       <- matrix(NA_real_, nrow = nrow(X), ncol = N)
  n_ok      <- 0L
  n_skipped <- 0L

  for (j in seq_len(N)) {
    dist   <- fits$best_dist[j]
    params <- fits$params[[j]]

    if (is.na(dist) || is.null(params) || !dist %in% .hm_supported_dists) {
      n_skipped <- n_skipped + 1L
      next
    }

    U[, j] <- .hm_pfun(dist, params)(X[, j])
    n_ok   <- n_ok + 1L
  }

  eps  <- 1e-6
  U[]  <- pmax(eps, pmin(1 - eps, U))

  Z <- NULL
  if (isTRUE(gaussian)) {
    Z <- stats::qnorm(U)
    Z <- sweep(Z, 2, apply(Z, 2, mean, na.rm = TRUE), FUN = "-")
  }

  list(U = U, Z = Z, n_cells_ok = n_ok, n_cells_skipped = n_skipped)
}
