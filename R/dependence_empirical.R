#' Compute pairwise haversine distances between grid cells
#'
#' Returns a symmetric \eqn{N \times N} distance matrix of great-circle
#' distances in kilometres between all pairs of land cells (Eq. 5 of
#' Clavijo Mesa et al., 2026, SF paper).
#'
#' @param coords A data.frame with \code{lat} and \code{lon} columns in
#'   decimal degrees. Typically \code{attr(X_land, "coords")} from
#'   [hm_filter_land_cells()].
#'
#' @return A symmetric numeric matrix with diagonal zero.
#'
#' @references
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @seealso [hm_empirical_correlogram()], [hm_filter_land_cells()]
#'
#' @examples
#' \dontrun{
#' X_land <- hm_filter_land_cells(X)
#' D      <- hm_distance_matrix(attr(X_land, "coords"))
#' range(D)
#' }
#'
#' @export
hm_distance_matrix <- function(coords) {

  if (!is.data.frame(coords))
    stop("`coords` must be a data.frame with `lat` and `lon` columns.")
  if (!all(c("lat", "lon") %in% names(coords)))
    stop("`coords` must contain `lat` and `lon` columns.")

  lat <- coords$lat
  lon <- coords$lon

  if (any(!is.finite(lat)) || any(!is.finite(lon)))
    stop("`lat` and `lon` must be finite numeric values.")

  N  <- nrow(coords)
  rE <- 6371

  phi <- lat * pi / 180
  psi <- lon * pi / 180

  phi_i <- matrix(phi, nrow = N, ncol = N, byrow = FALSE)
  phi_j <- matrix(phi, nrow = N, ncol = N, byrow = TRUE)
  psi_i <- matrix(psi, nrow = N, ncol = N, byrow = FALSE)
  psi_j <- matrix(psi, nrow = N, ncol = N, byrow = TRUE)

  a <- sin((phi_j - phi_i) / 2)^2 +
       cos(phi_i) * cos(phi_j) * sin((psi_j - psi_i) / 2)^2

  D <- matrix(2 * rE * asin(pmin(1, sqrt(a))), nrow = N, ncol = N)
  D[cbind(seq_len(N), seq_len(N))] <- 0

  if (!is.null(coords$cell_id))
    rownames(D) <- colnames(D) <- as.character(coords$cell_id)

  D
}



#' Compute the empirical pairwise Kendall correlation matrix
#'
#' Estimates Kendall's \eqn{\hat{\tau}} between every pair of grid cells
#' from the Gaussian scores matrix \code{Z} from [hm_pit_transform()].
#' Implements Eq. 20 of Clavijo Mesa et al. (2026, Copernicus paper).
#'
#' For large grids (N > 1000) supply a subsample via \code{cell_idx} to
#' keep computation feasible.
#'
#' @param Z Numeric matrix (\eqn{T \times N}) of centred Gaussian scores
#'   from [hm_pit_transform()].
#' @param cell_idx Optional integer vector of column indices. Defaults to
#'   all columns.
#'
#' @return A symmetric \eqn{N \times N} matrix of Kendall \eqn{\hat{\tau}}
#'   values in \eqn{[-1, 1]}.
#'
#' @references
#' Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026).
#' \emph{Reliability Engineering & System Safety}, 275, 112830.
#'
#' @seealso [hm_pit_transform()], [hm_empirical_correlogram()]
#'
#' @examples
#' \dontrun{
#' fits <- hm_fit_marginals(X)
#' pit  <- hm_pit_transform(X, fits)
#' tau  <- hm_empirical_kendall(pit$Z, cell_idx = 1:30)
#' range(tau)
#' }
#'
#' @export
hm_empirical_kendall <- function(Z, cell_idx = NULL) {

  if (!is.matrix(Z) || !is.numeric(Z))
    stop("`Z` must be a numeric matrix from `hm_pit_transform()`.")

  if (!is.null(cell_idx)) {
    if (!is.numeric(cell_idx) || any(cell_idx < 1) ||
        any(cell_idx > ncol(Z)))
      stop("`cell_idx` must be integer indices between 1 and ncol(Z).")
    Z <- Z[, as.integer(cell_idx), drop = FALSE]
  }

  N <- ncol(Z)

  if (N < 2L)
    stop("`Z` must have at least 2 columns.")

  tau <- matrix(1, nrow = N, ncol = N)

  for (i in seq_len(N - 1L)) {
    for (j in seq(i + 1L, N)) {
      xi <- Z[, i]
      xj <- Z[, j]
      ok <- is.finite(xi) & is.finite(xj)
      if (sum(ok) < 3L) {
        tau[i, j] <- tau[j, i] <- NA_real_
        next
      }
      tau[i, j] <- tau[j, i] <- stats::cor(xi[ok], xj[ok],
                                            method = "kendall")
    }
  }

  tau
}



#' Compute the empirical spatial correlogram
#'
#' Bins pairwise Kendall correlations by inter-site distance to give
#' average correlation as a function of distance. The result is the key
#' input for fitting a covariance model in [hm_fit_correlation_model()],
#' following Eq. 5--6 of Clavijo Mesa et al. (2026, SF paper).
#'
#' @param tau Symmetric numeric matrix of pairwise Kendall correlations
#'   from [hm_empirical_kendall()].
#' @param D Symmetric numeric distance matrix in km from
#'   [hm_distance_matrix()]. Must have the same dimensions as \code{tau}.
#' @param n_bins Number of equal-width distance bins. Defaults to 30.
#' @param max_dist Maximum distance in km to include. Defaults to the
#'   observed maximum.
#'
#' @return A data.frame with columns \code{bin_center}, \code{mean_tau},
#'   \code{n_pairs}, \code{bin_lo}, and \code{bin_hi}.
#'
#' @references
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @seealso [hm_empirical_kendall()], [hm_distance_matrix()],
#'   [hm_fit_correlation_model()]
#'
#' @examples
#' \dontrun{
#' tau  <- hm_empirical_kendall(pit$Z, cell_idx = 1:50)
#' D    <- hm_distance_matrix(coords[1:50, ])
#' corg <- hm_empirical_correlogram(tau, D, n_bins = 20)
#' head(corg)
#' }
#'
#' @export
hm_empirical_correlogram <- function(tau, D, n_bins = 30L, max_dist = NULL) {

  if (!is.matrix(tau) || !is.numeric(tau))
    stop("`tau` must be a numeric matrix from `hm_empirical_kendall()`.")
  if (!is.matrix(D) || !is.numeric(D))
    stop("`D` must be a numeric matrix from `hm_distance_matrix()`.")
  if (!identical(dim(tau), dim(D)))
    stop("`tau` and `D` must have the same dimensions.")
  if (!is.numeric(n_bins) || length(n_bins) != 1L || n_bins < 2L)
    stop("`n_bins` must be an integer >= 2.")

  n_bins <- as.integer(n_bins)

  upper  <- upper.tri(D)
  d_vals <- D[upper]
  t_vals <- tau[upper]

  ok     <- is.finite(d_vals) & is.finite(t_vals)
  d_vals <- d_vals[ok]
  t_vals <- t_vals[ok]

  if (is.null(max_dist))
    max_dist <- max(d_vals)

  in_range <- d_vals <= max_dist
  d_vals   <- d_vals[in_range]
  t_vals   <- t_vals[in_range]

  if (length(d_vals) == 0L)
    stop("No valid pairs within `max_dist` = ", max_dist, " km.")

  breaks     <- seq(0, max_dist, length.out = n_bins + 1L)
  bin_lo     <- breaks[-length(breaks)]
  bin_hi     <- breaks[-1L]
  bin_center <- (bin_lo + bin_hi) / 2

  mean_tau <- numeric(n_bins)
  n_pairs  <- integer(n_bins)

  for (b in seq_len(n_bins)) {
    in_bin      <- d_vals > bin_lo[b] & d_vals <= bin_hi[b]
    n_pairs[b]  <- sum(in_bin)
    mean_tau[b] <- if (n_pairs[b] > 0L)
                     mean(t_vals[in_bin], na.rm = TRUE)
                   else
                     NA_real_
  }

  keep <- n_pairs > 0L
  data.frame(
    bin_center       = bin_center[keep],
    mean_tau         = mean_tau[keep],
    n_pairs          = n_pairs[keep],
    bin_lo           = bin_lo[keep],
    bin_hi           = bin_hi[keep],
    stringsAsFactors = FALSE
  )
}



#' Compute spatial cluster sizes of extreme events
#'
#' For each time step, counts the land cells where the hazard value meets
#' or exceeds \code{threshold}. This is the cluster-size metric
#' \eqn{B(t)} (Eq. 21, Appendix C of Clavijo Mesa et al., 2026) used to
#' validate whether a dependence model reproduces the observed spatial
#' footprint.
#'
#' @param X Numeric matrix (\eqn{T \times N}) from [hm_filter_land_cells()].
#' @param threshold A single finite numeric value.
#' @param time Optional date vector of length \eqn{T}. Defaults to row
#'   indices.
#'
#' @return A data.frame with columns \code{date}, \code{cluster_size},
#'   and \code{cluster_fraction}.
#'
#' @references
#' Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026).
#' \emph{Reliability Engineering & System Safety}, 275, 112830.
#'
#' @seealso [hm_filter_land_cells()], [hm_select_extreme_events()]
#'
#' @examples
#' \dontrun{
#' X_land <- hm_filter_land_cells(X)
#' cs     <- hm_cluster_size(X_land, threshold = 25,
#'                            time = attr(X_land, "time"))
#' summary(cs$cluster_size)
#' }
#'
#' @export
hm_cluster_size <- function(X, threshold, time = NULL) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold))
    stop("`threshold` must be a single finite numeric value.")

  T_steps <- nrow(X)
  N_cells <- ncol(X)

  if (is.null(time)) {
    time <- seq_len(T_steps)
  } else if (length(time) != T_steps) {
    stop("`time` must have length nrow(X) (", T_steps, ").")
  }

  cluster_size     <- as.integer(rowSums(X >= threshold, na.rm = TRUE))
  cluster_fraction <- cluster_size / N_cells

  data.frame(
    date             = time,
    cluster_size     = cluster_size,
    cluster_fraction = cluster_fraction,
    stringsAsFactors = FALSE
  )
}
