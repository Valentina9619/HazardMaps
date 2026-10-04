#' Compute pairwise haversine distances between grid cells
#'
#' Returns a symmetric distance matrix giving the great-circle distance in
#' kilometres between every pair of grid cells, using the haversine formula.
#' This is the distance metric used throughout the spatial correlation
#' analysis in Clavijo Mesa et al. (2026, SF paper, Eq. 5).
#'
#' @param coords A data.frame with at least `lat` and `lon` columns in
#'   decimal degrees, as stored in `x$coords` after
#'   [hm_standardize_coords()]. Typically the `coords` attribute of the
#'   matrix returned by [hm_filter_land_cells()].
#'
#' @return A symmetric numeric matrix of dimension \eqn{N \times N} where
#'   entry \eqn{(l, o)} is the haversine distance in kilometres between
#'   grid cells \eqn{l} and \eqn{o}. The diagonal is zero.
#'
#' @references
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' Inoperability assessment of interdependent critical infrastructures
#' exposed to natural hazards considering climate change.
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @seealso [hm_empirical_correlogram()], [hm_filter_land_cells()]
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
#' D <- hm_distance_matrix(x$coords)
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
  rE <- 6371  # mean Earth radius in km

  # Convert to radians
  phi <- lat * pi / 180
  psi <- lon * pi / 180

  # Vectorised haversine: compute all N x N pairs at once
  phi_i <- matrix(phi, nrow = N, ncol = N, byrow = FALSE)
  phi_j <- matrix(phi, nrow = N, ncol = N, byrow = TRUE)
  psi_i <- matrix(psi, nrow = N, ncol = N, byrow = FALSE)
  psi_j <- matrix(psi, nrow = N, ncol = N, byrow = TRUE)

  a <- sin((phi_j - phi_i) / 2)^2 +
       cos(phi_i) * cos(phi_j) * sin((psi_j - psi_i) / 2)^2

  # Clamp before asin to avoid NaN from floating-point values > 1.
  # Use matrix() to guarantee dimensions are preserved after arithmetic.
  D <- matrix(2 * rE * asin(pmin(1, sqrt(a))), nrow = N, ncol = N)

  # Ensure diagonal is exactly zero
  D[cbind(seq_len(N), seq_len(N))] <- 0

  # Attach cell labels when available
  if (!is.null(coords$cell_id))
    rownames(D) <- colnames(D) <- as.character(coords$cell_id)

  D
}



#' Compute the empirical pairwise Kendall correlation matrix
#'
#' Estimates Kendall's \eqn{\hat{\tau}} between every pair of grid cells
#' from the Gaussian scores matrix produced by [hm_pit_transform()].
#' This is the spatial dependence measure used in Clavijo Mesa et al.
#' (2026, Copernicus paper, Appendix C, Eq. 20) to validate copula and
#' stochastic-field models.
#'
#' @details
#' Kendall's \eqn{\tau} is computed on the Gaussian scores `Z` rather than
#' on the raw or uniform values, following the practice in Clavijo Mesa et
#' al. (2026). Because \eqn{\tau} is invariant under strictly monotone
#' transformations, the result is identical to computing it on `U` or on
#' the original wind-gust values, but using `Z` avoids numerical issues at
#' the boundaries of the uniform scale.
#'
#' Computing all \eqn{N(N-1)/2} pairs is expensive for large grids. For
#' EURO-CORDEX-scale datasets (N > 1000) consider supplying a random
#' subsample of column indices via `cell_idx`.
#'
#' @param Z A numeric matrix of dimension \eqn{T \times N} of centred
#'   Gaussian scores, as returned in `pit$Z` by [hm_pit_transform()].
#' @param cell_idx Optional integer vector of column indices to use. When
#'   `NULL` (default), all columns are used.
#'
#' @return A symmetric numeric matrix of dimension \eqn{N \times N} (or
#'   \eqn{|}\code{cell_idx}\eqn{| \times |}\code{cell_idx}\eqn{|}) of
#'   pairwise Kendall \eqn{\hat{\tau}} values in \eqn{[-1, 1]}.
#'
#' @references
#' Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). A modeling framework
#' for the inoperability assessment of interdependent critical
#' infrastructures exposed to spatially distributed natural hazards.
#' \emph{Reliability Engineering & System Safety}, 275, 112830.
#'
#' @seealso [hm_pit_transform()], [hm_empirical_correlogram()],
#'   [hm_distance_matrix()]
#'
#' @examples
#' \dontrun{
#' fits <- hm_fit_marginals(X)
#' pit  <- hm_pit_transform(X, fits)
#'
#' # Use a subsample of 30 cells to keep computation fast
#' tau <- hm_empirical_kendall(pit$Z, cell_idx = 1:30)
#' range(tau)
#' }
#'
#' @export
hm_empirical_kendall <- function(Z, cell_idx = NULL) {

  if (!is.matrix(Z) || !is.numeric(Z))
    stop("`Z` must be a numeric matrix of Gaussian scores ",
         "from `hm_pit_transform()`.")

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

      # Drop rows where either column is NA
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
#' Bins pairwise Kendall correlations by inter-site distance to produce an
#' empirical correlogram: average correlation as a function of distance.
#' This is the key input for fitting a parametric spatial covariance model
#' (e.g. Matérn) in the stochastic field module, following Eq. 5-6 in
#' Clavijo Mesa et al. (2026, SF paper).
#'
#' @details
#' Distance bins are constructed by dividing the range
#' \eqn{[0, \texttt{max\_dist}]} into \code{n\_bins} equal-width intervals.
#' Each bin reports the mean Kendall \eqn{\hat{\tau}} over all cell pairs
#' whose inter-site distance falls within that interval, weighted by the
#' number of pairs in the bin (proportional weighting, as in Eq. 6 of
#' Clavijo Mesa et al., 2026).
#'
#' Only the upper triangle of the distance and correlation matrices is used
#' to avoid double-counting.
#'
#' @param tau A symmetric numeric matrix of pairwise Kendall correlations,
#'   as returned by [hm_empirical_kendall()].
#' @param D A symmetric numeric distance matrix in kilometres, as returned
#'   by [hm_distance_matrix()]. Must have the same dimensions as `tau`.
#' @param n_bins Integer. Number of distance bins. Defaults to 30.
#' @param max_dist Maximum distance in kilometres to include. Pairs beyond
#'   this distance are ignored. If `NULL` (default), the maximum observed
#'   distance is used.
#'
#' @return A data.frame with one row per non-empty bin and columns:
#'   \describe{
#'     \item{`bin_center`}{Mid-point of the distance bin in km.}
#'     \item{`mean_tau`}{Average Kendall \eqn{\hat{\tau}} in the bin.}
#'     \item{`n_pairs`}{Number of cell pairs contributing to the bin.}
#'     \item{`bin_lo`}{Lower edge of the bin in km.}
#'     \item{`bin_hi`}{Upper edge of the bin in km.}
#'   }
#'
#' @references
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' Inoperability assessment of interdependent critical infrastructures
#' exposed to natural hazards considering climate change.
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @seealso [hm_empirical_kendall()], [hm_distance_matrix()]
#'
#' @examples
#' \dontrun{
#' tau  <- hm_empirical_kendall(pit$Z, cell_idx = 1:50)
#' D    <- hm_distance_matrix(coords[1:50, ])
#' corg <- hm_empirical_correlogram(tau, D, n_bins = 20)
#' head(corg)
#'
#' # Quick plot
#' plot(corg$bin_center, corg$mean_tau, type = "b",
#'      xlab = "Distance (km)", ylab = "Kendall tau",
#'      main = "Empirical correlogram")
#' abline(h = 0, lty = 2)
#' }
#'
#' @export
hm_empirical_correlogram <- function(tau, D, n_bins = 30L,
                                      max_dist = NULL) {

  if (!is.matrix(tau) || !is.numeric(tau))
    stop("`tau` must be a numeric matrix from `hm_empirical_kendall()`.")
  if (!is.matrix(D) || !is.numeric(D))
    stop("`D` must be a numeric matrix from `hm_distance_matrix()`.")
  if (!identical(dim(tau), dim(D)))
    stop("`tau` and `D` must have the same dimensions.")
  if (!is.numeric(n_bins) || length(n_bins) != 1L || n_bins < 2L)
    stop("`n_bins` must be an integer >= 2.")

  n_bins <- as.integer(n_bins)
  N      <- nrow(tau)

  # Extract upper triangle (excluding diagonal)
  upper   <- upper.tri(D)
  d_vals  <- D[upper]
  t_vals  <- tau[upper]

  # Drop NA pairs
  ok     <- is.finite(d_vals) & is.finite(t_vals)
  d_vals <- d_vals[ok]
  t_vals <- t_vals[ok]

  if (is.null(max_dist))
    max_dist <- max(d_vals)

  in_range <- d_vals <= max_dist
  d_vals   <- d_vals[in_range]
  t_vals   <- t_vals[in_range]

  if (length(d_vals) == 0L)
    stop("No valid pairs found within `max_dist` = ", max_dist, " km.")

  # Build equal-width bins
  breaks     <- seq(0, max_dist, length.out = n_bins + 1L)
  bin_lo     <- breaks[-length(breaks)]
  bin_hi     <- breaks[-1L]
  bin_center <- (bin_lo + bin_hi) / 2

  mean_tau <- numeric(n_bins)
  n_pairs  <- integer(n_bins)

  for (b in seq_len(n_bins)) {
    in_bin       <- d_vals > bin_lo[b] & d_vals <= bin_hi[b]
    n_pairs[b]   <- sum(in_bin)
    mean_tau[b]  <- if (n_pairs[b] > 0L)
                      mean(t_vals[in_bin], na.rm = TRUE)
                    else
                      NA_real_
  }

  # Keep only non-empty bins
  keep <- n_pairs > 0L

  data.frame(
    bin_center = bin_center[keep],
    mean_tau   = mean_tau[keep],
    n_pairs    = n_pairs[keep],
    bin_lo     = bin_lo[keep],
    bin_hi     = bin_hi[keep],
    stringsAsFactors = FALSE
  )
}



#' Compute spatial cluster sizes of extreme events
#'
#' For each extreme-event day, counts the number of grid cells where the
#' hazard value meets or exceeds a given threshold. This is the cluster-size
#' metric \eqn{B(t)} defined in Appendix C of Clavijo Mesa et al. (2026,
#' Copernicus paper, Eq. 21), used to validate whether a spatial dependence
#' model reproduces the observed spatial footprint of extreme events.
#'
#' @details
#' The cluster size on day \eqn{t} is:
#' \deqn{B(t) = \sum_{d_l \in \mathcal{D}} \mathbf{1}\{H^*(d_l, t) \geq \tau\}}
#'
#' A large \eqn{B(t)} indicates a spatially widespread event; a small value
#' indicates a localised event. The distribution of \eqn{B(t)} over all
#' retained extreme-event days characterises the spatial tail behaviour of
#' the hazard field and is used to compare observed vs model-generated
#' scenarios (Table C3 in Clavijo Mesa et al., 2026).
#'
#' @param X A numeric matrix of dimension \eqn{T \times N} (time steps by
#'   land cells), as returned by [hm_filter_land_cells()] applied after
#'   [hm_to_matrix()].
#' @param threshold A single finite numeric value. Cells with values
#'   strictly below this threshold are not counted.
#' @param time Optional vector of dates (length \eqn{T}) to attach to the
#'   output. When `NULL`, row indices are used. Typically
#'   `attr(X, "time")`.
#'
#' @return A data.frame with one row per time step and columns:
#'   \describe{
#'     \item{`date`}{Date or row index of the time step.}
#'     \item{`cluster_size`}{Number of cells \eqn{\geq} threshold on that
#'       day.}
#'     \item{`cluster_fraction`}{`cluster_size` as a fraction of the total
#'       number of land cells.}
#'   }
#'
#' @references
#' Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). A modeling framework
#' for the inoperability assessment of interdependent critical
#' infrastructures exposed to spatially distributed natural hazards.
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
#' hist(cs$cluster_size, main = "Spatial cluster size distribution",
#'      xlab = "Number of cells >= threshold")
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

  if (!is.null(time)) {
    if (length(time) != T_steps)
      stop("`time` must have length equal to `nrow(X)` (", T_steps, ").")
  } else {
    time <- seq_len(T_steps)
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
