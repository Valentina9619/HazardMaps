if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "lon", "lat", "group_id", "poly_id", "exceedance_value",
    "period_label", "best_dist_label", "count", "dist_label",
    "bin_center", "mean_tau", "n_pairs", "pct", "density",
    "pdf", "y", "cdf", "theoretical", "empirical"
  ))
}

# @noRd
.hm_validate_bbox <- function(bbox, coords) {

  if (is.null(coords) || !is.data.frame(coords))
    stop("`coords` must be a data.frame.")
  if (!all(c("lon", "lat") %in% names(coords)))
    stop("`coords` must contain `lon` and `lat` columns.")

  if (is.null(bbox)) {
    bbox <- c(
      lon_min = min(coords$lon, na.rm = TRUE),
      lon_max = max(coords$lon, na.rm = TRUE),
      lat_min = min(coords$lat, na.rm = TRUE),
      lat_max = max(coords$lat, na.rm = TRUE)
    )
  }

  required <- c("lon_min", "lon_max", "lat_min", "lat_max")
  if (!is.numeric(bbox) || !all(required %in% names(bbox)))
    stop("`bbox` must be a named numeric vector with `lon_min`, `lon_max`, ",
         "`lat_min`, and `lat_max`.")

  bbox <- bbox[required]

  if (!all(is.finite(bbox)))
    stop("All values in `bbox` must be finite.")
  if (bbox[["lon_min"]] >= bbox[["lon_max"]])
    stop("`bbox['lon_min']` must be smaller than `bbox['lon_max']`.")
  if (bbox[["lat_min"]] >= bbox[["lat_max"]])
    stop("`bbox['lat_min']` must be smaller than `bbox['lat_max']`.")

  bbox
}



# @noRd
.hm_regular_spacing <- function(z) {
  z <- sort(unique(as.vector(z)))
  d <- diff(z)
  d <- d[is.finite(d) & d > 0]
  if (length(d) == 0L)
    stop("Cannot infer grid spacing from coordinate values.")
  stats::median(d, na.rm = TRUE)
}

# @noRd
.hm_regular_sequence <- function(zmin, zmax, dz) {
  if (!is.finite(dz) || dz <= 0)
    stop("`dz` must be a positive finite number.")
  z <- seq(zmin, zmax + dz / 1000, by = dz)
  z[z <= zmax + dz / 1000]
}

# @noRd
.hm_label_lon <- function(x) {
  out        <- paste0(abs(x), "\u00B0", ifelse(x < 0, "W", "E"))
  out[x == 0] <- "0\u00B0"
  out
}

# @noRd
.hm_label_lat <- function(x) {
  out        <- paste0(abs(x), "\u00B0", ifelse(x < 0, "S", "N"))
  out[x == 0] <- "0\u00B0"
  out
}



# @noRd
.hm_resolve_plot_type <- function(x, dataset = "auto") {

  dataset <- match.arg(dataset,
                       c("auto", "copernicus", "eurocordex", "regular", "rotated"))

  if (dataset %in% c("copernicus", "regular")) return("regular")
  if (dataset %in% c("eurocordex", "rotated")) return("rotated")

  grid_type <- x$meta$grid_type %||% NA_character_

  if (!is.na(grid_type) && grid_type == "regular") return("regular")
  if (!is.na(grid_type) && grid_type == "rotated") return("rotated")

  if (!is.null(x$coords) && is.data.frame(x$coords) &&
      all(c("rlon", "rlat") %in% names(x$coords)))
    return("rotated")

  "regular"
}



# @noRd
.hm_build_regular_grid_lines <- function(x, bbox = NULL, grid_res = NULL) {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!all(c("lon", "lat") %in% names(x$coords)))
    stop("`x$coords` must contain `lon` and `lat` columns.")

  coords <- x$coords
  bbox   <- .hm_validate_bbox(bbox, coords)

  if (is.null(grid_res)) {
    dx <- .hm_regular_spacing(coords$lon)
    dy <- .hm_regular_spacing(coords$lat)
  } else {
    if (!is.numeric(grid_res) || length(grid_res) != 1L ||
        !is.finite(grid_res) || grid_res <= 0)
      stop("`grid_res` must be a positive numeric value.")
    dx <- grid_res
    dy <- grid_res
  }

  lon_u <- .hm_regular_sequence(bbox[["lon_min"]], bbox[["lon_max"]], dx)
  lat_u <- .hm_regular_sequence(bbox[["lat_min"]], bbox[["lat_max"]], dy)

  lines_lon <- do.call(rbind, lapply(seq_along(lon_u), function(i) {
    data.frame(line_type   = "constant_lon",
               line_id     = i,
               group_id    = paste0("constant_lon_", i),
               lon         = lon_u[i],
               lat         = lat_u,
               point_order = seq_along(lat_u))
  }))

  lines_lat <- do.call(rbind, lapply(seq_along(lat_u), function(j) {
    data.frame(line_type   = "constant_lat",
               line_id     = j,
               group_id    = paste0("constant_lat_", j),
               lon         = lon_u,
               lat         = lat_u[j],
               point_order = seq_along(lon_u))
  }))

  rbind(lines_lon, lines_lat)
}



# @noRd
.hm_build_rotated_grid_lines <- function(x, bbox = NULL) {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!all(c("lon", "lat", "rlon", "rlat") %in% names(x$coords)))
    stop("Rotated grids require `lon`, `lat`, `rlon`, and `rlat` in `x$coords`.")

  coords <- x$coords
  bbox   <- .hm_validate_bbox(bbox, coords)

  in_bbox <- coords$lon >= bbox[["lon_min"]] & coords$lon <= bbox[["lon_max"]] &
             coords$lat >= bbox[["lat_min"]] & coords$lat <= bbox[["lat_max"]]

  rlon_keep <- sort(unique(coords$rlon[in_bbox]))
  rlat_keep <- sort(unique(coords$rlat[in_bbox]))

  if (length(rlon_keep) == 0L || length(rlat_keep) == 0L)
    stop("No rotated-grid rows or columns intersect `bbox`.")

  lines_rlat <- do.call(rbind, lapply(seq_along(rlat_keep), function(j) {
    z <- coords[coords$rlat == rlat_keep[j] & coords$rlon %in% rlon_keep, ,
                drop = FALSE]
    z <- z[order(z$rlon), , drop = FALSE]
    data.frame(line_type   = "constant_rlat",
               line_id     = j,
               group_id    = paste0("constant_rlat_", j),
               lon         = z$lon,
               lat         = z$lat,
               point_order = seq_len(nrow(z)))
  }))

  lines_rlon <- do.call(rbind, lapply(seq_along(rlon_keep), function(i) {
    z <- coords[coords$rlon == rlon_keep[i] & coords$rlat %in% rlat_keep, ,
                drop = FALSE]
    z <- z[order(z$rlat), , drop = FALSE]
    data.frame(line_type   = "constant_rlon",
               line_id     = i,
               group_id    = paste0("constant_rlon_", i),
               lon         = z$lon,
               lat         = z$lat,
               point_order = seq_len(nrow(z)))
  }))

  rbind(lines_rlat, lines_rlon)
}



# @noRd
.hm_base_domain_plot <- function(bbox, show_map = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")

  p <- ggplot2::ggplot()

  if (isTRUE(show_map)) {
    if (requireNamespace("maps", quietly = TRUE)) {
      p <- p + ggplot2::borders(database = "world",
                                fill      = "grey95",
                                colour    = "black",
                                linewidth = 0.25)
    } else {
      warning("Package `maps` not available; plotting without map background.",
              call. = FALSE)
    }
  }
  p
}

# @noRd
.hm_format_domain_plot <- function(p, bbox, title = NULL, coord = "fixed") {

  coord <- match.arg(coord, c("fixed", "quickmap"))

  xlim <- c(bbox[["lon_min"]], bbox[["lon_max"]])
  ylim <- c(bbox[["lat_min"]], bbox[["lat_max"]])

  if (coord == "quickmap") {
    p <- p + ggplot2::coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE)
  } else {
    p <- p + ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE)
  }

  p +
    ggplot2::scale_x_continuous(labels = .hm_label_lon) +
    ggplot2::scale_y_continuous(labels = .hm_label_lat) +
    ggplot2::labs(title = title, x = "Longitude", y = "Latitude") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "aliceblue", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0.5)
    )
}



# @noRd
.hm_plot_regular_domain <- function(x, bbox = NULL, title = NULL, show_map = TRUE,
                                    show_points = FALSE, grid_res = NULL,
                                    grid_colour = "grey55", grid_linewidth = 0.35,
                                    point_size = 0.5, coord = "fixed") {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!all(c("lon", "lat") %in% names(x$coords)))
    stop("`x$coords` must contain `lon` and `lat` columns.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")

  coords     <- x$coords
  bbox       <- .hm_validate_bbox(bbox, coords)
  grid_lines <- .hm_build_regular_grid_lines(x, bbox = bbox, grid_res = grid_res)

  if (is.null(title))
    title <- "Spatial domain grid (regular grid)"

  p <- .hm_base_domain_plot(bbox, show_map = show_map)
  p <- p + ggplot2::geom_path(
    data      = grid_lines,
    ggplot2::aes(x = lon, y = lat, group = group_id),
    colour    = grid_colour,
    linewidth = grid_linewidth
  )

  if (isTRUE(show_points)) {
    point_coords <- coords[coords$lon >= bbox[["lon_min"]] &
                           coords$lon <= bbox[["lon_max"]] &
                           coords$lat >= bbox[["lat_min"]] &
                           coords$lat <= bbox[["lat_max"]], , drop = FALSE]
    p <- p + ggplot2::geom_point(
      data        = point_coords,
      ggplot2::aes(x = lon, y = lat),
      inherit.aes = FALSE,
      size        = point_size
    )
  }

  .hm_format_domain_plot(p, bbox = bbox, title = title, coord = coord)
}



# @noRd
.hm_plot_rotated_domain <- function(x, bbox = NULL, title = NULL, show_map = TRUE,
                                    show_points = FALSE, grid_colour = "grey55",
                                    grid_linewidth = 0.35, point_size = 0.5,
                                    coord = "fixed") {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!all(c("lon", "lat", "rlon", "rlat") %in% names(x$coords)))
    stop("Rotated grids require `lon`, `lat`, `rlon`, and `rlat` in `x$coords`.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")

  coords     <- x$coords
  bbox       <- .hm_validate_bbox(bbox, coords)
  grid_lines <- .hm_build_rotated_grid_lines(x, bbox = bbox)

  if (is.null(title))
    title <- "Spatial domain grid (rotated grid)"

  p <- .hm_base_domain_plot(bbox, show_map = show_map)
  p <- p + ggplot2::geom_path(
    data      = grid_lines,
    ggplot2::aes(x = lon, y = lat, group = group_id),
    colour    = grid_colour,
    linewidth = grid_linewidth
  )

  if (isTRUE(show_points)) {
    point_coords <- coords[coords$lon >= bbox[["lon_min"]] &
                           coords$lon <= bbox[["lon_max"]] &
                           coords$lat >= bbox[["lat_min"]] &
                           coords$lat <= bbox[["lat_max"]], , drop = FALSE]
    p <- p + ggplot2::geom_point(
      data        = point_coords,
      ggplot2::aes(x = lon, y = lat),
      inherit.aes = FALSE,
      size        = point_size
    )
  }

  .hm_format_domain_plot(p, bbox = bbox, title = title, coord = coord)
}



#' Plot spatial domain grid
#'
#' Plots the spatial grid of a hazard dataset. Detects regular
#' (Copernicus/ERA5) and rotated (EURO-CORDEX) grids automatically or from
#' \code{dataset}.
#'
#' @param x An \code{hm_hazard} object after [hm_standardize_coords()].
#' @param dataset Grid type: \code{"auto"}, \code{"copernicus"},
#'   \code{"eurocordex"}, \code{"regular"}, or \code{"rotated"}.
#' @param bbox Optional named numeric vector with \code{lon_min},
#'   \code{lon_max}, \code{lat_min}, \code{lat_max}.
#' @param title Optional plot title.
#' @param show_map Logical. Add world map background (requires \code{maps}).
#' @param show_points Logical. Add grid-cell centre points.
#' @param grid_res Grid resolution in degrees (regular grids only).
#' @param grid_colour Grid-line colour.
#' @param grid_linewidth Grid-line width.
#' @param point_size Point size when \code{show_points = TRUE}.
#' @param coord \code{"fixed"} or \code{"quickmap"}.
#'
#' @return A \code{ggplot} object.
#' @export
hm_plot_spatial_domain <- function(x, dataset = "auto", bbox = NULL, title = NULL,
                                   show_map = TRUE, show_points = FALSE,
                                   grid_res = NULL, grid_colour = "grey55",
                                   grid_linewidth = 0.35, point_size = 0.5,
                                   coord = "fixed") {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!all(c("lon", "lat") %in% names(x$coords)))
    stop("`x$coords` must contain `lon` and `lat` columns.")

  plot_type <- .hm_resolve_plot_type(x, dataset = dataset)

  if (plot_type == "rotated") {
    .hm_plot_rotated_domain(x = x, bbox = bbox, title = title,
                            show_map = show_map, show_points = show_points,
                            grid_colour = grid_colour,
                            grid_linewidth = grid_linewidth,
                            point_size = point_size, coord = coord)
  } else {
    .hm_plot_regular_domain(x = x, bbox = bbox, title = title,
                            show_map = show_map, show_points = show_points,
                            grid_res = grid_res, grid_colour = grid_colour,
                            grid_linewidth = grid_linewidth,
                            point_size = point_size, coord = coord)
  }
}





# @noRd
.hm_exceedance_metrics <- function(x, threshold) {

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!all(c("lon", "lat") %in% names(x$coords)))
    stop("`x$coords` must contain `lon` and `lat` columns.")
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("`threshold` must be one finite numeric value.")

  if (!inherits(x$time, c("Date", "POSIXct", "POSIXt")))
    x <- hm_decode_time(x)

  X      <- hm_to_matrix(x)
  coords <- attr(X, "coords")

  if (is.null(coords) || !is.data.frame(coords))
    stop("No coordinate metadata found in `hm_to_matrix()` output.")

  valid_n <- colSums(!is.na(X))
  count   <- colSums(X >= threshold, na.rm = TRUE)

  coords$exceedance_count     <- as.integer(count)
  coords$exceedance_frequency <- ifelse(valid_n > 0, count / valid_n, NA_real_)
  coords
}



# @noRd
.hm_exceedance_regular_polygons <- function(coords, bbox, grid_res = NULL,
                                            metric = "count") {

  if (is.null(grid_res)) {
    dx <- .hm_regular_spacing(coords$lon)
    dy <- .hm_regular_spacing(coords$lat)
  } else {
    if (!is.numeric(grid_res) || length(grid_res) != 1L ||
        !is.finite(grid_res) || grid_res <= 0)
      stop("`grid_res` must be a positive numeric value.")
    dx <- grid_res
    dy <- grid_res
  }

  z <- coords[coords$lon >= bbox[["lon_min"]] &
              coords$lon <= bbox[["lon_max"]] &
              coords$lat >= bbox[["lat_min"]] &
              coords$lat <= bbox[["lat_max"]], , drop = FALSE]

  if (nrow(z) == 0L)
    stop("No grid cells intersect `bbox`.")

  value <- if (metric == "frequency") z$exceedance_frequency else z$exceedance_count

  do.call(rbind, lapply(seq_len(nrow(z)), function(k) {
    lon0 <- z$lon[k]
    lat0 <- z$lat[k]
    data.frame(
      poly_id          = paste0("cell_", k),
      lon              = c(lon0 - dx / 2, lon0 + dx / 2,
                           lon0 + dx / 2, lon0 - dx / 2, lon0 - dx / 2),
      lat              = c(lat0 - dy / 2, lat0 - dy / 2,
                           lat0 + dy / 2, lat0 + dy / 2, lat0 - dy / 2),
      exceedance_value = value[k]
    )
  }))
}



# @noRd
.hm_exceedance_corner_matrix <- function(z) {

  nr <- nrow(z)
  nc <- ncol(z)

  if (nr < 2L || nc < 2L)
    stop("At least 2 rows and 2 columns required for rotated-grid polygons.")

  e <- matrix(NA_real_, nrow = nr + 2L, ncol = nc + 2L)
  e[2:(nr + 1L), 2:(nc + 1L)] <- z

  e[1L, 2:(nc + 1L)]      <- 2 * z[1L, ]  - z[2L, ]
  e[nr + 2L, 2:(nc + 1L)] <- 2 * z[nr, ]  - z[nr - 1L, ]
  e[2:(nr + 1L), 1L]      <- 2 * z[, 1L]  - z[, 2L]
  e[2:(nr + 1L), nc + 2L] <- 2 * z[, nc]  - z[, nc - 1L]

  e[1L, 1L]             <- 2 * e[1L, 2L]       - e[1L, 3L]
  e[1L, nc + 2L]        <- 2 * e[1L, nc + 1L]  - e[1L, nc]
  e[nr + 2L, 1L]        <- 2 * e[nr + 2L, 2L]  - e[nr + 2L, 3L]
  e[nr + 2L, nc + 2L]   <- 2 * e[nr + 2L, nc + 1L] - e[nr + 2L, nc]

  out <- matrix(NA_real_, nrow = nr + 1L, ncol = nc + 1L)
  for (i in seq_len(nr + 1L)) {
    for (j in seq_len(nc + 1L)) {
      out[i, j] <- mean(e[i:(i + 1L), j:(j + 1L)], na.rm = TRUE)
    }
  }
  out
}



# @noRd
.hm_exceedance_rotated_polygons <- function(coords, bbox, metric = "count") {

  required <- c("lon", "lat", "rlon", "rlat",
                "exceedance_count", "exceedance_frequency")
  if (!all(required %in% names(coords)))
    stop("Rotated grids require `lon`, `lat`, `rlon`, `rlat`, and exceedance metrics.")

  rlon_u <- sort(unique(coords$rlon))
  rlat_u <- sort(unique(coords$rlat))
  nr     <- length(rlon_u)
  nc     <- length(rlat_u)

  lon_mat <- matrix(NA_real_, nrow = nr, ncol = nc)
  lat_mat <- matrix(NA_real_, nrow = nr, ncol = nc)

  ii <- match(coords$rlon, rlon_u)
  jj <- match(coords$rlat, rlat_u)

  lon_mat[cbind(ii, jj)] <- coords$lon
  lat_mat[cbind(ii, jj)] <- coords$lat

  lon_corner <- .hm_exceedance_corner_matrix(lon_mat)
  lat_corner <- .hm_exceedance_corner_matrix(lat_mat)

  keep <- coords$lon >= bbox[["lon_min"]] & coords$lon <= bbox[["lon_max"]] &
          coords$lat >= bbox[["lat_min"]] & coords$lat <= bbox[["lat_max"]]

  if (!any(keep))
    stop("No grid cells intersect `bbox`.")

  value <- if (metric == "frequency") coords$exceedance_frequency else coords$exceedance_count
  cells <- data.frame(i = ii[keep], j = jj[keep], value = value[keep])

  do.call(rbind, lapply(seq_len(nrow(cells)), function(k) {
    i <- cells$i[k]
    j <- cells$j[k]
    data.frame(
      poly_id          = paste0("cell_", k),
      lon              = c(lon_corner[i, j],     lon_corner[i + 1L, j],
                           lon_corner[i + 1L, j + 1L], lon_corner[i, j + 1L],
                           lon_corner[i, j]),
      lat              = c(lat_corner[i, j],     lat_corner[i + 1L, j],
                           lat_corner[i + 1L, j + 1L], lat_corner[i, j + 1L],
                           lat_corner[i, j]),
      exceedance_value = cells$value[k]
    )
  }))
}



# @noRd
.hm_exceedance_polygon_bbox <- function(cells) {

  if (is.null(cells) || !is.data.frame(cells))
    stop("`cells` must be a data.frame.")
  if (!all(c("lon", "lat") %in% names(cells)))
    stop("`cells` must contain `lon` and `lat` columns.")

  c(lon_min = min(cells$lon, na.rm = TRUE),
    lon_max = max(cells$lon, na.rm = TRUE),
    lat_min = min(cells$lat, na.rm = TRUE),
    lat_max = max(cells$lat, na.rm = TRUE))
}

# @noRd
.hm_exceedance_union_bbox <- function(bbox_a, bbox_b) {

  keys   <- c("lon_min", "lon_max", "lat_min", "lat_max")
  bbox_a <- bbox_a[keys]
  bbox_b <- bbox_b[keys]

  c(lon_min = min(bbox_a[["lon_min"]], bbox_b[["lon_min"]]),
    lon_max = max(bbox_a[["lon_max"]], bbox_b[["lon_max"]]),
    lat_min = min(bbox_a[["lat_min"]], bbox_b[["lat_min"]]),
    lat_max = max(bbox_a[["lat_max"]], bbox_b[["lat_max"]]))
}



#' Plot threshold exceedances
#'
#' Maps exceedance count or frequency at each grid cell. Cells are drawn as
#' filled polygons for both regular and rotated grids.
#'
#' @param x An \code{hm_hazard} object after [hm_standardize_coords()].
#' @param threshold Numeric threshold defining exceedances.
#' @param metric \code{"count"} (default) or \code{"frequency"}.
#' @param dataset Grid type: \code{"auto"}, \code{"copernicus"},
#'   \code{"eurocordex"}, \code{"regular"}, or \code{"rotated"}.
#' @param bbox Optional named numeric vector with \code{lon_min},
#'   \code{lon_max}, \code{lat_min}, \code{lat_max}.
#' @param title Optional plot title.
#' @param show_map Logical. Add world map background.
#' @param show_grid Logical. Draw cell borders.
#' @param grid_res Grid resolution in degrees (regular grids only).
#' @param colours Colour ramp for filled cells.
#' @param grid_colour Cell border colour.
#' @param grid_linewidth Cell border width.
#' @param cell_alpha Cell transparency.
#' @param legend_title Optional legend title.
#' @param coord \code{"fixed"} or \code{"quickmap"}.
#'
#' @return A \code{ggplot} object.
#' @export
hm_plot_exceedance <- function(x, threshold, metric = "count",
                               dataset = "auto", bbox = NULL, title = NULL,
                               show_map = TRUE, show_grid = TRUE,
                               grid_res = NULL,
                               colours = c("white", "gold", "orange", "firebrick"),
                               grid_colour = "grey55", grid_linewidth = 0.25,
                               cell_alpha = 0.9, legend_title = NULL,
                               coord = "fixed") {

  metric <- match.arg(metric, c("count", "frequency"))

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!all(c("lon", "lat") %in% names(x$coords)))
    stop("`x$coords` must contain `lon` and `lat` columns.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")

  coords    <- .hm_exceedance_metrics(x, threshold)
  bbox      <- .hm_validate_bbox(bbox, coords)
  plot_type <- .hm_resolve_plot_type(x, dataset = dataset)

  cells <- if (plot_type == "regular") {
    .hm_exceedance_regular_polygons(coords, bbox = bbox,
                                    grid_res = grid_res, metric = metric)
  } else {
    .hm_exceedance_rotated_polygons(coords, bbox = bbox, metric = metric)
  }

  plot_bbox <- .hm_exceedance_union_bbox(bbox, .hm_exceedance_polygon_bbox(cells))

  if (is.null(title))
    title <- if (metric == "count")
               paste0("Exceedance count (threshold = ", threshold, ")")
             else
               paste0("Exceedance frequency (threshold = ", threshold, ")")

  if (is.null(legend_title))
    legend_title <- if (metric == "count") "Exceedance days" else "Frequency"

  p <- .hm_base_domain_plot(plot_bbox, show_map = show_map)
  p <- p + ggplot2::geom_polygon(
    data      = cells,
    ggplot2::aes(x = lon, y = lat, group = poly_id, fill = exceedance_value),
    colour    = if (isTRUE(show_grid)) grid_colour else NA,
    linewidth = if (isTRUE(show_grid)) grid_linewidth else 0,
    alpha     = cell_alpha
  )

  if (isTRUE(show_map) && requireNamespace("maps", quietly = TRUE))
    p <- p + ggplot2::borders(database = "world", fill = NA,
                               colour = "black", linewidth = 0.25)

  if (metric == "frequency") {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = colours, name = legend_title,
      labels  = function(z) paste0(round(100 * z, 2), "%")
    )
  } else {
    p <- p + ggplot2::scale_fill_gradientn(colours = colours, name = legend_title)
  }

  .hm_format_domain_plot(p, bbox = plot_bbox, title = title, coord = coord)
}



# @noRd
.hm_exceedance_time_summary <- function(x, threshold, by = "year") {

  by <- match.arg(by, c("year", "2_years", "5_years", "decade"))

  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold))
    stop("`threshold` must be one finite numeric value.")

  if (!inherits(x$time, c("Date", "POSIXct", "POSIXt")))
    x <- hm_decode_time(x)

  X    <- hm_to_matrix(x)
  time <- as.Date(x$time)
  year <- as.integer(format(time, "%Y"))

  step <- switch(by,
    year    = 1L,
    `2_years` = 2L,
    `5_years` = 5L,
    decade  = 10L
  )

  period_start <- floor(year / step) * step
  period_end   <- period_start + step - 1L

  event_day <- as.integer(
    rowSums(!is.na(X)) > 0 & rowSums(X >= threshold, na.rm = TRUE) > 0
  )

  out <- stats::aggregate(event_day,
                          by  = list(period_start = period_start,
                                     period_end   = period_end),
                          FUN = sum, na.rm = TRUE)

  names(out)[names(out) == "x"] <- "exceedance_value"
  out <- out[order(out$period_start), , drop = FALSE]

  out$period_label <- if (step == 1L) {
    as.character(out$period_start)
  } else {
    paste0(out$period_start, "-", out$period_end)
  }

  out$period_label <- factor(out$period_label, levels = out$period_label)
  out
}



#' Plot extreme-event days over time
#'
#' Bar chart of extreme-event days aggregated by year, two-year, five-year
#' or decade periods. A day is counted once when at least one cell exceeds
#' \code{threshold}.
#'
#' @param x An \code{hm_hazard} object after [hm_standardize_coords()].
#' @param threshold Numeric threshold defining exceedances.
#' @param by Aggregation period: \code{"year"}, \code{"2_years"},
#'   \code{"5_years"}, or \code{"decade"}.
#' @param title Optional plot title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param bar_fill Bar fill colour.
#' @param bar_colour Bar border colour.
#'
#' @return A \code{ggplot} object.
#' @export
hm_plot_exceedance_time <- function(x, threshold, by = "year",
                                    title = NULL, xlab = NULL, ylab = NULL,
                                    bar_fill = "grey70", bar_colour = "grey30") {

  by <- match.arg(by, c("year", "2_years", "5_years", "decade"))

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")

  dat <- .hm_exceedance_time_summary(x = x, threshold = threshold, by = by)

  if (is.null(title))
    title <- paste0("Extreme-event days (threshold = ", threshold, ")")

  if (is.null(xlab))
    xlab <- switch(by,
      year      = "Year",
      `2_years` = "Two-year period",
      `5_years` = "Five-year period",
      decade    = "Decade"
    )

  if (is.null(ylab))
    ylab <- "Extreme-event days"

  ggplot2::ggplot(dat, ggplot2::aes(x = period_label, y = exceedance_value)) +
    ggplot2::geom_col(fill = bar_fill, colour = bar_colour, linewidth = 0.25) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5),
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1)
    )
}



#' Plot best-fitting marginal distribution by grid cell
#'
#' Maps each cell by the marginal family selected by [hm_fit_marginals()].
#' Reproduces Fig. 6 of Clavijo Mesa et al. (2026).
#'
#' @param x An \code{hm_hazard} object after [hm_standardize_coords()].
#' @param fits The data.frame returned by [hm_fit_marginals()].
#' @param bbox Optional named numeric vector with `lon_min`, `lon_max`,
#'   `lat_min`, and `lat_max` to control the plot extent. If `NULL`
#'   (default), the extent is derived from the coordinate table with a
#'   small padding so edge points are not clipped.
#' @param title Optional plot title.
#' @param show_map Logical. If `TRUE`, adds a world map background.
#' @param palette Optional named character vector mapping distribution labels
#'   to colours. Names must match the values in `fits$best_dist_label`
#'   (e.g. `"Weibull"`, `"Gamma"`, `"Lognormal"`, `"Gumbel"`,
#'   `"Zero-trunc. Gaussian"`, `"Zero-trunc. Laplace"`, `"No fit"`).
#'   If `NULL` (default), a built-in muted palette is used.
#' @param coord Coordinate display method. Either `"fixed"` or `"quickmap"`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [hm_fit_marginals()], [hm_plot_marginal_fit()],
#'   [hm_plot_dist_frequency()]
#'
#' @examples
#' \dontrun{
#' x     <- hm_read_netcdf(f, var = "i10fg")
#' x     <- hm_standardize_coords(x)
#' x     <- hm_decode_time(x)
#' x_ext <- hm_select_extreme_events(x, threshold = 25)
#' X     <- hm_to_matrix(x_ext)
#' fits  <- hm_fit_marginals(X)
#'
#' hm_plot_best_dist(x_ext, fits)
#' }
#'
#' @export
hm_plot_best_dist <- function(x, fits, bbox = NULL, title = NULL,
                               show_map = TRUE, palette = NULL,
                               coord = "fixed") {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")
  if (is.null(x) || !"hm_hazard" %in% class(x))
    stop("`x` must be an object of class `hm_hazard`.")
  if (is.null(x$coords) || !is.data.frame(x$coords))
    stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame returned by `hm_fit_marginals()`.")
  if (nrow(fits) != nrow(x$coords))
    stop("`fits` and `x$coords` must have the same number of rows.")

  coords <- x$coords
  coords$best_dist_label <- fits$best_dist_label

  coords$best_dist_label[is.na(coords$best_dist_label)] <- "No fit"

  if (is.null(bbox)) {
    pad_lon <- .hm_regular_spacing(coords$lon) / 2
    pad_lat <- .hm_regular_spacing(coords$lat) / 2
    bbox <- c(
      lon_min = min(coords$lon, na.rm = TRUE) - pad_lon,
      lon_max = max(coords$lon, na.rm = TRUE) + pad_lon,
      lat_min = min(coords$lat, na.rm = TRUE) - pad_lat,
      lat_max = max(coords$lat, na.rm = TRUE) + pad_lat
    )
  } else {
    bbox <- .hm_validate_bbox(bbox, coords)
  }

  # Colour palette — muted tones that are distinguishable without being
  # visually overwhelming. Users can override via the `palette` argument.
  default_colours <- c(
    "Weibull"              = "#7BAFD4",   # soft steel blue
    "Gamma"                = "#8FC98F",   # soft green
    "Lognormal"            = "#F4A95A",   # soft amber
    "Gumbel"               = "#C8A0C8",   # soft lavender
    "Zero-trunc. Gaussian" = "#E88080",   # soft rose
    "Zero-trunc. Laplace"  = "#C4956A",   # soft tan
    "No fit"               = "grey85"
  )

  families <- sort(unique(coords$best_dist_label))

  if (is.null(palette)) {
    palette <- default_colours[families]
    palette[is.na(palette)] <- "grey85"
  } else {
    # User-supplied palette: fill any missing families with grey
    missing_fam <- setdiff(families, names(palette))
    palette[missing_fam] <- "grey85"
    palette <- palette[families]
  }

  if (is.null(title))
    title <- "Best-fitting marginal distribution by grid cell"

  # Build filled cell polygons reusing the same helpers as hm_plot_exceedance().
  # This ensures regular and rotated grids are handled identically — regular
  # grids get axis-aligned rectangles, rotated grids get true curvilinear
  # quadrilaterals via the corner-matrix interpolation already implemented.
  #
  # The only change vs hm_plot_exceedance() is that cells are coloured by
  # distribution family (discrete) instead of exceedance value (continuous).

  plot_type <- .hm_resolve_plot_type(x)

  # Temporarily attach a numeric dummy column so we can reuse the polygon
  # builders, then swap in the distribution label afterwards
  coords$exceedance_count     <- seq_len(nrow(coords))
  coords$exceedance_frequency <- seq_len(nrow(coords))

  if (plot_type == "regular") {
    cells_raw <- .hm_exceedance_regular_polygons(
      coords   = coords,
      bbox     = bbox,
      grid_res = NULL,
      metric   = "count"
    )
  } else {
    cells_raw <- .hm_exceedance_rotated_polygons(
      coords = coords,
      bbox   = bbox,
      metric = "count"
    )
  }

  # Map poly_id back to distribution label.
  # poly_id is "cell_k" where k is the row index within the bbox subset.
  # We recover k and look up the label from coords.
  poly_ids <- unique(cells_raw$poly_id)
  k_index  <- as.integer(sub("cell_", "", poly_ids))

  # For regular grids the subset is bbox-filtered; rebuild the mapping
  if (plot_type == "regular") {
    coords_in_bbox <- coords[
      coords$lon >= bbox[["lon_min"]] & coords$lon <= bbox[["lon_max"]] &
      coords$lat >= bbox[["lat_min"]] & coords$lat <= bbox[["lat_max"]], ,
      drop = FALSE
    ]
  } else {
    coords_in_bbox <- coords[
      coords$lon >= bbox[["lon_min"]] & coords$lon <= bbox[["lon_max"]] &
      coords$lat >= bbox[["lat_min"]] & coords$lat <= bbox[["lat_max"]], ,
      drop = FALSE
    ]
  }

  label_map <- data.frame(
    poly_id        = poly_ids,
    best_dist_label = coords_in_bbox$best_dist_label[k_index],
    stringsAsFactors = FALSE
  )

  cells <- merge(cells_raw[, c("poly_id", "lon", "lat")],
                 label_map, by = "poly_id", sort = FALSE)

  p <- .hm_base_domain_plot(bbox, show_map = show_map)

  p <- p +
    ggplot2::geom_polygon(
      data = cells,
      ggplot2::aes(x = lon, y = lat,
                   group = poly_id,
                   fill  = best_dist_label),
      colour    = "grey80",
      linewidth = 0.1
    ) +
    ggplot2::scale_fill_manual(
      values   = palette,
      name     = "Distribution",
      na.value = "grey85"
    )

  # Draw map borders on top of the filled cells
  if (isTRUE(show_map) && requireNamespace("maps", quietly = TRUE)) {
    p <- p + ggplot2::borders(database = "world",
                               fill      = NA,
                               colour    = "black",
                               linewidth = 0.25)
  }

  .hm_format_domain_plot(p, bbox = bbox, title = title, coord = coord)
}



#' Plot marginal fit diagnostics for a single grid cell
#'
#' Produces a three-panel diagnostic figure for the best-fitting marginal
#' distribution at a chosen grid cell: (1) a histogram of the observed values
#' with the fitted PDF overlaid, (2) the empirical vs fitted CDF, and (3) a
#' quantile-quantile plot. Reproduces the style of Fig. 7 and Appendix B
#' (Fig. B1) in Clavijo Mesa et al. (2026).
#'
#' @param X A numeric matrix (same one passed to [hm_fit_marginals()]).
#' @param fits The data.frame returned by [hm_fit_marginals()].
#' @param cell_id Integer. The cell to diagnose (column index of \code{X}).
#' @param title Optional overall plot title. Defaults to the cell coordinates
#'   when `x` is supplied, or to the cell index otherwise.
#' @param x Optional `hm_hazard` object after [hm_standardize_coords()].
#'   When supplied, the cell's geographic coordinates are shown in the title.
#' @param n_bins Number of histogram bins.
#' @param colour_fit Colour used for the fitted PDF and CDF lines.
#'
#' @return A `ggplot` object (three panels arranged with
#'   [ggplot2::facet_wrap()] on a long-format data frame).
#'
#' @seealso [hm_fit_marginals()], [hm_marginal_gof()],
#'   [hm_plot_best_dist()], [hm_plot_dist_frequency()]
#'
#' @examples
#' \dontrun{
#' fits <- hm_fit_marginals(X)
#' hm_plot_marginal_fit(X, fits, cell_id = 1, x = x_ext)
#' }
#'
#' @export
hm_plot_marginal_fit <- function(X, fits, cell_id, title = NULL, x = NULL,
                                  n_bins = 30, colour_fit = "firebrick") {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")
  if (!is.matrix(X) || !is.numeric(X))
    stop("`X` must be a numeric matrix.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame returned by `hm_fit_marginals()`.")
  if (!is.numeric(cell_id) || length(cell_id) != 1L ||
      cell_id < 1L || cell_id > ncol(X))
    stop("`cell_id` must be a single integer between 1 and ncol(X).")

  cell_id <- as.integer(cell_id)

  # Retrieve GoF data (KS, QQ, params) via the existing function
  gof <- hm_marginal_gof(X, fits, cell_id = cell_id)

  xj    <- X[, cell_id]
  xj    <- xj[is.finite(xj) & xj > 0]
  dist  <- gof$dist
  label <- gof$dist_label

  # Build title from coordinates if hm_hazard object is supplied
  if (is.null(title)) {
    if (!is.null(x) && !is.null(x$coords) && is.data.frame(x$coords) &&
        cell_id <= nrow(x$coords)) {
      lat_c <- round(x$coords$lat[cell_id], 2)
      lon_c <- round(x$coords$lon[cell_id], 2)
      title <- paste0("[", lat_c, "\u00B0N, ", lon_c, "\u00B0E]  \u2014  ",
                      label)
    } else {
      title <- paste0(label, "  \u2014  cell ", cell_id)
    }
  }

  params <- gof$params

  # Fitted PDF and CDF functions
  pfun <- switch(dist,
    weibull  = function(q) stats::pweibull(q, params$shape, params$scale),
    gamma    = function(q) stats::pgamma(q, params$shape, rate = params$rate),
    lnorm    = function(q) stats::plnorm(q, params$meanlog, params$sdlog),
    gumbel   = function(q) .pgumbel(q, params$mu, params$beta),
    tnorm    = function(q) .ptnorm(q, params$mean, params$sd),
    tlaplace = function(q) .ptlaplace(q, params$mu, params$b)
  )

  dfun <- switch(dist,
    weibull  = function(q) stats::dweibull(q, params$shape, params$scale),
    gamma    = function(q) stats::dgamma(q, params$shape, rate = params$rate),
    lnorm    = function(q) stats::dlnorm(q, params$meanlog, params$sdlog),
    gumbel   = function(q) .dgumbel(q, params$mu, params$beta),
    tnorm    = function(q) .dtnorm(q, params$mean, params$sd),
    tlaplace = function(q) .dtlaplace(q, params$mu, params$b)
  )

  # Grid for smooth PDF and CDF curves
  x_grid <- seq(min(xj) * 0.95, max(xj) * 1.05, length.out = 300)
  curve_df <- data.frame(
    x   = x_grid,
    pdf = dfun(x_grid),
    cdf = pfun(x_grid)
  )

  # Empirical CDF
  n     <- length(xj)
  ecdf_df <- data.frame(
    x = sort(xj),
    y = seq_len(n) / n
  )

  # --- Panel 1: histogram + fitted PDF ------------------------------------
  p1 <- ggplot2::ggplot(data.frame(x = xj), ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins    = n_bins,
                            fill    = "grey80",
                            colour  = "white",
                            linewidth = 0.25) +
    ggplot2::geom_line(data = curve_df,
                       ggplot2::aes(x = x, y = pdf),
                       colour    = colour_fit,
                       linewidth = 0.8) +
    ggplot2::labs(x = "Wind gust [m/s]", y = "Density",
                  title = "PDF") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # --- Panel 2: empirical vs fitted CDF ----------------------------------
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_step(data = ecdf_df,
                       ggplot2::aes(x = x, y = y),
                       colour    = "grey40",
                       linewidth = 0.6) +
    ggplot2::geom_line(data = curve_df,
                       ggplot2::aes(x = x, y = cdf),
                       colour    = colour_fit,
                       linewidth = 0.8) +
    ggplot2::labs(x = "Wind gust [m/s]", y = "F(x)",
                  title = "CDF") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # --- Panel 3: QQ plot --------------------------------------------------
  p3 <- ggplot2::ggplot(gof$qq,
                         ggplot2::aes(x = theoretical, y = empirical)) +
    ggplot2::geom_point(size = 0.8, colour = "grey40") +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = colour_fit, linewidth = 0.8) +
    ggplot2::labs(x = "Theoretical quantiles",
                  y = "Empirical quantiles",
                  title = "QQ Plot") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # Combine with patchwork if available, otherwise return a list
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(p1, p2, p3, nrow = 1) +
      patchwork::plot_annotation(title = title,
                                 theme = ggplot2::theme(
                                   plot.title = ggplot2::element_text(
                                     hjust = 0.5, size = 12,
                                     face = "bold")))
    return(combined)
  }

  # Fallback: return list of three ggplot objects
  message("Install the `patchwork` package to display all three panels ",
          "in one figure. Returning a list of three ggplot objects instead.")
  list(pdf = p1, cdf = p2, qq = p3)
}



#' Plot frequency of best-fitting marginal distribution families
#'
#' Draws a bar chart showing how many grid cells were assigned each
#' distribution family by [hm_fit_marginals()]. Reproduces the style of
#' Fig. A1 in Clavijo Mesa et al. (2026, Appendix A).
#'
#' @param fits The data.frame returned by [hm_fit_marginals()].
#' @param title Optional plot title.
#' @param bar_fill Bar fill colour. If `NULL` (default), each bar is coloured
#'   by distribution family using the same palette as [hm_plot_best_dist()].
#' @param bar_colour Bar border colour.
#' @param show_pct Logical. If `TRUE`, adds percentage labels above each bar.
#'
#' @return A `ggplot` object.
#'
#' @seealso [hm_fit_marginals()], [hm_plot_best_dist()],
#'   [hm_plot_marginal_fit()]
#'
#' @examples
#' \dontrun{
#' fits <- hm_fit_marginals(X)
#' hm_plot_dist_frequency(fits)
#' }
#'
#' @export
hm_plot_dist_frequency <- function(fits, title = NULL, bar_fill = NULL,
                                    bar_colour = "white", show_pct = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")
  if (!is.data.frame(fits))
    stop("`fits` must be the data.frame returned by `hm_fit_marginals()`.")

  if (is.null(title))
    title <- "Frequency of best-fitting marginal distributions"

  # Build summary table ordered by frequency
  labels <- fits$best_dist_label
  labels[is.na(labels)] <- "No fit"

  counts   <- sort(table(labels), decreasing = TRUE)
  n_total  <- sum(counts)

  dat <- data.frame(
    dist_label = factor(names(counts), levels = names(counts)),
    count      = as.integer(counts),
    pct        = round(100 * as.integer(counts) / n_total, 1),
    stringsAsFactors = FALSE
  )

  # Same muted palette as hm_plot_best_dist() for visual consistency
  base_colours <- c(
    "Weibull"              = "#7BAFD4",
    "Gamma"                = "#8FC98F",
    "Lognormal"            = "#F4A95A",
    "Gumbel"               = "#C8A0C8",
    "Zero-trunc. Gaussian" = "#E88080",
    "Zero-trunc. Laplace"  = "#C4956A",
    "No fit"               = "grey85"
  )

  p <- ggplot2::ggplot(dat,
                        ggplot2::aes(x = dist_label, y = count,
                                     fill = dist_label)) +
    ggplot2::geom_col(colour = bar_colour, linewidth = 0.25)

  if (is.null(bar_fill)) {
    palette <- base_colours[levels(dat$dist_label)]
    palette[is.na(palette)] <- "grey70"
    p <- p + ggplot2::scale_fill_manual(values = palette, guide = "none")
  } else {
    p <- p + ggplot2::scale_fill_manual(
      values = rep(bar_fill, nlevels(dat$dist_label)),
      guide  = "none"
    )
  }

  if (isTRUE(show_pct)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = paste0(pct, "%")),
      vjust   = -0.4,
      size    = 3.5,
      colour  = "grey30"
    )
  }

  p +
    ggplot2::labs(title = title,
                  x     = "Distribution",
                  y     = "Number of grid cells") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
    )
}



#' Plot the empirical spatial correlogram
#'
#' Draws the empirical correlogram produced by [hm_empirical_correlogram()]:
#' average Kendall \eqn{\hat{\tau}} as a function of inter-site distance,
#' with point size proportional to the number of cell pairs in each bin.
#' Reproduces the style of Fig. 8 in Clavijo Mesa et al. (2026, SF paper).
#'
#' @details
#' A horizontal dashed reference line is drawn at \eqn{\tau = 0} to help
#' identify the distance at which spatial dependence becomes negligible.
#' Point size is scaled by `sqrt(n_pairs)` so bins with more pairs are
#' visually more prominent without dominating the plot.
#'
#' @param corg A data.frame returned by [hm_empirical_correlogram()], with
#'   columns `bin_center`, `mean_tau`, and `n_pairs`.
#' @param title Optional plot title.
#' @param colour Line and point colour. Defaults to `"#2166AC"` (steel blue).
#' @param show_points Logical. If `TRUE` (default), draws points at each
#'   bin centre sized by the number of pairs.
#' @param show_zero Logical. If `TRUE` (default), adds a horizontal dashed
#'   line at \eqn{\tau = 0}.
#' @param ylim Optional numeric vector of length 2 for the y-axis range.
#'   If `NULL` (default), the range is set automatically.
#'
#' @return A `ggplot` object.
#'
#' @seealso [hm_empirical_correlogram()], [hm_empirical_kendall()],
#'   [hm_distance_matrix()]
#'
#' @references
#' Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
#' Inoperability assessment of interdependent critical infrastructures
#' exposed to natural hazards considering climate change.
#' \emph{International Journal of Disaster Risk Reduction}, 141, 106172.
#'
#' @examples
#' \dontrun{
#' X_land <- hm_filter_land_cells(X)
#' fits   <- hm_fit_marginals(X_land)
#' pit    <- hm_pit_transform(X_land, fits)
#' tau    <- hm_empirical_kendall(pit$Z)
#' D      <- hm_distance_matrix(attr(X_land, "coords"))
#' corg   <- hm_empirical_correlogram(tau, D)
#'
#' hm_plot_correlogram(corg,
#'   title = "Empirical correlogram — Copernicus wind gust")
#' }
#'
#' @export
hm_plot_correlogram <- function(corg, title = NULL,
                                 colour = "#2166AC",
                                 show_points = TRUE,
                                 show_zero   = TRUE,
                                 ylim = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package `ggplot2` is required.")
  if (!is.data.frame(corg))
    stop("`corg` must be a data.frame from `hm_empirical_correlogram()`.")
  if (!all(c("bin_center", "mean_tau", "n_pairs") %in% names(corg)))
    stop("`corg` must contain `bin_center`, `mean_tau`, and `n_pairs`.")

  if (is.null(title))
    title <- "Empirical spatial correlogram"

  p <- ggplot2::ggplot(corg,
                        ggplot2::aes(x = bin_center, y = mean_tau))

  # Zero reference line drawn first so it sits behind everything else
  if (isTRUE(show_zero))
    p <- p + ggplot2::geom_hline(yintercept = 0,
                                  linetype  = "dashed",
                                  colour    = "grey60",
                                  linewidth = 0.5)

  # Connecting line
  p <- p + ggplot2::geom_line(colour = colour, linewidth = 0.8)

  # Points sized by number of pairs in each bin
  if (isTRUE(show_points))
    p <- p + ggplot2::geom_point(
      ggplot2::aes(size = sqrt(n_pairs)),
      colour = colour,
      show.legend = FALSE
    )

  p <- p +
    ggplot2::labs(
      title = title,
      x     = "Distance (km)",
      y     = "Kendall \u03C4"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  if (!is.null(ylim)) {
    if (!is.numeric(ylim) || length(ylim) != 2L)
      stop("`ylim` must be a numeric vector of length 2.")
    p <- p + ggplot2::ylim(ylim[1], ylim[2])
  }

  p
}
