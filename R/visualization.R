if (getRversion() >= "2.15.1") {utils::globalVariables(c("lon", "lat", "group_id", "poly_id", "exceedance_value"))}

#' Validate spatial bounding box
#'
#' @param bbox Optional named numeric vector with `lon_min`, `lon_max`, `lat_min`, and `lat_max`.
#' @param coords Coordinate data.frame with `lon` and `lat`.
#'
#' @return A named numeric vector.
#' @noRd
.hm_validate_bbox <- function(bbox, coords)
                     {if (is.null(coords) || !is.data.frame(coords)) stop("`coords` must be a data.frame.")
                      if (!all(c("lon","lat") %in% names(coords))) stop("`coords` must contain `lon` and `lat` columns.")

                      if (is.null(bbox)) {bbox <- c(lon_min = min(coords$lon, na.rm = TRUE),
                                                                  lon_max = max(coords$lon, na.rm = TRUE),
                                                                  lat_min = min(coords$lat, na.rm = TRUE),
                                                                  lat_max = max(coords$lat, na.rm = TRUE))}

                      if (!is.numeric(bbox) || is.null(names(bbox)) || !all(c("lon_min","lon_max","lat_min","lat_max") %in% names(bbox)))
                      stop("`bbox` must be a named numeric vector with `lon_min`, `lon_max`, `lat_min`, and `lat_max`.")

                      bbox <- bbox[c("lon_min","lon_max","lat_min","lat_max")]

                      if (!all(is.finite(bbox))) stop("All values in `bbox` must be finite.")
                      if (bbox[["lon_min"]] >= bbox[["lon_max"]]) stop("`bbox['lon_min']` must be smaller than `bbox['lon_max']`.")
                      if (bbox[["lat_min"]] >= bbox[["lat_max"]]) stop("`bbox['lat_min']` must be smaller than `bbox['lat_max']`.")

                      bbox}



#' Infer coordinate spacing
#'
#' @param z Numeric coordinate vector.
#'
#' @return Numeric spacing.
#' @noRd
.hm_regular_spacing <- function(z)
                       {z <- sort(unique(as.vector(z)))
                        d <- diff(z)
                        d <- d[is.finite(d) & d > 0]
                        if (length(d) == 0L) stop("Cannot infer grid spacing from coordinate values.")
                            stats::median(d, na.rm = TRUE)}



#' Build regular coordinate sequence
#'
#' @param zmin Minimum coordinate.
#' @param zmax Maximum coordinate.
#' @param dz Coordinate spacing.
#'
#' @return Numeric vector.
#' @noRd
.hm_regular_sequence <- function(zmin, zmax, dz)
                        {if (!is.finite(dz) || dz <= 0) stop("`dz` must be a positive finite number.")
                         z <- seq(zmin, zmax + dz / 1000, by = dz)
                         z[z <= zmax + dz / 1000]}



#' Format longitude labels
#'
#' @param x Numeric longitude values.
#'
#' @return Character vector.
#' @noRd
.hm_label_lon <- function(x) {out <- paste0(abs(x), "\u00B0", ifelse(x < 0, "W", "E"))
                              out[x == 0] <- paste0("0", "\u00B0")
                              out}



#' Format latitude labels
#'
#' @param x Numeric latitude values.
#'
#' @return Character vector.
#' @noRd
.hm_label_lat <- function(x) {out <- paste0(abs(x), "\u00B0", ifelse(x < 0, "S", "N"))
                              out[x == 0] <- paste0("0", "\u00B0")
                              out}



#' Resolve visualization dataset type
#'
#' @param x An `hm_hazard` object.
#' @param dataset Dataset/grid type. One of `"auto"`, `"copernicus"`,
#'   `"eurocordex"`, `"regular"`, or `"rotated"`.
#'
#' @return Either `"regular"` or `"rotated"`.
#' @noRd
.hm_resolve_plot_type <- function(x, dataset = "auto")
                         {dataset <- match.arg(dataset, c("auto","copernicus","eurocordex","regular","rotated"))

                          if (dataset %in% c("copernicus","regular")) return("regular")
                          if (dataset %in% c("eurocordex","rotated")) return("rotated")

                          grid_type <- if (!is.null(x$meta$grid_type)) x$meta$grid_type else NA_character_

                          if (!is.na(grid_type) && grid_type == "regular") return("regular")
                          if (!is.na(grid_type) && grid_type == "rotated") return("rotated")

                          if (!is.null(x$coords) && is.data.frame(x$coords) && all(c("rlon","rlat") %in% names(x$coords))) return("rotated")
                          "regular"}



#' Build regular-grid lines
#'
#' Builds straight longitude and latitude grid lines for regular gridded datasets,
#' such as Copernicus or ERA5-like products.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#' @param bbox Optional named numeric vector with `lon_min`, `lon_max`,
#'   `lat_min`, and `lat_max`.
#' @param grid_res Optional numeric grid spacing in degrees. If `NULL`, spacing is
#'   inferred from the coordinate table.
#'
#' @return A data.frame with grid-line coordinates.
#' @noRd
.hm_build_regular_grid_lines <- function(x, bbox = NULL, grid_res = NULL)
                                {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                                 if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                                 if (!all(c("lon","lat") %in% names(x$coords))) stop("`x$coords` must contain `lon` and `lat` columns.")

                                 coords <- x$coords
                                 bbox <- .hm_validate_bbox(bbox, coords)

                                 if (is.null(grid_res)) {dx <- .hm_regular_spacing(coords$lon)
                                                        dy <- .hm_regular_spacing(coords$lat)}
                                 else {if (!is.numeric(grid_res) || length(grid_res) != 1L || !is.finite(grid_res) || grid_res <= 0) stop("`grid_res` must be a positive numeric value.")
                                       dx <- grid_res
                                       dy <- grid_res}

                                 lon_u <- .hm_regular_sequence(bbox[["lon_min"]], bbox[["lon_max"]], dx)
                                 lat_u <- .hm_regular_sequence(bbox[["lat_min"]], bbox[["lat_max"]], dy)

                                 lines_lon <- do.call(rbind,
                                                      lapply(seq_along(lon_u),
                                                             function(i) {data.frame(line_type = "constant_lon",
                                                                          line_id = i,
                                                                          group_id = paste0("constant_lon_", i),
                                                                          lon = lon_u[i],
                                                                          lat = lat_u,
                                                                          point_order = seq_along(lat_u))}
                                                              ))

                                 lines_lat <- do.call(rbind,
                                                      lapply(seq_along(lat_u),
                                                             function(j) {data.frame(line_type = "constant_lat",
                                                                          line_id = j,
                                                                          group_id = paste0("constant_lat_", j),
                                                                          lon = lon_u,
                                                                          lat = lat_u[j],
                                                                          point_order = seq_along(lon_u))}
                                                             ))

                                 rbind(lines_lon, lines_lat)}



#' Build rotated-grid lines
#'
#' Builds curvilinear grid lines for rotated or curvilinear gridded datasets,
#' such as EURO-CORDEX products.
#'
#' @details
#' The function identifies the rotated-grid rows and columns intersecting the
#' requested geographic bounding box. It then connects cells with constant `rlat`
#' and constant `rlon`. This reproduces the logic normally used when plotting
#' EURO-CORDEX grids from two-dimensional `lon` and `lat` coordinates.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#' @param bbox Optional named numeric vector with `lon_min`, `lon_max`,
#'   `lat_min`, and `lat_max`.
#'
#' @return A data.frame with grid-line coordinates.
#' @noRd
.hm_build_rotated_grid_lines <- function(x, bbox = NULL)
                                {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                                 if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                                 if (!all(c("lon","lat","rlon","rlat") %in% names(x$coords))) stop("Rotated grids require `lon`, `lat`, `rlon`, and `rlat` columns in `x$coords`.")

                                 coords <- x$coords
                                 bbox <- .hm_validate_bbox(bbox, coords)

                                 in_bbox <- coords$lon >= bbox[["lon_min"]] &
                                            coords$lon <= bbox[["lon_max"]] &
                                            coords$lat >= bbox[["lat_min"]] &
                                            coords$lat <= bbox[["lat_max"]]

                                 rlon_keep <- sort(unique(coords$rlon[in_bbox]))
                                 rlat_keep <- sort(unique(coords$rlat[in_bbox]))

                                 if (length(rlon_keep) == 0L || length(rlat_keep) == 0L) stop("No rotated-grid rows or columns intersect `bbox`.")

                                 lines_rlat <- do.call(rbind,
                                                       lapply(seq_along(rlat_keep),
                                                              function(j) {z <- coords[coords$rlat == rlat_keep[j] & coords$rlon %in% rlon_keep, , drop = FALSE]
                                                                           z <- z[order(z$rlon), , drop = FALSE]

                                                                           data.frame(line_type = "constant_rlat",
                                                                                      line_id = j,
                                                                                      group_id = paste0("constant_rlat_", j),
                                                                                      lon = z$lon,
                                                                                      lat = z$lat,
                                                                                      point_order = seq_len(nrow(z)))}
                                                              ))

                                 lines_rlon <- do.call(rbind,
                                                       lapply(seq_along(rlon_keep),
                                                              function(i) {z <- coords[coords$rlon == rlon_keep[i] & coords$rlat %in% rlat_keep, , drop = FALSE]
                                                                           z <- z[order(z$rlat), , drop = FALSE]

                                                                           data.frame(line_type = "constant_rlon",
                                                                                      line_id = i,
                                                                                      group_id = paste0("constant_rlon_", i),
                                                                                      lon = z$lon,
                                                                                      lat = z$lat,
                                                                                      point_order = seq_len(nrow(z)))}
                                                              ))

                                 rbind(lines_rlat, lines_rlon)}



#' Create base spatial-domain plot
#'
#' @param bbox Named bounding box.
#' @param show_map Logical. If `TRUE`, adds world borders using the suggested
#'   package `maps`.
#'
#' @return A `ggplot` object.
#' @noRd
.hm_base_domain_plot <- function(bbox, show_map = TRUE)
                        {if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")

                         p <- ggplot2::ggplot()

                         if (isTRUE(show_map)) {if (requireNamespace("maps", quietly = TRUE)) {p <- p + ggplot2::borders(database = "world",
                                                                                                                         fill = "grey95",
                                                                                                                         colour = "black",
                                                                                                                         linewidth = 0.25)}
                         else {warning("Package `maps` is not available; plotting without map background.", call. = FALSE)}}

                        p}



#' Apply common spatial-domain plot formatting
#'
#' @param p A `ggplot` object.
#' @param bbox Named bounding box.
#' @param title Plot title.
#' @param coord Coordinate display method. Either `"fixed"` or `"quickmap"`.
#'
#' @return A `ggplot` object.
#' @noRd
.hm_format_domain_plot <- function(p, bbox, title = NULL, coord = "fixed")
                          {coord <- match.arg(coord, c("fixed","quickmap"))

                           if (coord == "quickmap") {p <- p + ggplot2::coord_quickmap(xlim = c(bbox[["lon_min"]], bbox[["lon_max"]]),
                                                                                      ylim = c(bbox[["lat_min"]], bbox[["lat_max"]]),
                                                                                      expand = FALSE)}
                           else {p <- p + ggplot2::coord_fixed(xlim = c(bbox[["lon_min"]], bbox[["lon_max"]]),
                                                               ylim = c(bbox[["lat_min"]], bbox[["lat_max"]]),
                                                               expand = FALSE)}

                           p +
                            ggplot2::scale_x_continuous(labels = .hm_label_lon) +
                            ggplot2::scale_y_continuous(labels = .hm_label_lat) +
                            ggplot2::labs(title = title,
                                          x = "Longitude",
                                          y = "Latitude") +
                            ggplot2::theme_minimal(base_size = 12) +
                            ggplot2::theme(panel.background = ggplot2::element_rect(fill = "aliceblue", colour = NA),
                                           panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.25),
                                           panel.grid.minor = ggplot2::element_blank(),
                                           plot.title = ggplot2::element_text(hjust = 0.5))}



#' Plot regular spatial domain
#'
#' Internal plotting routine for regular longitude-latitude grids.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#' @param bbox Optional named bounding box.
#' @param title Optional plot title.
#' @param show_map Logical. If `TRUE`, adds a map background.
#' @param show_points Logical. If `TRUE`, adds grid-cell centres.
#' @param grid_res Optional grid resolution in degrees.
#' @param grid_colour Grid-line colour.
#' @param grid_linewidth Grid-line width.
#' @param point_size Point size.
#' @param coord Coordinate display method.
#'
#' @return A `ggplot` object.
#' @noRd
.hm_plot_regular_domain <- function(x, bbox = NULL, title = NULL, show_map = TRUE,
                                    show_points = FALSE, grid_res = NULL,
                                    grid_colour = "grey55", grid_linewidth = 0.35,
                                    point_size = 0.5, coord = "fixed")
                          {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")

                           if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                           if (!all(c("lon","lat") %in% names(x$coords))) stop("`x$coords` must contain `lon` and `lat` columns.")
                           if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")

                           coords <- x$coords
                           bbox <- .hm_validate_bbox(bbox, coords)
                           grid_lines <- .hm_build_regular_grid_lines(x, bbox = bbox, grid_res = grid_res)

                           if (is.null(title)) title <- "Spatial domain grid (regular grid)"

                                               p <- .hm_base_domain_plot(bbox, show_map = show_map)

                                               p <- p + ggplot2::geom_path(data = grid_lines,
                                                                           ggplot2::aes(x = lon, y = lat, group = group_id),
                                                                           colour = grid_colour,
                                                                           linewidth = grid_linewidth)

                          if (isTRUE(show_points)) {point_coords <- coords[coords$lon >= bbox[["lon_min"]] &
                                                                    coords$lon <= bbox[["lon_max"]] &
                                                                    coords$lat >= bbox[["lat_min"]] &
                                                                    coords$lat <= bbox[["lat_max"]], , drop = FALSE]

                                                   p <- p + ggplot2::geom_point(data = point_coords,
                                                                                ggplot2::aes(x = lon, y = lat),
                                                                                inherit.aes = FALSE,
                                                                                size = point_size)}

                          .hm_format_domain_plot(p, bbox = bbox, title = title, coord = coord)}



#' Plot rotated spatial domain
#'
#' Internal plotting routine for rotated or curvilinear grids.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#' @param bbox Optional named bounding box.
#' @param title Optional plot title.
#' @param show_map Logical. If `TRUE`, adds a map background.
#' @param show_points Logical. If `TRUE`, adds grid-cell centres.
#' @param grid_colour Grid-line colour.
#' @param grid_linewidth Grid-line width.
#' @param point_size Point size.
#' @param coord Coordinate display method.
#'
#' @return A `ggplot` object.
#' @noRd
.hm_plot_rotated_domain <- function(x, bbox = NULL, title = NULL, show_map = TRUE,
                                    show_points = FALSE, grid_colour = "grey55",
                                    grid_linewidth = 0.35, point_size = 0.5,
                                    coord = "fixed")
                          {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                           if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                           if (!all(c("lon","lat","rlon","rlat") %in% names(x$coords))) stop("Rotated grids require `lon`, `lat`, `rlon`, and `rlat` columns in `x$coords`.")
                           if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")

                           coords <- x$coords
                           bbox <- .hm_validate_bbox(bbox, coords)
                           grid_lines <- .hm_build_rotated_grid_lines(x, bbox = bbox)

                           if (is.null(title)) title <- "Spatial domain grid (rotated grid)"

                                               p <- .hm_base_domain_plot(bbox, show_map = show_map)

                                                p <- p + ggplot2::geom_path(data = grid_lines,
                                                                            ggplot2::aes(x = lon, y = lat, group = group_id),
                                                                            colour = grid_colour,
                                                                            linewidth = grid_linewidth)

                           if (isTRUE(show_points)) {point_coords <- coords[coords$lon >= bbox[["lon_min"]] &
                                                     coords$lon <= bbox[["lon_max"]] &
                                                     coords$lat >= bbox[["lat_min"]] &
                                                     coords$lat <= bbox[["lat_max"]], , drop = FALSE]

                                                     p <- p + ggplot2::geom_point(data = point_coords,
                                                                                  ggplot2::aes(x = lon, y = lat),
                                                                                  inherit.aes = FALSE,
                                                                                  size = point_size)}

                          .hm_format_domain_plot(p, bbox = bbox, title = title, coord = coord)}



#' Plot spatial domain grid
#'
#' Plots the spatial domain of a hazard dataset as a grid in geographic
#' longitude-latitude coordinates.
#'
#' @details
#' This function is a public dispatcher. It selects a specialized plotting
#' routine according to the dataset/grid type.
#'
#' \itemize{
#'   \item `dataset = "copernicus"` or `"regular"` uses a regular longitude-latitude
#'   plotting routine.
#'   \item `dataset = "eurocordex"` or `"rotated"` uses a rotated-grid plotting
#'   routine based on `rlon` and `rlat`.
#'   \item `dataset = "auto"` uses `x$meta$grid_type` when available.
#' }
#'
#' Regular grids are plotted as straight longitude-latitude grid lines. Rotated
#' grids are plotted by connecting grid rows and columns in geographic coordinates,
#' which is required for EURO-CORDEX-like datasets where `lat` and `lon` are
#' two-dimensional coordinate matrices.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#' @param dataset Dataset/grid type. One of `"auto"`, `"copernicus"`,
#'   `"eurocordex"`, `"regular"`, or `"rotated"`.
#' @param bbox Optional named numeric vector with `lon_min`, `lon_max`, `lat_min`,
#'   and `lat_max`.
#' @param title Optional plot title.
#' @param show_map Logical. If `TRUE`, a world map is added using the suggested
#'   package `maps`.
#' @param show_points Logical. If `TRUE`, grid-cell centres are added as points.
#' @param grid_res Optional grid resolution in degrees. Used only for regular
#'   grids. If `NULL`, it is inferred from the coordinate table.
#' @param grid_colour Colour used for grid lines.
#' @param grid_linewidth Line width used for grid lines.
#' @param point_size Point size used when `show_points = TRUE`.
#' @param coord Coordinate display method. Either `"fixed"` or `"quickmap"`.
#'
#' @return A `ggplot` object.
#' @export
hm_plot_spatial_domain <- function(x, dataset = "auto", bbox = NULL, title = NULL,
                                   show_map = TRUE, show_points = FALSE,
                                   grid_res = NULL, grid_colour = "grey55",
                                   grid_linewidth = 0.35, point_size = 0.5,
                                   coord = "fixed")
                          {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                           if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                           if (!all(c("lon","lat") %in% names(x$coords))) stop("`x$coords` must contain `lon` and `lat` columns.")

                           plot_type <- .hm_resolve_plot_type(x, dataset = dataset)

                           if (plot_type == "rotated") {.hm_plot_rotated_domain(x = x,
                                                                                bbox = bbox,
                                                                                title = title,
                                                                                show_map = show_map,
                                                                                show_points = show_points,
                                                                                grid_colour = grid_colour,
                                                                                grid_linewidth = grid_linewidth,
                                                                                point_size = point_size,
                                                                                coord = coord)}
                           else {.hm_plot_regular_domain(x = x,
                                                        bbox = bbox,
                                                        title = title,
                                                        show_map = show_map,
                                                        show_points = show_points,
                                                        grid_res = grid_res,
                                                        grid_colour = grid_colour,
                                                        grid_linewidth = grid_linewidth,
                                                        point_size = point_size,
                                                        coord = coord)}}





#' Compute exceedance metrics by grid cell
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#' @param threshold Numeric threshold used to define exceedances.
#'
#' @return A data.frame with coordinates and exceedance metrics.
#' @noRd
.hm_exceedance_metrics <- function(x, threshold)
                          {if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                           if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                           if (!all(c("lon","lat") %in% names(x$coords))) stop("`x$coords` must contain `lon` and `lat` columns.")
                           if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) stop("`threshold` must be one finite numeric value.")

                           if (is.null(x$time) || !inherits(x$time, c("Date","POSIXct","POSIXt"))) {x <- hm_decode_time(x)}

                           X <- hm_to_matrix(x)
                           coords <- attr(X, "coords")

                           if (is.null(coords) || !is.data.frame(coords)) stop("No coordinate metadata found in `hm_to_matrix()` output.")
                           if (ncol(X) != nrow(coords)) stop("Number of matrix columns and coordinate rows are not consistent.")

                           valid_n <- colSums(!is.na(X))
                           count <- colSums(X >= threshold, na.rm = TRUE)

                           coords$exceedance_count <- as.integer(count)
                           coords$exceedance_frequency <- ifelse(valid_n > 0, count / valid_n, NA_real_)

                           coords}



#' Build regular-grid exceedance polygons
#'
#' @param coords Coordinate data.frame with exceedance metrics.
#' @param bbox Named bounding box used to select cell centres.
#' @param grid_res Optional grid resolution in degrees.
#' @param metric Metric to plot.
#'
#' @return A data.frame of cell polygons.
#' @noRd
.hm_exceedance_regular_polygons <- function(coords, bbox, grid_res = NULL, metric = "count")
                                   {if (is.null(grid_res)) {dx <- .hm_regular_spacing(coords$lon)
                                                            dy <- .hm_regular_spacing(coords$lat)}
                                    else {if (!is.numeric(grid_res) || length(grid_res) != 1L || !is.finite(grid_res) || grid_res <= 0) stop("`grid_res` must be a positive numeric value.")
                                          dx <- grid_res
                                          dy <- grid_res}

                                    z <- coords[coords$lon >= bbox[["lon_min"]] &
                                         coords$lon <= bbox[["lon_max"]] &
                                         coords$lat >= bbox[["lat_min"]] &
                                         coords$lat <= bbox[["lat_max"]], , drop = FALSE]

                                    if (nrow(z) == 0L) stop("No grid cells intersect `bbox`.")

                                    value <- if (metric == "frequency") z$exceedance_frequency else z$exceedance_count

                                    do.call(rbind,
                                            lapply(seq_len(nrow(z)),
                                            function(k) {lon0 <- z$lon[k]
                                                         lat0 <- z$lat[k]

                                            data.frame(poly_id = paste0("cell_", k),
                                                       lon = c(lon0 - dx / 2, lon0 + dx / 2,
                                                               lon0 + dx / 2, lon0 - dx / 2,
                                                               lon0 - dx / 2),
                                                       lat = c(lat0 - dy / 2, lat0 - dy / 2,
                                                               lat0 + dy / 2, lat0 + dy / 2,
                                                               lat0 - dy / 2),
                                                       exceedance_value = value[k])}))}



#' Build a corner matrix from cell-centre coordinates
#'
#' @param z Numeric matrix of cell-centre coordinates.
#'
#' @return Numeric matrix of cell-corner coordinates.
#' @noRd
.hm_exceedance_corner_matrix <- function(z)
                                {nr <- nrow(z)
                                 nc <- ncol(z)

                                 if (nr < 2L || nc < 2L) stop("At least two rows and two columns are required to build rotated-grid polygons.")

                                 e <- matrix(NA_real_, nrow = nr + 2L, ncol = nc + 2L)
                                 e[2:(nr + 1L), 2:(nc + 1L)] <- z

                                 e[1L, 2:(nc + 1L)] <- 2 * z[1L, ] - z[2L, ]
                                 e[nr + 2L, 2:(nc + 1L)] <- 2 * z[nr, ] - z[nr - 1L, ]

                                 e[2:(nr + 1L), 1L] <- 2 * z[, 1L] - z[, 2L]
                                 e[2:(nr + 1L), nc + 2L] <- 2 * z[, nc] - z[, nc - 1L]

                                 e[1L, 1L] <- 2 * e[1L, 2L] - e[1L, 3L]
                                 e[1L, nc + 2L] <- 2 * e[1L, nc + 1L] - e[1L, nc]
                                 e[nr + 2L, 1L] <- 2 * e[nr + 2L, 2L] - e[nr + 2L, 3L]
                                 e[nr + 2L, nc + 2L] <- 2 * e[nr + 2L, nc + 1L] - e[nr + 2L, nc]

                                 out <- matrix(NA_real_, nrow = nr + 1L, ncol = nc + 1L)

                                 for (i in seq_len(nr + 1L)) {for (j in seq_len(nc + 1L)) {out[i, j] <- mean(e[i:(i + 1L), j:(j + 1L)], na.rm = TRUE)}}

                                 out}



#' Build rotated-grid exceedance polygons
#'
#' @param coords Coordinate data.frame with rotated-grid coordinates and exceedance metrics.
#' @param bbox Named bounding box used to select cell centres.
#' @param metric Metric to plot.
#'
#' @return A data.frame of cell polygons.
#' @noRd
.hm_exceedance_rotated_polygons <- function(coords, bbox, metric = "count")
                                   {if (!all(c("lon","lat","rlon","rlat","exceedance_count","exceedance_frequency") %in% names(coords))) stop("Rotated grids require `lon`, `lat`, `rlon`, `rlat`, and exceedance metrics.")

                                    rlon_u <- sort(unique(coords$rlon))
                                    rlat_u <- sort(unique(coords$rlat))

                                    nr <- length(rlon_u)
                                    nc <- length(rlat_u)

                                    lon_mat <- matrix(NA_real_, nrow = nr, ncol = nc)
                                    lat_mat <- matrix(NA_real_, nrow = nr, ncol = nc)

                                    ii <- match(coords$rlon, rlon_u)
                                    jj <- match(coords$rlat, rlat_u)

                                    lon_mat[cbind(ii, jj)] <- coords$lon
                                    lat_mat[cbind(ii, jj)] <- coords$lat

                                    lon_corner <- .hm_exceedance_corner_matrix(lon_mat)
                                    lat_corner <- .hm_exceedance_corner_matrix(lat_mat)

                                    keep <- coords$lon >= bbox[["lon_min"]] &
                                            coords$lon <= bbox[["lon_max"]] &
                                            coords$lat >= bbox[["lat_min"]] &
                                            coords$lat <= bbox[["lat_max"]]

                                    if (!any(keep)) stop("No grid cells intersect `bbox`.")

                                    value <- if (metric == "frequency") coords$exceedance_frequency else coords$exceedance_count

                                    cells <- data.frame(i = ii[keep],
                                                        j = jj[keep],
                                                        value = value[keep])

                                    do.call(rbind,
                                            lapply(seq_len(nrow(cells)),
                                    function(k) {i <- cells$i[k]
                                                 j <- cells$j[k]

                                    data.frame(poly_id = paste0("cell_", k),
                                               lon = c(lon_corner[i, j],
                                                       lon_corner[i + 1L, j],
                                                       lon_corner[i + 1L, j + 1L],
                                                       lon_corner[i, j + 1L],
                                                       lon_corner[i, j]),
                                               lat = c(lat_corner[i, j],
                                                       lat_corner[i + 1L, j],
                                                       lat_corner[i + 1L, j + 1L],
                                                       lat_corner[i, j + 1L],
                                                       lat_corner[i, j]),
                                               exceedance_value = cells$value[k])}))}



#' Build bounding box from polygon coordinates
#'
#' @param cells Polygon data.frame with `lon` and `lat`.
#'
#' @return A named numeric bounding box.
#' @noRd
.hm_exceedance_polygon_bbox <- function(cells)
                               {if (is.null(cells) || !is.data.frame(cells)) stop("`cells` must be a data.frame.")
                                if (!all(c("lon","lat") %in% names(cells))) stop("`cells` must contain `lon` and `lat` columns.")

                                c(lon_min = min(cells$lon, na.rm = TRUE),
                                  lon_max = max(cells$lon, na.rm = TRUE),
                                  lat_min = min(cells$lat, na.rm = TRUE),
                                  lat_max = max(cells$lat, na.rm = TRUE))}



#' Merge two bounding boxes
#'
#' @param bbox_a First named bounding box.
#' @param bbox_b Second named bounding box.
#'
#' @return A named numeric bounding box covering both inputs.
#' @noRd
.hm_exceedance_union_bbox <- function(bbox_a, bbox_b)
                             {bbox_a <- bbox_a[c("lon_min","lon_max","lat_min","lat_max")]
                              bbox_b <- bbox_b[c("lon_min","lon_max","lat_min","lat_max")]

                              c(lon_min = min(bbox_a[["lon_min"]], bbox_b[["lon_min"]]),
                                lon_max = max(bbox_a[["lon_max"]], bbox_b[["lon_max"]]),
                                lat_min = min(bbox_a[["lat_min"]], bbox_b[["lat_min"]]),
                                lat_max = max(bbox_a[["lat_max"]], bbox_b[["lat_max"]]))}



#' Plot threshold exceedances
#'
#' Plots threshold exceedance count or frequency at each spatial grid cell.
#'
#' @details
#' Exceedances are defined as time steps for which the hazard value is greater
#' than or equal to `threshold`.
#'
#' By default, the function plots absolute exceedance count, i.e. the number of
#' valid time steps exceeding the threshold at each grid cell. Set
#' `metric = "frequency"` to plot the fraction of valid time steps exceeding the
#' threshold.
#'
#' Filled cells are drawn using cell-boundary polygons. Therefore, the colour
#' layer and the visible grid correspond to the same cell geometry.
#'
#' @param x An `hm_hazard` object after [hm_standardize_coords()].
#' @param threshold Numeric threshold used to define exceedances.
#' @param metric Metric to plot. Either `"count"` or `"frequency"`.
#' @param dataset Dataset/grid type. One of `"auto"`, `"copernicus"`,
#'   `"eurocordex"`, `"regular"`, or `"rotated"`.
#' @param bbox Optional named numeric vector with `lon_min`, `lon_max`, `lat_min`,
#'   and `lat_max`.
#' @param title Optional plot title.
#' @param show_map Logical. If `TRUE`, a world map is added when package `maps`
#'   is available.
#' @param show_grid Logical. If `TRUE`, draws cell borders.
#' @param grid_res Optional grid resolution in degrees. Used only for regular grids.
#' @param colours Colour vector used for the filled cells.
#' @param grid_colour Colour used for cell borders.
#' @param grid_linewidth Line width used for cell borders.
#' @param cell_alpha Transparency of coloured cells.
#' @param legend_title Optional legend title.
#' @param coord Coordinate display method. Either `"fixed"` or `"quickmap"`.
#'
#' @return A `ggplot` object.
#' @export
hm_plot_exceedance <- function(x, threshold, metric = "count",
                               dataset = "auto", bbox = NULL, title = NULL,
                               show_map = TRUE, show_grid = TRUE,
                               grid_res = NULL,
                               colours = c("white", "gold", "orange", "firebrick"),
                               grid_colour = "grey55", grid_linewidth = 0.25,
                               cell_alpha = 0.9, legend_title = NULL,
                               coord = "fixed")
                      {metric <- match.arg(metric, c("count","frequency"))

                       if (is.null(x) || !"hm_hazard" %in% class(x)) stop("`x` must be an object of class `hm_hazard`.")
                       if (is.null(x$coords) || !is.data.frame(x$coords)) stop("`x$coords` must be a data.frame. Run `hm_standardize_coords()` first.")
                       if (!all(c("lon","lat") %in% names(x$coords))) stop("`x$coords` must contain `lon` and `lat` columns.")
                       if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.")

                       coords <- .hm_exceedance_metrics(x, threshold)
                       bbox <- .hm_validate_bbox(bbox, coords)
                       plot_type <- .hm_resolve_plot_type(x, dataset = dataset)

                       cells <- if (plot_type == "regular") {.hm_exceedance_regular_polygons(coords, bbox = bbox, grid_res = grid_res, metric = metric)}
                                else {.hm_exceedance_rotated_polygons(coords, bbox = bbox, metric = metric)}

                       cell_bbox <- .hm_exceedance_polygon_bbox(cells)
                       plot_bbox <- .hm_exceedance_union_bbox(bbox, cell_bbox)

                       if (is.null(title)) {title <- if (metric == "count") paste0("Exceedance count (threshold = ", threshold, ")")
                       else paste0("Exceedance frequency (threshold = ", threshold, ")")}

                       if (is.null(legend_title)) {legend_title <- if (metric == "count") "Exceedance days" else "Frequency"}

                       p <- .hm_base_domain_plot(plot_bbox, show_map = show_map)

                       p <- p + ggplot2::geom_polygon(data = cells,
                                                      ggplot2::aes(x = lon,
                                                                   y = lat,
                                                                   group = poly_id,
                                                                   fill = exceedance_value),
                                                      colour = if (isTRUE(show_grid)) grid_colour else NA,
                                                      linewidth = if (isTRUE(show_grid)) grid_linewidth else 0,
                                                      alpha = cell_alpha)

                      if (isTRUE(show_map) && requireNamespace("maps", quietly = TRUE)) {p <- p + ggplot2::borders(database = "world",
                                                                                                                   fill = NA,
                                                                                                                   colour = "black",
                                                                                                                   linewidth = 0.25)}

                      if (metric == "frequency") {p <- p + ggplot2::scale_fill_gradientn(colours = colours,
                                                                                         name = legend_title,
                                                                                         labels = function(z) paste0(round(100 * z, 2), "%"))}
                      else {p <- p + ggplot2::scale_fill_gradientn(colours = colours,
                                                                   name = legend_title)}

                      .hm_format_domain_plot(p, bbox = plot_bbox, title = title, coord = coord)}
