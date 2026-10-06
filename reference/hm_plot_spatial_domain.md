# Plot spatial domain grid

Plots the spatial grid of a hazard dataset. Detects regular
(Copernicus/ERA5) and rotated (EURO-CORDEX) grids automatically or from
`dataset`.

## Usage

``` r
hm_plot_spatial_domain(
  x,
  dataset = "auto",
  bbox = NULL,
  title = NULL,
  show_map = TRUE,
  show_points = FALSE,
  grid_res = NULL,
  grid_colour = "grey55",
  grid_linewidth = 0.35,
  point_size = 0.5,
  coord = "fixed"
)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\].

- dataset:

  Grid type: `"auto"`, `"copernicus"`, `"eurocordex"`, `"regular"`, or
  `"rotated"`.

- bbox:

  Optional named numeric vector with `lon_min`, `lon_max`, `lat_min`,
  `lat_max`.

- title:

  Optional plot title.

- show_map:

  Logical. Add world map background (requires `maps`).

- show_points:

  Logical. Add grid-cell centre points.

- grid_res:

  Grid resolution in degrees (regular grids only).

- grid_colour:

  Grid-line colour.

- grid_linewidth:

  Grid-line width.

- point_size:

  Point size when `show_points = TRUE`.

- coord:

  `"fixed"` or `"quickmap"`.

## Value

A `ggplot` object.
