# Plot threshold exceedances

Maps exceedance count or frequency at each grid cell. Cells are drawn as
filled polygons for both regular and rotated grids.

## Usage

``` r
hm_plot_exceedance(
  x,
  threshold,
  metric = "count",
  dataset = "auto",
  bbox = NULL,
  title = NULL,
  show_map = TRUE,
  show_grid = TRUE,
  grid_res = NULL,
  colours = c("white", "gold", "orange", "firebrick"),
  grid_colour = "grey55",
  grid_linewidth = 0.25,
  cell_alpha = 0.9,
  legend_title = NULL,
  coord = "fixed"
)
```

## Arguments

- x:

  An `hm_hazard` object after \[hm_standardize_coords()\].

- threshold:

  Numeric threshold defining exceedances.

- metric:

  `"count"` (default) or `"frequency"`.

- dataset:

  Grid type: `"auto"`, `"copernicus"`, `"eurocordex"`, `"regular"`, or
  `"rotated"`.

- bbox:

  Optional named numeric vector with `lon_min`, `lon_max`, `lat_min`,
  `lat_max`.

- title:

  Optional plot title.

- show_map:

  Logical. Add world map background.

- show_grid:

  Logical. Draw cell borders.

- grid_res:

  Grid resolution in degrees (regular grids only).

- colours:

  Colour ramp for filled cells.

- grid_colour:

  Cell border colour.

- grid_linewidth:

  Cell border width.

- cell_alpha:

  Cell transparency.

- legend_title:

  Optional legend title.

- coord:

  `"fixed"` or `"quickmap"`.

## Value

A `ggplot` object.
