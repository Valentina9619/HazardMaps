# Plot one simulated hazard scenario on the spatial grid

Maps a single row of a scenario matrix produced by
\[hm_simulate_copula()\] or \[hm_simulate_sf_kle()\] onto the spatial
grid. Useful for visually inspecting the spatial coherence of generated
pseudoscenarios.

## Usage

``` r
hm_plot_scenario(
  scenarios,
  coords,
  scenario_id = 1L,
  title = NULL,
  show_map = TRUE,
  legend_title = "m/s",
  colours = c("white", "#7BAFD4", "#2166AC", "#1A1A2E"),
  coord = "fixed"
)
```

## Arguments

- scenarios:

  A numeric matrix (\\S \times N\\) returned by \[hm_simulate_copula()\]
  or \[hm_simulate_sf_kle()\].

- coords:

  A data.frame with `lat`, `lon`, and optionally `rlon`/`rlat` columns.
  Typically `attr(X_land, "coords")` from \[hm_filter_land_cells()\].
  Must have `N` rows.

- scenario_id:

  Integer. Row index of the scenario to plot. Defaults to `1L`.

- title:

  Optional plot title.

- show_map:

  Logical. Add world map background (requires `maps`).

- legend_title:

  Optional legend title. Defaults to `"m/s"`.

- colours:

  Colour ramp for cell values.

- coord:

  `"fixed"` or `"quickmap"`.

## Value

A `ggplot` object.

## See also

\[hm_simulate_copula()\], \[hm_simulate_sf_kle()\],
\[hm_plot_exceedance()\]

## Examples

``` r
if (FALSE) { # \dontrun{
X_land    <- hm_filter_land_cells(X)
fits      <- hm_fit_marginals(X_land)
pit       <- hm_pit_transform(X_land, fits)
cop       <- hm_fit_copula(pit$U, type = "cvine")
scenarios <- hm_simulate_copula(cop, fits, n_scenarios = 100)

hm_plot_scenario(scenarios, attr(X_land, "coords"), scenario_id = 1L)
} # }
```
