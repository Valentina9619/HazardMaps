# Package index

## Data preparation

Read, standardize, and reshape gridded NetCDF datasets.

- [`hm_read_netcdf()`](https://valentina9619.github.io/HazardMaps/reference/hm_read_netcdf.md)
  : Read a hazard NetCDF dataset
- [`hm_standardize_coords()`](https://valentina9619.github.io/HazardMaps/reference/hm_standardize_coords.md)
  : Standardize spatial coordinates
- [`hm_get_grid_geometry()`](https://valentina9619.github.io/HazardMaps/reference/hm_get_grid_geometry.md)
  : Extract grid geometry
- [`hm_decode_time()`](https://valentina9619.github.io/HazardMaps/reference/hm_decode_time.md)
  : Decode the NetCDF time axis
- [`hm_to_matrix()`](https://valentina9619.github.io/HazardMaps/reference/hm_to_matrix.md)
  : Convert hazard array to a time-by-cell matrix
- [`hm_filter_land_cells()`](https://valentina9619.github.io/HazardMaps/reference/hm_filter_land_cells.md)
  : Remove sea and no-data cells from a hazard matrix

## Event selection

Filter extreme-event days and summarise their properties.

- [`hm_select_extreme_events()`](https://valentina9619.github.io/HazardMaps/reference/hm_select_extreme_events.md)
  : Select extreme-event days by threshold
- [`hm_summarise_events()`](https://valentina9619.github.io/HazardMaps/reference/hm_summarise_events.md)
  : Summarise retained extreme-event days

## Marginal distributions

Fit, diagnose, and transform marginal distributions at each grid cell.

- [`hm_fit_marginals()`](https://valentina9619.github.io/HazardMaps/reference/hm_fit_marginals.md)
  : Fit marginal distributions to each grid cell
- [`hm_marginal_gof()`](https://valentina9619.github.io/HazardMaps/reference/hm_marginal_gof.md)
  : Goodness-of-fit diagnostics for a single grid cell
- [`hm_pit_transform()`](https://valentina9619.github.io/HazardMaps/reference/hm_pit_transform.md)
  : Apply the Probability Integral Transform

## Spatial dependence

Estimate empirical spatial correlation and cluster-size statistics.

- [`hm_distance_matrix()`](https://valentina9619.github.io/HazardMaps/reference/hm_distance_matrix.md)
  : Compute pairwise haversine distances between grid cells
- [`hm_empirical_kendall()`](https://valentina9619.github.io/HazardMaps/reference/hm_empirical_kendall.md)
  : Compute the empirical pairwise Kendall correlation matrix
- [`hm_empirical_correlogram()`](https://valentina9619.github.io/HazardMaps/reference/hm_empirical_correlogram.md)
  : Compute the empirical spatial correlogram
- [`hm_cluster_size()`](https://valentina9619.github.io/HazardMaps/reference/hm_cluster_size.md)
  : Compute spatial cluster sizes of extreme events

## Copula models

Fit vine and Gaussian copulas and simulate pseudoscenarios.

- [`hm_fit_copula()`](https://valentina9619.github.io/HazardMaps/reference/hm_fit_copula.md)
  : Fit a spatial copula model to the uniform PIT data
- [`hm_simulate_copula()`](https://valentina9619.github.io/HazardMaps/reference/hm_simulate_copula.md)
  : Simulate spatially coherent hazard scenarios from a fitted copula

## Stochastic field models

Fit parametric covariance models and simulate via KLE.

- [`hm_fit_correlation_model()`](https://valentina9619.github.io/HazardMaps/reference/hm_fit_correlation_model.md)
  : Fit a parametric spatial correlation model to the empirical
  correlogram
- [`hm_simulate_sf_kle()`](https://valentina9619.github.io/HazardMaps/reference/hm_simulate_sf_kle.md)
  : Simulate hazard scenarios via Karhunen-Loève Expansion

## Visualization

Plot spatial domains, exceedances, marginals, correlograms, and
scenarios.

- [`hm_plot_spatial_domain()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_spatial_domain.md)
  : Plot spatial domain grid
- [`hm_plot_exceedance()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_exceedance.md)
  : Plot threshold exceedances
- [`hm_plot_exceedance_time()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_exceedance_time.md)
  : Plot extreme-event days over time
- [`hm_plot_best_dist()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_best_dist.md)
  : Plot best-fitting marginal distribution by grid cell
- [`hm_plot_marginal_fit()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_marginal_fit.md)
  : Plot marginal fit diagnostics for a single grid cell
- [`hm_plot_dist_frequency()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_dist_frequency.md)
  : Plot frequency of best-fitting marginal distribution families
- [`hm_plot_correlogram()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_correlogram.md)
  : Plot the empirical spatial correlogram
- [`hm_plot_scenario()`](https://valentina9619.github.io/HazardMaps/reference/hm_plot_scenario.md)
  : Plot one simulated hazard scenario on the spatial grid
