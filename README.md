# HazardMaps

<img src="man/figures/logo.png" align="right" height="160" alt="HazardMaps logo" />

HazardMaps is an R package for reproducible spatial dependence analysis of
natural hazards. It provides a unified workflow to ingest gridded climate
datasets, fit and compare spatial dependence models, and generate spatially
coherent pseudoscenarios for uncertainty quantification.

The package implements the methodology of Clavijo Mesa et al. (2026),
validated on wind-gust extremes over Northeast Italy using Copernicus and
EURO-CORDEX reanalysis data.

## Installation

```r
remotes::install_github("Valentina9619/HazardMaps")
```

## Workflow

```r
library(HazardMaps)

# 1. Ingest
x <- hm_read_netcdf("Copernicus_dailymax.nc", var = "i10fg")
x <- hm_standardize_coords(x)
x <- hm_decode_time(x)

# 2. Filter extreme events
x_ext <- hm_select_extreme_events(x, threshold = 25)
X     <- hm_to_matrix(x_ext)
X     <- hm_filter_land_cells(X)

# 3. Fit marginal distributions
fits <- hm_fit_marginals(X, criterion = "bic")
pit  <- hm_pit_transform(X, fits)

# 4. Estimate spatial dependence
D    <- hm_distance_matrix(attr(X, "coords"))
tau  <- hm_empirical_kendall(pit$Z)
corg <- hm_empirical_correlogram(tau, D)

# 5a. Copula model (C-vine, D-vine, or Gaussian)
cop       <- hm_fit_copula(pit$U, type = "cvine")
scenarios <- hm_simulate_copula(cop, fits, n_scenarios = 1000)

# 5b. Stochastic field model (KLE-based)
model     <- hm_fit_correlation_model(corg)
scenarios <- hm_simulate_sf_kle(model, D, fits, n_scenarios = 1000)
```

## Key functions

| Module | Functions |
|---|---|
| Data ingest | `hm_read_netcdf`, `hm_standardize_coords`, `hm_decode_time`, `hm_to_matrix`, `hm_filter_land_cells` |
| Event selection | `hm_select_extreme_events`, `hm_summarise_events` |
| Marginals | `hm_fit_marginals`, `hm_marginal_gof`, `hm_pit_transform` |
| Dependence | `hm_distance_matrix`, `hm_empirical_kendall`, `hm_empirical_correlogram`, `hm_cluster_size` |
| Copulas | `hm_fit_copula`, `hm_simulate_copula` |
| Stochastic fields | `hm_fit_correlation_model`, `hm_simulate_sf_kle` |
| Visualization | `hm_plot_spatial_domain`, `hm_plot_exceedance`, `hm_plot_best_dist`, `hm_plot_correlogram` |

## References

Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). A modeling framework for
the inoperability assessment of interdependent critical infrastructures
exposed to spatially distributed natural hazards. *Reliability Engineering &
System Safety*, 275, 112830.

Clavijo Mesa, M.V., Broggi, M., Di Maio, F. & Zio, E. (2026).
Inoperability assessment of interdependent critical infrastructures exposed
to natural hazards considering climate change. *International Journal of
Disaster Risk Reduction*, 141, 106172.
