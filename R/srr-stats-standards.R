#' srr_stats
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} Primary references: Clavijo Mesa et al. (2026a, 2026b) are cited in the DESCRIPTION and README.
#' @srrstats {G1.1} HazardMaps is the first R implementation of this spatial dependence pipeline combining marginal fitting, vine copulas, and KLE-based stochastic fields.
#' @srrstats {G1.2} A Life Cycle Statement is included in README describing the current maturing stage, available functionality, and anticipated future development.
#' @srrstats {G1.3} All statistical terminology (vine copula, KLE, Matern covariance, PIT, empirical Kendall correlation) is defined in function documentation and vignettes.
#' @srrstats {G1.4} All exported functions are documented with roxygen2.
#' @srrstats {G1.4a} All internal functions have @noRd tags.
#' @srrstats {G2.0} All functions assert input lengths (e.g. threshold must be length 1, cell_id must be length 1).
#' @srrstats {G2.0a} Input length expectations are documented in @param descriptions for all functions.
#' @srrstats {G2.1} All functions assert input types via is.numeric(), is.matrix(), inherits() checks.
#' @srrstats {G2.1a} Data type expectations are documented in @param descriptions.
#' @srrstats {G2.2} Univariate parameters (threshold, criterion, seed, n_scenarios) are checked for length == 1.
#' @srrstats {G2.3} Character inputs use match.arg() (criterion in hm_fit_marginals, by in hm_exceedance_time_summary).
#' @srrstats {G2.3a} match.arg() is used for criterion and by parameters.
#' @srrstats {G2.8} All exported functions validate inputs and dispatch to internal helpers that receive single defined types.
#' @srrstats {G2.13} Missing data checks are implemented in hm_filter_land_cells, hm_fit_marginals, and hm_to_matrix.
#' @srrstats {G2.15} NA values are explicitly handled with is.finite() and na.rm = TRUE where applicable.
#' @srrstats {G3.1} Covariance calculations use custom Matern, exponential, Gaussian and spherical kernel functions. Users select the kernel family via the fam parameter in hm_fit_correlation_model.
#' @srrstats {G3.1a} The choice of covariance family is documented in hm_fit_correlation_model and demonstrated in the vignette.
#' @srrstats {G5.0} Tests use two real-world datasets included as package data (Copernicus ERA5 and EURO-CORDEX).
#' @srrstats {G5.1} Test datasets are included as extdata and documented.
#' @srrstats {G5.2} Error and warning behaviour is tested for all major input validation checks.
#' @srrstats {G5.2a} All stop() messages are unique across the package.
#' @srrstats {G5.2b} Tests explicitly trigger each error condition and verify the message.
#' @srrstats {G5.3} Return objects are tested for absence of NA values where applicable.
#' @srrstats {G5.4} Correctness tests verify outputs against known properties of the input datasets.
#' @srrstats {G5.4a} As a novel pipeline, correctness is verified against simple analytical cases and known distributional properties.
#' @srrstats {G5.5} All stochastic tests use set.seed() for reproducibility.
#' @srrstats {G5.6} Parameter recovery tests verify that hm_fit_marginals recovers physically plausible Weibull parameters from both Copernicus and EURO-CORDEX datasets.
#' @srrstats {G5.6a} Recovery is verified using tolerance-based range checks on shape and scale parameters.
#' @srrstats {G5.6b} Convergence is verified via conv_code == 0 across all fitted cells in both datasets.
#' @srrstats {G5.8} Edge condition tests are implemented for zero-length data, wrong types, and out-of-range inputs.
#' @srrstats {G5.8a} Zero-length data is tested.
#' @srrstats {G5.8b} Wrong data types are tested.
#' @srrstats {G5.8c} All-NA columns are handled in hm_fit_marginals and tested.
#' @srrstats {G5.9} Noise susceptibility tests verify stability of hm_fit_marginals under trivial noise in both datasets.
#' @srrstats {G5.9a} Tests confirm that adding noise at .Machine$double.eps * 1e6 scale does not change best_dist selection.
#' @srrstats {G5.9b} Tests confirm that hm_simulate_copula returns identical results given the same seed.
#' @srrstats {SP1.0} The coordinate reference system is documented and validated in hm_standardize_coords.
#' @srrstats {SP1.1} The package explicitly documents that inputs must be regular or rotated longitude/latitude grids from NetCDF files.
#' @srrstats {SP1.2} Spatial data classes used (hm_hazard) are documented with their coordinate structure.
#' @srrstats {SP1.3} Coordinate transformations are documented in hm_standardize_coords.
#' @srrstats {SP2.0} Spatial input validation is implemented in hm_standardize_coords and hm_get_grid_geometry.
#' @srrstats {SP2.3} Grid geometry (regular vs rotated) is detected and documented in hm_get_grid_geometry.
#' @srrstats {SP3.0} Spatial operations (distance matrix, correlogram) are documented with their assumptions.
#' @srrstats {SP4.0} Visualisation functions (hm_plot_spatial_domain, hm_plot_exceedance) correctly handle both regular and rotated grids.
#' @srrstats {SP4.2} All spatial plots include coordinate axis labels and titles.
#' @srrstats {PD1.0} References for all distributions are provided in hm_fit_marginals documentation and the associated publications.
#' @srrstats {PD3.0} All distribution manipulations are analytic (closed-form PDF, CDF, quantile functions).
#' @srrstats {PD3.1} Distribution operations are contained in .hm_pfun and .hm_qfun which accept distribution names as parameters.
#' @srrstats {PD3.2} Optimisation in hm_fit_marginals uses stats::optim with Nelder-Mead; all parameters and defaults are documented.
#' @srrstats {PD3.3} hm_fit_marginals return object includes conv_code (optim convergence code) and n_iter (number of function evaluations) for each fitted cell.
#' @srrstats {PD4.0} Numeric outputs of distribution functions are tested against known values.
#' @srrstats {PD4.1} Tests compare distribution function outputs with explicitly derived reference values.
#' @srrstats {PD4.2} hm_fit_marginals is tested with multiple distributions (Weibull, Gamma, Gumbel, etc.).
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G1.5} The package implements a novel pipeline; no prior implementation exists to reproduce against. Vignettes demonstrate the full workflow on included datasets.
#' @srrstatsNA {G1.6} No alternative R implementation of the full pipeline exists. Partial comparisons with rvinecopulib and gstat are noted in the README.
#' @srrstatsNA {G2.3b} No case-sensitive character parameters exist; match.arg() handles case normalisation implicitly.
#' @srrstatsNA {G2.4} The package does not perform type conversions between data types as part of its workflow.
#' @srrstatsNA {G2.4a} No integer conversion is performed on user inputs.
#' @srrstatsNA {G2.4b} No explicit numeric conversion is performed on user inputs.
#' @srrstatsNA {G2.4c} No character conversion is performed on user inputs.
#' @srrstatsNA {G2.4d} No factor conversion is performed.
#' @srrstatsNA {G2.4e} No factor inputs are accepted.
#' @srrstatsNA {G2.5} No factor inputs are accepted by any function.
#' @srrstatsNA {G2.6} One-dimensional inputs are validated for length == 1 before use.
#' @srrstatsNA {G2.7} The package uses a custom hm_hazard class for spatial data, not generic tabular formats, as this encapsulates the spatial metadata required for the workflow.
#' @srrstatsNA {G2.9} No type conversion with information loss is performed on user inputs.
#' @srrstatsNA {G2.10} No single-column extraction from tabular inputs is performed in user-facing functions.
#' @srrstatsNA {G2.11} The package does not accept generic data.frame inputs with non-standard columns.
#' @srrstatsNA {G2.12} List columns in the fits data.frame (params column) are handled explicitly and documented.
#' @srrstatsNA {G2.14} Missing data handling is fixed by design: NA values are excluded from fitting, appropriate for gridded climate data where NA indicates land/sea mask.
#' @srrstatsNA {G2.14a} Errors on missing data are not appropriate; NA cells are skipped with a count returned.
#' @srrstatsNA {G2.14b} NA cells are silently skipped; n_failed is returned to inform the user.
#' @srrstatsNA {G2.14c} Imputation is not appropriate for gridded climate extremes data.
#' @srrstatsNA {G2.16} Inf and NaN values are excluded via is.finite() in hm_fit_marginals and hm_to_matrix.
#' @srrstatsNA {G3.0} No floating point equality comparisons are made; all comparisons use >= or <=.
#' @srrstatsNA {G4.0} The package does not write outputs to local files.
#' @srrstatsNA {G5.4b} No prior R implementation exists for comparison.
#' @srrstatsNA {G5.4c} No published numerical outputs are available for direct comparison.
#' @srrstatsNA {G5.7} Algorithm performance tests are not applicable as HazardMaps implements a novel pipeline with no baseline implementation to compare convergence times against.
#' @srrstatsNA {G5.8d} The package operates on matrices where rows are observations and columns are grid cells; no square-matrix constraint exists.
#' @srrstatsNA {G5.10} Extended tests are not applicable at bronze badge level. All tests run within standard testthat framework using included datasets.
#' @srrstatsNA {G5.11} No additional data downloads are required for testing; all test data is included in the package.
#' @srrstatsNA {G5.11a} Not applicable; no downloads are required for tests.
#' @srrstatsNA {G5.12} Not applicable; no extended tests require special platform requirements or memory constraints.
#' @srrstatsNA {SP1.4} The package does not perform distance calculations in projected coordinates; distances are computed in km from lon/lat using hm_distance_matrix.
#' @srrstatsNA {SP2.1} The package does not use sf or sp objects; it uses a custom hm_hazard class optimised for gridded NetCDF data.
#' @srrstatsNA {SP2.2} No transformation between spatial classes is performed.
#' @srrstatsNA {SP3.1} The package does not perform polygon operations.
#' @srrstatsNA {SP3.2} The package does not perform point-in-polygon operations.
#' @srrstatsNA {SP4.1} The package does not produce interactive maps.
#' @srrstatsNA {SP4.3} The package does not use external tile services.
#' @srrstatsNA {PD2.0} Distributions are implemented directly or via stats/fitdistrplus. A general representation package is not used as these distributions require custom truncation.
#' @srrstatsNA {PD3.4} No numerical integration is used; all CDF/quantile operations are analytic.
#' @srrstatsNA {PD3.5} No discrete summation is used to approximate integrals.
#' @srrstatsNA {PD3.5a} Not applicable; no discrete summation is used.
#' @srrstatsNA {PD4.3} Only Nelder-Mead via stats::optim is appropriate for these distributions; alternative parameter comparisons are covered by the existing distribution-specific tests.
#' @srrstatsNA {PD4.4} No alternative optimisation algorithm is applicable for the custom truncated distributions implemented in HazardMaps.
#' @noRd
NULL
