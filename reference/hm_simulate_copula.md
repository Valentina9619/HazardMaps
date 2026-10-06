# Simulate spatially coherent hazard scenarios from a fitted copula

Draws `n_scenarios` realisations from a fitted copula and
back-transforms to physical units via the inverse marginal CDFs.
Implements Eq. 7–8 of Clavijo Mesa et al. (2026).

## Usage

``` r
hm_simulate_copula(copula, fits, n_scenarios = 1000L, seed = NULL)
```

## Arguments

- copula:

  An `hm_copula` object from \[hm_fit_copula()\].

- fits:

  data.frame from \[hm_fit_marginals()\] for the same land cells used to
  fit `copula`.

- n_scenarios:

  Number of scenarios. Defaults to 1000.

- seed:

  Optional integer random seed.

## Value

A numeric matrix (\\S \times N\\) of simulated hazard intensities.
Attributes `type`, `n_scenarios`, and `n_cells` are attached.

## Details

Vine copulas are sampled via
[`rvinecopulib::rvinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop_methods.html).
The Gaussian copula uses
[`MASS::mvrnorm()`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) followed
by the standard normal CDF. Both packages are listed in `Suggests`.

## References

Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). *Reliability
Engineering & System Safety*, 275, 112830.

## See also

\[hm_fit_copula()\], \[hm_fit_marginals()\]

## Examples

``` r
if (FALSE) { # \dontrun{
cop       <- hm_fit_copula(pit$U, type = "cvine")
scenarios <- hm_simulate_copula(cop, fits, n_scenarios = 1000)
dim(scenarios)
range(scenarios, na.rm = TRUE)
} # }
```
