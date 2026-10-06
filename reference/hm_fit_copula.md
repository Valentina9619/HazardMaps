# Fit a spatial copula model to the uniform PIT data

Fits a C-vine, D-vine, or Gaussian copula to the uniform PIT matrix from
\[hm_pit_transform()\], following Eq. 6 of Clavijo Mesa et al. (2026).
The three structures correspond to the comparison in Table IV of that
paper.

## Usage

``` r
hm_fit_copula(U, type = "cvine", ...)
```

## Arguments

- U:

  Numeric matrix (\\T \times N\\) of uniform PIT values from
  \[hm_pit_transform()\].

- type:

  Copula structure: `"cvine"` (default), `"dvine"`, or `"gaussian"`.

- ...:

  Additional arguments passed to
  [`rvinecopulib::vinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.html).

## Value

An `hm_copula` list with `type`, `model`, `n_cells`, and `n_obs`.

## Details

Vine copulas are fitted via
[`rvinecopulib::vinecop()`](https://vinecopulib.github.io/rvinecopulib/reference/vinecop.html).
The Gaussian copula uses the sin-transform \\\rho = \sin(\pi \hat{\tau}
/ 2)\\ to convert Kendall \\\tau\\ to a Pearson correlation matrix.
Requires `rvinecopulib` (listed in `Suggests`).

## References

Clavijo Mesa, M.V., Di Maio, F. & Zio, E. (2026). *Reliability
Engineering & System Safety*, 275, 112830.

## See also

\[hm_pit_transform()\], \[hm_simulate_copula()\]

## Examples

``` r
if (FALSE) { # \dontrun{
X_land <- hm_filter_land_cells(X)
fits   <- hm_fit_marginals(X_land)
pit    <- hm_pit_transform(X_land, fits)

cop_c <- hm_fit_copula(pit$U, type = "cvine")
cop_d <- hm_fit_copula(pit$U, type = "dvine")
cop_g <- hm_fit_copula(pit$U, type = "gaussian")
} # }
```
