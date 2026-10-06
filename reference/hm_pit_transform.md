# Apply the Probability Integral Transform

Transforms `X` to uniform margins using each cell's fitted CDF, then
optionally applies the probit transform to obtain centred Gaussian
scores. Implements Equations (3)–(4) from Clavijo Mesa et al. (2026, SF
paper).

## Usage

``` r
hm_pit_transform(X, fits, gaussian = TRUE)
```

## Arguments

- X:

  Numeric matrix (\\T \times N\\) from \[hm_to_matrix()\].

- fits:

  data.frame from \[hm_fit_marginals()\].

- gaussian:

  If `TRUE` (default), also returns centred Gaussian scores \\Z =
  \Phi^{-1}(U)\\.

## Value

A list with `U` (uniform matrix), `Z` (Gaussian scores or `NULL`),
`n_cells_ok`, and `n_cells_skipped`.

## See also

\[hm_fit_marginals()\], \[hm_marginal_gof()\]
