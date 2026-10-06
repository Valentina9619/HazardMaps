# Decode the NetCDF time axis

Converts CF-style time values (e.g. `"days since 1970-01-01"`) to an R
`Date` or `POSIXct` vector.

## Usage

``` r
hm_decode_time(x, tz = "UTC")
```

## Arguments

- x:

  An `hm_hazard` object returned by \[hm_read_netcdf()\].

- tz:

  Time zone string for sub-daily data. Defaults to `"UTC"`.

## Value

The same object with `x$time` replaced by a `Date` vector (daily data)
or `POSIXct` vector (sub-daily data).

## See also

\[hm_read_netcdf()\], \[hm_to_matrix()\]
