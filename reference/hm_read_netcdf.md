# Read a hazard NetCDF dataset

Reads a NetCDF file containing a hazard variable on a spatial grid
(regular or rotated) and returns an \`hm_hazard\` object. CRS is not
required.

## Usage

``` r
hm_read_netcdf(file, var = NULL)
```

## Arguments

- file:

  Path to a NetCDF file.

- var:

  Name of the hazard variable. If `NULL` and the file contains exactly
  one variable, it is selected automatically.

## Value

An object of class `hm_hazard`.

## Details

Supports CF-convention NetCDF files with regular grids (e.g. Copernicus,
ERA5) and rotated or curvilinear grids (e.g. EURO-CORDEX). Coordinate
names `lat`/`lon` and `latitude`/`longitude` are detected automatically.

## See also

\[hm_standardize_coords()\], \[hm_decode_time()\]
