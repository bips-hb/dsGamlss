# Local function that creates the basis matrix for p-splines

This is an internal function called by the functions 'gamlssDS3',
'gamlssDS4', 'gamlssDS5' and 'gamlssDS6'.

## Usage

``` r
bbase(x, xl, xr, ndx = 20, deg = 3)
```

## Arguments

- x:

  A vector of x-values.

- xl:

  Left limit for the knots.

- xr:

  Right limit for the knots.

- ndx:

  Number of equal space intervals in x.

- deg:

  Degree of the polynomial.

## Value

A B-basis matrix for p-splines.

## Details

The function is based on the internal bbase-function defined in
gamlss::pb (version 5.4-12). Note that unlike in the original function
the use of quantiles is currently not implemented in bbase

## Author

Mikis Stasinopoulos, Bob Rigby and Paul Eilers
