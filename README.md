
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsGamlss

<!-- badges: start -->
<!-- badges: end -->

The `dsGamlss` package is a [DataSHIELD](https://www.datashield.org)
server-side package that includes the server-side functions to fit
Generalized Additive Models for Location, Scale and Shape (GAMLSS) \[1\]
using DataSHIELD. It is based on the original
[gamlss](https://cran.r-project.org/package=gamlss) implementation \[1\]
and the [dsBase](https://github.com/datashield/dsBase) package \[2\].

### DataSHIELD

DataSHIELD is a software infrastructure which allows you to do
non-disclosive federated analysis on sensitive data. The
[DataSHIELD](https://www.datashield.org) website has in depth
descriptions of what it is, how it works and how to install it. A key
point to highlight is that DataSHIELD has a client-server
infrastructure, so the `dsGamlss` package needs to be used in
conjunction with the
[dsGamlssClient](https://github.com/bips-hb/dsGamlssClient) package -
trying to use one without the other makes no sense. Detailed
instructions on how to install DataSHIELD can be found at the DataSHIELD
[Wiki](https://www.datashield.org/wiki). Discussion and help with using
DataSHIELD can be obtained from the DataSHIELD
[Forum](https://datashield.discourse.group/)

## Installation

You can install the development version of `dsGamlss` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bips-hb/dsGamlss")
```

## Example

Direct usage of the functions in this package is not recommended.
Therefore, no example is provided. Instead, you should use the functions
from the corresponding client-side package
[dsGamlssClient](https://github.com/bips-hb/dsGamlssClient) to fit
GAMLSS models with DataSHIELD and the `dsGamlssClient` package provides
examples for that.

## References

1.  Rigby RA, Stasinopoulos DM. Generalized additive models for
    location, scale and shape. Journal of the Royal Statistical Society:
    Series C (Applied Statistics). 2005;54(3):507-54.
2.  DataSHIELD Developers (2023). *dsBaseClient: DataSHIELD Client
    Functions*. R package version 6.3.0.
