# gamlssDS2 called by ds.gamlss

This is the second server-side aggregate function called by `ds.gamlss`.

## Usage

``` r
gamlssDS2(
  parameter = parameter,
  family = family,
  data = data,
  mu.beta.vect = mu.beta.vect,
  sigma.beta.vect = sigma.beta.vect,
  nu.beta.vect = nu.beta.vect,
  tau.beta.vect = tau.beta.vect,
  control = control,
  i.control = i.control
)
```

## Arguments

- parameter:

  A string specifing for which of the distribution parameters
  `c('mu', 'sigma', 'nu', 'tau')` the model fitting should be performed.

- family:

  A family string in the legal transmission format for DataSHIELD, which
  is used to define the distribution of the response variable. The
  DataSHIELD legal transmission format means that special characters,
  like '(' are replaced with the corresponding verbal descriptions, e.g.
  'left_parenthesis'. Currently, only the following families are
  supported:
  `family=c('NOleft_parenthesisright_parenthesis', 'NO2left_parenthesisright_parenthesis', 'BCCGleft_parenthesisright_parenthesis', 'BCPEleft_parenthesisright_parenthesis')`.

- data:

  A character string specifying a data.frame object holding the data to
  be analysed under the specified model.

- mu.beta.vect:

  A comma-separated string created by the client-side function
  specifying the vector of regression coefficients for mu at the current
  iteration.

- sigma.beta.vect:

  A comma-separated string created by the client-side function
  specifying the vector of regression coefficients for sigma at the
  current iteration.

- nu.beta.vect:

  A comma-separated string created by the client-side function
  specifying the vector of regression coefficients for nu at the current
  iteration.

- tau.beta.vect:

  A comma-separated string created by the client-side function
  specifying the vector of regression coefficients for tau at the
  current iteration.

- control:

  This sets the control parameters of the outer iterations algorithm
  using the
  [`gamlss.control`](https://rdrr.io/pkg/gamlss/man/gamlss.control.html)
  function. This is a comma-separated string of 7 numeric values: (i)
  c.crit (the convergence criterion for the algorithm), (ii) n.cyc (the
  number of cycles of the algorithm), (iii) mu.step (the step length for
  the parameter mu), (iv) sigma.step (the step length for the parameter
  sigma), (v) nu.step (the step length for the parameter nu), (vi)
  tau.step (the step length for the parameter tau), (vii) gd.tol (global
  deviance tolerance level). The default values for these 7 parameters
  are set to `control='0.001,20,1,1,1,1,Inf'`.

- i.control:

  This sets the control parameters of the inner iterations of the RS
  algorithm using the
  [`glim.control`](https://rdrr.io/pkg/gamlss/man/glim.control.html)
  function. This is a comma-separated string of 4 numeric values: (i) cc
  (the convergence criterion for the algorithm), (ii) cyc (the number of
  cycles of the algorithm), (iii) bf.cyc (the number of cycles of the
  backfitting algorithm), (iv) bf.tol (the convergence criterion
  (tolerance level) for the backfitting algorithm). The default values
  for these 4 parameters are set to `i.control='0.001,50,30,0.001'`.

## Value

A list with the following elements.

- `matrix`:

  Numeric matrix that can be aggregated on the client-side to obtain the
  updated weighted least squares estimator for the distribution
  parameter.

- `vector`:

  Numeric vector that can be aggregated on the client-side to obtain the
  updated weighted least squares estimator for the distribution
  parameter.

- `dv`:

  Numeric value for the new deviance on the server.

- `disclosure.risk`:

  Numeric value, either `0` or `1`, whereby `1` indicates a disclosure
  risk.

- `errorMessage2`:

  String for the disclosure risk. `errorMessage='No errors'` means that
  no disclosure risk was identified.

## Details

It is an aggregation function that uses the model structure and starting
parameter vectors constructed by `gamlssDS1` to iteratively obtain the
weighted least squares estimator for beta. The function `gamlssDS2` also
carries out a series of disclosure checks and if the arguments or data
fail any of those tests, model construction is blocked and an
appropriate server-side error message is created and returned to
`ds.gamlss` on the client-side. This function is not intended for direct
use by the user. For more details please see the extensive header of
`ds.gamlss`.

## Author

Annika Swenne
