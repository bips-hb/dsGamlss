# gamlssDS6 called by ds.gamlss

This is the sixth server-side aggregate function called by `ds.gamlss`.

## Usage

``` r
gamlssDS6(
  parameter = parameter,
  family = family,
  data = data,
  mu.beta.vect = mu.beta.vect,
  sigma.beta.vect = sigma.beta.vect,
  nu.beta.vect = nu.beta.vect,
  tau.beta.vect = tau.beta.vect,
  mu.gamma.vect = mu.gamma.vect,
  sigma.gamma.vect = sigma.gamma.vect,
  nu.gamma.vect = nu.gamma.vect,
  tau.gamma.vect = tau.gamma.vect,
  smoother.names = smoother.names,
  smoother.xl = smoother.xl,
  smoother.xr = smoother.xr,
  control = control,
  i.control = i.control,
  autostep = autostep,
  inner.iteration.count = inner.iteration.count,
  autostep.count = autostep.count
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

- mu.gamma.vect:

  A comma-separated string created by the client-side function
  specifying the vector of smoothing regression coefficients for mu at
  the current iteration.

- sigma.gamma.vect:

  A comma-separated string created by the client-side function
  specifying the vector of smoothing regression coefficients for sigma
  at the current iteration.

- nu.gamma.vect:

  A comma-separated string created by the client-side function
  specifying the vector of smoothing regression coefficients for nu at
  the current iteration.

- tau.gamma.vect:

  A comma-separated string created by the client-side function
  specifying the vector of smoothing regression coefficients for tau at
  the current iteration.

- smoother.names:

  A string vector specifying the unique variable names for the smoother.

- smoother.xl:

  A comma-separated string created by the client-side function
  specifying the left boundary for the knots for the smoother in
  `smoother.names`.

- smoother.xr:

  A comma-separated string created by the client-side function
  specifying the right boundary for the knots for the smoother in
  `smoother.names`.

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

- autostep:

  Logical, indicating whether the steps should be halved automatically
  if the new global deviance is greater than the old one. The default is
  `autostep=TRUE`.

- inner.iteration.count:

  An integer with the iteration count for the inner iteration.

- autostep.count:

  An integer with the iteration count for the autostep iteration.

## Value

A list with the following elements.

- `dv`:

  Numeric value for the new deviance on the server.

- `pen`:

  Numeric value with the sum of the quadratic penalties for the
  `parameter`.

- `errorMessage3`:

  String for the validity of the fitted values in the inner iteration.
  `errorMessage3='Fitted values in the inner iteration out of range.'`
  indicates a problem, whereas `errorMessage3=NULL` means no problem
  with the fitted values.

## Details

It is an aggregation function that updates the distribution parameters
and returns the current deviance. The deviance can then be used on the
client side to check whether the inner iteration converged. This
function is not intended for direct use by the user. For more details
please see the extensive header of `ds.gamlss`.

## Author

Annika Swenne
