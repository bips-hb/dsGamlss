# gamlssDS1 called by ds.gamlss

This is the first server-side aggregate function called by `ds.gamlss`.

## Usage

``` r
gamlssDS1(
  formula = formula,
  sigma.formula = sigma.formula,
  nu.formula = nu.formula,
  tau.formula = tau.formula,
  family = family,
  data = data,
  mu.fix = mu.fix,
  sigma.fix = sigma.fix,
  nu.fix = nu.fix,
  tau.fix = tau.fix,
  global.mean = global.mean,
  global.sd = global.sd,
  control = control,
  i.control = i.control,
  autostep = autostep
)
```

## Arguments

- formula:

  A formula string in the legal transmission format for DataSHIELD,
  specifying the model for the mu distribution parameter. The DataSHIELD
  legal transmission format means that special characters, like '(' are
  replaced with the corresponding verbal descriptions, e.g.
  'left_parenthesis'.

- sigma.formula:

  A formula string in the legal transmission format for DataSHIELD,
  specifying the model for the sigma distribution parameter. The
  DataSHIELD legal transmission format means that special characters,
  like '(' are replaced with the corresponding verbal descriptions, e.g.
  'left_parenthesis'.

- nu.formula:

  A formula string in the legal transmission format for DataSHIELD,
  specifying the model for the nu distribution parameter. The DataSHIELD
  legal transmission format means that special characters, like '(' are
  replaced with the corresponding verbal descriptions, e.g.
  'left_parenthesis'.

- tau.formula:

  A formula string in the legal transmission format for DataSHIELD,
  specifying the model for the tau distribution parameter. The
  DataSHIELD legal transmission format means that special characters,
  like '(' are replaced with the corresponding verbal descriptions, e.g.
  'left_parenthesis'.

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

- mu.fix:

  Logical, indicating whether the mu parameter should be kept fixed in
  the fitting processes. Default `mu.fix=FALSE`.

- sigma.fix:

  Logical, indicating whether the sigma parameter should be kept fixed
  in the fitting processes. Default `sigma.fix=FALSE`.

- nu.fix:

  Logical, indicating whether the nu parameter should be kept fixed in
  the fitting processes. Default `nu.fix=FALSE`.

- tau.fix:

  Logical, indicating whether the tau parameter should be kept fixed in
  the fitting processes. Default `tau.fix=FALSE`.

- global.mean:

  Numeric value, giving the global mean of the outcome variable, which
  is necessary to initialize the distribution parameters for some
  families. Otherwise `global.mean=NULL`.

- global.sd:

  Numeric value, giving the global sd of the outcome variable which is
  necessary to initialize the distribution parameters for some families.
  Otherwise `global.sd=NULL`.

- control:

  This sets the control parameters of the outer iterations algorithm
  using the using the
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
  algorithm using the using the
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

## Value

A list with the following elements.

- `mod.gamlss.ds`:

  A `gamlss` object with all components as in the native R
  [`gamlss`](https://rdrr.io/pkg/gamlss/man/gamlss.html) function.
  Individual-level information like the components `y` (the response)
  and `residuals` (the normalised quantile residuals of the model) are
  not disclosed to the client-side.

- `G.dev`:

  Numeric value for the initial deviance on the server.

- `dim.mu.x`:

  Numeric vector with two elements, specifying the dimension of the
  design matrix for mu.

- `dim.sigma.x`:

  Numeric vector with two elements, specifying the dimension of the
  design matrix for sigma.

- `dim.nu.x`:

  Numeric vector with two elements, specifying the dimension of the
  design matrix for nu.

- `dim.tau.x`:

  Numeric vector with two elements, specifying the dimension of the
  design matrix for tau.

- `smoother.names`:

  String vector with the unique variable names that are used for
  smoothing.

- `smoother.xmin`:

  Numeric vector with anononymized minima for the variables in
  `smoother.names`.

- `smoother.xmax`:

  Numeric vector with anononymized maxima for the variables in
  `smoother.names`.

- `y.invalid`:

  Numeric value, either `0` or `1`, whereby `1` indicates a disclosure
  risk in the response variable.

- `mu.par.invalid`:

  Numeric vector with elements `0` or `1`, whereby `1` indicates a
  disclosure risk in the corresponding explanatory variable for mu.

- `sigma.par.invalid`:

  Numeric vector with elements `0` or `1`, whereby `1` indicates a
  disclosure risk in the corresponding explanatory variable for sigma.

- `nu.par.invalid`:

  Numeric vector with elements `0` or `1`, whereby `1` indicates a
  disclosure risk in the corresponding explanatory variable for nu.

- `tau.par.invalid`:

  Numeric vector with elements `0` or `1`, whereby `1` indicates a
  disclosure risk in the corresponding explanatory variable for tau.

- `gamlss.saturation.invalid`:

  Numeric value, either `0` or `1`, whereby `1` indicates a disclosure
  risk from an oversaturated model.

- `errorMessage`:

  String for the disclosure risk.
  `errorMessage='Study data or applied model invalid for this source'`
  indicates a disclosure risk, whereas `errorMessage='No errors'` means
  that no disclosure risk was identified.

## Details

It is an aggregation function that sets up the appropriate model
structure and dimensions to fit a `ds.gamlss` model. This function is
not intended for direct use by the user. For more details please see the
extensive header of `ds.gamlss`.

## Author

Annika Swenne
