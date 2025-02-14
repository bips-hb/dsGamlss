#'
#' @title gamlssDS3 called by ds.gamlss
#' @description This is the third server-side aggregate function called by \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @details It is an aggregation function that uses the model structure and starting
#' parameter vectors constructed by \code{gamlssDS1} to iteratively obtain the penalized weighted least squares estimator for gamma.
#' This function is not intended for direct use by the user. For more details please see the extensive header of \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @param parameter A string specifing for which of the distribution parameters \code{c('mu', 'sigma', 'nu', 'tau')}
#' the model fitting should be performed.
#' @param smoother An integer indicating the number of the smoother for the \code{parameter} that should be fitted.
#' @param family A family string in the legal transmission format for DataSHIELD, which
#' is used to define the distribution of the response variable. The DataSHIELD legal transmission
#' format means that special characters, like '(' are replaced with the corresponding verbal
#' descriptions, e.g. 'left_parenthesis'. Currently, only the following families are supported:
#' \code{family=c('NOleft_parenthesisright_parenthesis', 'NO2left_parenthesisright_parenthesis',
#' 'BCCGleft_parenthesisright_parenthesis', 'BCPEleft_parenthesisright_parenthesis')}.
#' @param data A character string specifying a data.frame object holding
#' the data to be analysed under the specified model.
#' @param mu.beta.vect A comma-separated string created by the client-side function specifying the
#' vector of regression coefficients for mu at the current iteration.
#' @param sigma.beta.vect A comma-separated string created by the client-side function specifying the
#' vector of regression coefficients for sigma at the current iteration.
#' @param nu.beta.vect A comma-separated string created by the client-side function specifying the
#' vector of regression coefficients for nu at the current iteration.
#' @param tau.beta.vect A comma-separated string created by the client-side function specifying the
#' vector of regression coefficients for tau at the current iteration.
#' @param mu.gamma.vect A comma-separated string created by the client-side function specifying the
#' vector of smoothing regression coefficients for mu at the current iteration.
#' @param sigma.gamma.vect A comma-separated string created by the client-side function specifying the
#' vector of smoothing regression coefficients for sigma at the current iteration.
#' @param nu.gamma.vect A comma-separated string created by the client-side function specifying the
#' vector of smoothing regression coefficients for nu at the current iteration.
#' @param tau.gamma.vect A comma-separated string created by the client-side function specifying the
#' vector of smoothing regression coefficients for tau at the current iteration.
#' @param smoother.names A string vector specifying the unique variable names for the smoother.
#' @param smoother.xl A comma-separated string created by the client-side function specifying the left
#' boundary for the knots for the smoother in \code{smoother.names}.
#' @param smoother.xr A comma-separated string created by the client-side function specifying the right
#' boundary for the knots for the smoother in \code{smoother.names}.
#' @param control This sets the control parameters of the outer iterations algorithm
#' using the gamlss.control function. This is a comma-separated string of 7 numeric values:
#' (i) c.crit (the convergence criterion for the algorithm), (ii) n.cyc (the number of cycles of
#' the algorithm), (iii) mu.step (the step length for the parameter mu), (iv) sigma.step
#' (the step length for the parameter sigma), (v) nu.step (the step length for the
#' parameter nu), (vi) tau.step (the step length for the parameter tau), (vii) gd.tol
#' (global deviance tolerance level). The default values for these 7 parameters are
#' set to \code{control='0.001,20,1,1,1,1,Inf'}.
#' @param i.control This sets the control parameters of the inner iterations of the
#' RS algorithm using the glim.control function. This is a comma-separated string of 4 numeric values:
#' (i) cc (the convergence criterion for the algorithm), (ii) cyc (the number of
#' cycles of the algorithm), (iii) bf.cyc (the number of cycles of the backfitting
#' algorithm), (iv) bf.tol (the convergence criterion (tolerance level) for the
#' backfitting algorithm). The default values for these 4 parameters are set to
#' \code{i.control='0.001,50,30,0.001'}.
#' @return A list with the following elements.
#' \describe{
#'  \item{\code{matrix}}{Numeric matrix that can be aggregated on the client-side to obtain the updated penalized weighted least squares estimator for the distribution parameter.}
#'  \item{\code{vector}}{Numeric vector that can be aggregated on the client-side to obtain the updated penalized weighted least squares estimator for the distribution parameter.}
#'  \item{\code{sumofsquares}}{Numeric value for the sum of squares for the old and new smoothing fitted values. This is needed on the client-side to determine the stopping criterion for backfitting.}
#'  \item{\code{sumofweights}}{Numeric value for the sum of weights. This is needed on the client-side to determine the stopping criterion for backfitting.}
#' }
#' @author Annika Swenne
#' @import gamlss.dist
#' @export
#'

gamlssDS3 <- function(parameter = parameter, smoother = smoother, family = family, data = data,
                      mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
                      nu.beta.vect = nu.beta.vect, tau.beta.vect = tau.beta.vect,
                      mu.gamma.vect = mu.gamma.vect, sigma.gamma.vect = sigma.gamma.vect,
                      nu.gamma.vect = nu.gamma.vect, tau.gamma.vect = tau.gamma.vect,
                      smoother.names = smoother.names,
                      smoother.xl = smoother.xl, smoother.xr = smoother.xr,
                      control = control, i.control = i.control) {
  #**************************************************************************
  # I) Preparation ----
  # Reconvert the transfer strings into required variable types
  # Extract the varnames
  # Get function to calculate deviance increment
  #**************************************************************************

  ## Get the value of the 'data' parameter provided as character on the client side
  dataname <- data
  if (is.null(dataname)) {
    data <- NULL
  } else {
    data <- eval(parse(text = dataname), envir = parent.frame())
  }
  Ntotal <- dim(data)[1]

  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text = family), envir = environment()))

  c1 <- as.numeric(unlist(strsplit(control, split = ",")))
  c2 <- as.numeric(unlist(strsplit(i.control, split = ",")))

  # Convert parameter vectors from transmittable (character) format to numeric
  # beta
  if (!(is.null(mu.beta.vect))) {
    mu.beta.vect <- as.numeric(unlist(strsplit(mu.beta.vect, split = ",")))
  } else {
    mu.beta.vect <- NULL
  }

  if (!(is.null(sigma.beta.vect))) {
    sigma.beta.vect <- as.numeric(unlist(strsplit(sigma.beta.vect, split = ",")))
  } else {
    sigma.beta.vect <- NULL
  }

  if (!(is.null(nu.beta.vect))) {
    nu.beta.vect <- as.numeric(unlist(strsplit(nu.beta.vect, split = ",")))
  } else {
    nu.beta.vect <- NULL
  }

  if (!(is.null(tau.beta.vect))) {
    tau.beta.vect <- as.numeric(unlist(strsplit(tau.beta.vect, split = ",")))
  } else {
    tau.beta.vect <- NULL
  }

  # gamma
  if (!(is.null(mu.gamma.vect))) {
    mu.gamma.vect <- as.numeric(unlist(strsplit(mu.gamma.vect, split = ",")))
  } else {
    mu.gamma.vect <- NULL
  }

  if (!(is.null(sigma.gamma.vect))) {
    sigma.gamma.vect <- as.numeric(unlist(strsplit(sigma.gamma.vect, split = ",")))
  } else {
    sigma.gamma.vect <- NULL
  }

  if (!(is.null(nu.gamma.vect))) {
    nu.gamma.vect <- as.numeric(unlist(strsplit(nu.gamma.vect, split = ",")))
  } else {
    nu.gamma.vect <- NULL
  }

  if (!(is.null(tau.gamma.vect))) {
    tau.gamma.vect <- as.numeric(unlist(strsplit(tau.gamma.vect, split = ",")))
  } else {
    tau.gamma.vect <- NULL
  }

  # get the beta and gammma vectors for the respective parameter
  beta.vect <- eval(parse(text = paste(parameter, ".beta.vect", sep = "")), envir = environment())
  gamma.vect <- eval(parse(text = paste(parameter, ".gamma.vect", sep = "")), envir = environment())

  # Convert knot boundaries from transmittable (character) format to numeric
  if (!(is.null(smoother.xl))) {
    smoother.xl <- as.numeric(unlist(strsplit(smoother.xl, split = ",")))
  } else {
    smoother.xl <- NULL
  }
  if (!(is.null(smoother.xr))) {
    smoother.xr <- as.numeric(unlist(strsplit(smoother.xr, split = ",")))
  } else {
    smoother.xr <- NULL
  }

  #**************************************************************************
  # II) Calculate matrix & vector to return to client ----
  #**************************************************************************

  #* A) Get fitted model ----
  mod.gamlss.ds <- base::get("temp_mod.gamlss.ds", env = parent.frame())

  ## get design matrix for the parameter
  X.mat <- as.matrix(eval(parse(text = paste("mod.gamlss.ds$", parameter, ".x", sep = "")), envir = environment()))
  y <- as.vector(mod.gamlss.ds$y)

  ## get design matrix for the smoothers
  # get the control parameters for the smoothers
  coefficients <- eval(parse(text = paste("names(mod.gamlss.ds$", parameter, ".coefficients)", sep = "")), envir = environment())
  smoother.coef <- coefficients[grep(pattern = "pb(", x = tolower(coefficients), fixed = TRUE)]
  # only keep the arguments for the pb() function
  pb.args <- substr(smoother.coef, start = 4, stop = nchar(smoother.coef) - 1)
  pb.args <- strsplit(pb.args, split = ",", fixed = TRUE)

  # create design matrices for current & previous smoother if possible
  name <- eval(parse(text = paste("mod.gamlss.ds$", parameter, ".coefSmo[[", smoother, "]]$name", sep = "")), envir = environment())
  if (length(grep(pattern = "pb.control", x = pb.args[[smoother]], fixed = TRUE)) > 0) {
    # control parameters specified
    pb.control <- eval(parse(text = pb.args[[smoother]][grep(pattern = "pb.control", x = pb.args[[smoother]], fixed = TRUE)]))
  } else {
    # no control parameters specified - use default
    pb.control <- eval(parse(text = "pb.control()"))
  }
  x <- eval(parse(text = name), envir = parent.frame())
  Z.mat <- bbase(
    x = x, xl = smoother.xl[which(smoother.names == name)], xr = smoother.xr[which(smoother.names == name)],
    ndx = pb.control$inter, deg = pb.control$degree
  )
  # get design matrix for previous smoother if possible
  if (smoother > 1) {
    name <- eval(parse(text = paste("mod.gamlss.ds$", parameter, ".coefSmo[[", smoother - 1, "]]$name", sep = "")), envir = environment())
    if (length(grep(pattern = "pb.control", x = pb.args[[smoother - 1]], fixed = TRUE)) > 0) {
      # control parameters specified
      pb.control <- eval(parse(text = pb.args[[smoother - 1]][grep(pattern = "pb.control", x = pb.args[[smoother - 1]], fixed = TRUE)]))
    } else {
      # no control parameters specified - use default
      pb.control <- eval(parse(text = "pb.control()"))
    }
    x <- eval(parse(text = name), envir = parent.frame())
    Z.mat.old <- bbase(
      x = x, xl = smoother.xl[which(smoother.names == name)], xr = smoother.xr[which(smoother.names == name)],
      ndx = pb.control$inter, deg = pb.control$degree
    )
  }

  ## Get fitted values for all distribution parameters
  # necessary for all parameters to calculate deviance
  if ("mu" %in% names(family$parameters)) {
    mu <- base::get("temp_mu", env = parent.frame())
  }
  if ("sigma" %in% names(family$parameters)) {
    sigma <- base::get("temp_sigma", env = parent.frame())
  }
  if ("nu" %in% names(family$parameters)) {
    nu <- base::get("temp_nu", env = parent.frame())
  }
  if ("tau" %in% names(family$parameters)) {
    tau <- base::get("temp_tau", env = parent.frame())
  }

  ## Calculate smoothing fitted value matrix s
  gamma.start <- 1
  coefSmo <- eval(parse(text = paste("mod.gamlss.ds$", parameter, ".coefSmo", sep = "")), envir = environment())
  s.old <- base::get(paste("temp_", parameter, ".s", sep = ""), env = parent.frame())
  # get the gamma vectors for the respective parameter & multiply them with the matrices
  for (i in 1:length(coefSmo)) {
    gamma.length <- dim(coefSmo[[i]]$coef)[1]
    gamma.end <- gamma.start + gamma.length - 1
    gamma <- gamma.vect[gamma.start:gamma.end]
    # calculate new smoothing fitted values for previous smoother
    if (i == (smoother - 1)) {
      s.update <- as.vector(Z.mat.old %*% gamma)
    }
    gamma.start <- gamma.end + 1
  }
  s <- s.old
  # update smoothing fitted value matrix for previous smoother
  if (smoother > 1) {
    s[, (smoother - 1)] <- s.update
    base::assign(paste("temp_", parameter, ".s", sep = ""), s, env = parent.frame())
  }

  #* B) Update vectors----
  ## Calculate predictor vector eta for the parameter
  eta <- eval(parse(text = paste("family$", parameter, ".linkfun(", parameter, ")", sep = "")), envir = environment())

  ## Calculate score and weights (for Fisher-scoring algorithm)
  dr <- eval(parse(text = paste("family$", parameter, ".dr(eta)", sep = "")), envir = environment()) # dparameter/ deta
  dr <- 1 / dr # deta/ dparameter = 1/ (dparameter/ deta)

  # get the first and second derivatives of the log-likelihood
  if (parameter == "mu") {
    fv <- mu
    # first derivative of log-likelihood
    dldp.function <- family$dldm
    formals(dldp.function, envir = new.env()) <- alist(mu = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldm2
    formals(d2ldp2.function, envir = new.env()) <- alist(mu = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter == "sigma") {
    fv <- sigma
    # first derivative of log-likelihood
    dldp.function <- family$dldd
    formals(dldp.function, envir = new.env()) <- alist(sigma = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldd2
    formals(d2ldp2.function, envir = new.env()) <- alist(sigma = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter == "nu") {
    fv <- nu
    # first derivative of log-likelihood
    dldp.function <- family$dldv
    formals(dldp.function, envir = new.env()) <- alist(nu = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldv2
    formals(d2ldp2.function, envir = new.env()) <- alist(nu = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter == "tau") {
    fv <- tau
    # first derivative of log-likelihood
    dldp.function <- family$dldt
    formals(dldp.function, envir = new.env()) <- alist(tau = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldt2
    formals(d2ldp2.function, envir = new.env()) <- alist(tau = fv) # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }

  dldp <- dldp.function(fv) # first derivative of log-likelihood with respect to parameter
  d2ldp2 <- d2ldp2.function(fv) # expected  second derivative of log-Likelihood with respect to parameter
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2 / (dr * dr)) # weights for Fisher-scoring algorithm =-(d2l/d2p)(dparameter/deta)^2
  # we need to stop the weights to go to Infty
  wt <- ifelse(wt > 1e+10, 1e+10, wt)
  wt <- ifelse(wt < 1e-10, 1e-10, wt)

  ## stopping criterion for backfitting
  # the sums are necessary to calculate deltaf on the server which is needed to determine
  # the stoppping criterion for backfitting
  if (smoother == 1) {
    sumofsquares <- 0
    sumofweights <- 1
  } else {
    sumofsquares <- sum((s[, (smoother - 1)] - s.old[, (smoother - 1)])^2 * wt)
    sumofweights <- sum(wt)
  }

  ## Update working variable vector wv
  wv <- eta + dldp / (dr * wt)
  if (family$type == "Mixed") {
    wv <- ifelse(is.nan(wv), 0, wv)
  }

  #* C) Calculate matrix & vectors ----

  ## Calculate partial residuals
  partial.residuals <- wv - X.mat %*% beta.vect - base::rowSums(as.matrix(s[, -smoother]))

  ## Calculate matrix and vector to return to the client
  vector <- t(Z.mat) %*% (wt * partial.residuals)
  matrix <- t(Z.mat) %*% diag(wt) %*% Z.mat
  # remove the dimnames attributes
  attr(vector, "dimnames") <- NULL
  attr(matrix, "dimnames") <- NULL


  #**************************************************************************
  # III) Output ----
  #**************************************************************************

  return(list(
    matrix = matrix, vector = vector, sumofsquares = sumofsquares,
    sumofweights = sumofweights
  ))
}
# AGGREGATE FUNCTION
# gamlssDS3
