#'
#' @title gamlssDS4 called by ds.gamlss
#' @description This is the fourth serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that returns the required inner product
#' to estimate the smoothing parameter lambda with the ML method.
#' For more details please see the extensive header of ds.gamlss and also the
#' gamlss function in native R gamlss package.
#' @param parameter a string specifing for which of the model parameters "mu", "sigma", "nu"
#' or "tau" the model fitting should be performed
#' @param smoother an integer indicating the number of the smoother that should be fitted
#' @param family a gamlss.family object, which is used to define the distribution 
#' and the link functions of the various parameters. The distribution families 
#' supported by gamlss() can be found in gamlss.family. Functions such as BI() 
#' (binomial) produce a family object. Also can be given without the parentheses
#' i.e. BI. Family functions can take arguments, as in BI(mu.link=probit).
#' @param data an optional character string specifying a data.frame object holding 
#' the data to be analysed under the specified model 
#' @param mu.beta.vect a numeric vector created by the clientside function specifying the
#' vector of regression coefficients for mu at the current iteration.
#' @param sigma.beta.vect a numeric vector created by the clientside function specifying the
#' vector of regression coefficients for sigma at the current iteration.
#' @param nu.beta.vect a numeric vector created by the clientside function specifying the
#' vector of regression coefficients for nu at the current iteration.
#' @param tau.beta.vect a numeric vector created by the clientside function specifying the
#' vector of regression coefficients for tau at the current iteration.
#' @param mu.gamma.vect a numeric vector created by the clientside function specifying the
#' vector of smoothing regression coefficients for mu at the current iteration.
#' @param sigma.gamma.vect a numeric vector created by the clientside function specifying the
#' vector of smoothing regression coefficients for sigma at the current iteration.
#' @param nu.gamma.vect a numeric vector created by the clientside function specifying the
#' vector of smoothing regression coefficients for nu at the current iteration.
#' @param tau.gamma.vect a numeric vector created by the clientside function specifying the
#' vector of smoothing regression coefficients for tau at the current iteration.
#' @param smoother.names a character vector specifying the unique variable names for the smoother.
#' @param smoother.xl a numeric vector created by the clientside function specifying the left
#' boundary for the knots for the pb-smoother.
#' @param smoother.xr a numeric vector created by the clientside function specifying the right
#' boundary for the knots for the pb-smoother.
#' @param control this sets the control parameters of the outer iterations algorithm 
#' using the gamlss.control function. This is a vector of 7 numeric values: (i) c.crit 
#' (the convergence criterion for the algorithm), (ii) n.cyc (the number of cycles of 
#' the algorithm), (iii) mu.step (the step length for the parameter mu), (iv) sigma.step 
#' (the step length for the parameter sigma), (v) nu.step (the step length for the
#' parameter nu), (vi) tau.step (the step length for the parameter tau), (vii) gd.tol
#' (global deviance tolerance level). The default values for these 7 parameters are 
#' set to c(0.001, 20, 1, 1, 1, 1, Inf).
#' @param i.control this sets the control parameters of the inner iterations of the 
#' RS algorithm using the glim.control function. This is a vector of 4 numeric values: 
#' (i) cc (the convergence criterion for the algorithm), (ii) cyc (the number of 
#' cycles of the algorithm), (iii) bf.cyc (the number of cycles of the backfitting 
#' algorithm), (iv) bf.tol (the convergence criterion (tolerance level) for the 
#' backfitting algorithm). The default values for these 4 parameters are set to 
#' c(0.001, 50, 30, 0.001).
#' @return a gamlss object with all components as in the native R gamlss function. 
#' Individual-level information like the components y (the response response) and 
#' residuals (the normalised quantile residuals of the model) are not disclosed to 
#' the client-side.
#' @author Annika Swenne
#' @import gamlss.dist
#' @export
#'

gamlssDS4 <- function(parameter = parameter, smoother = smoother, family = family, data = data, 
                      mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
                      nu.beta.vect = nu.beta.vect, tau.beta.vect = tau.beta.vect,
                      mu.gamma.vect = mu.gamma.vect, sigma.gamma.vect = sigma.gamma.vect,
                      nu.gamma.vect = nu.gamma.vect, tau.gamma.vect = tau.gamma.vect,
                      smoother.names = smoother.names,
                      smoother.xl = smoother.xl, smoother.xr = smoother.xr,
                      control = control, i.control = i.control){
  
  #**************************************************************************
  # I) Preparation ---- 
  # Reconvert the transfer strings into required variable types
  # Extract the varnames
  # Get function to calculate deviance increment
  #**************************************************************************
  
  ## Get the value of the 'data' parameter provided as character on the client side
  dataname <- data
  if(is.null(dataname)){
    data <- NULL 
  }else{
    data <- eval(parse(text=dataname), envir=parent.frame())
  }
  Ntotal <- dim(data)[1]
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family), envir=environment()))
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  # Convert parameter vectors from transmittable (character) format to numeric 
  # beta
  if (!(is.null(mu.beta.vect))){
    mu.beta.vect <- as.numeric(unlist(strsplit(mu.beta.vect, split=",")))
  } else {
    mu.beta.vect <- NULL
  }
  
  if (!(is.null(sigma.beta.vect))){
    sigma.beta.vect <- as.numeric(unlist(strsplit(sigma.beta.vect, split=",")))
  } else {
    sigma.beta.vect <- NULL
  }
  
  if (!(is.null(nu.beta.vect))){
    nu.beta.vect <- as.numeric(unlist(strsplit(nu.beta.vect, split=",")))
  } else {
    nu.beta.vect <- NULL
  }
  
  if (!(is.null(tau.beta.vect))){
    tau.beta.vect <- as.numeric(unlist(strsplit(tau.beta.vect, split=",")))
  } else {
    tau.beta.vect <- NULL
  }
  
  #gamma
  if (!(is.null(mu.gamma.vect))){
    mu.gamma.vect <- as.numeric(unlist(strsplit(mu.gamma.vect, split=",")))
  } else {
    mu.gamma.vect <- NULL
  }
  
  if (!(is.null(sigma.gamma.vect))){
    sigma.gamma.vect <- as.numeric(unlist(strsplit(sigma.gamma.vect, split=",")))
  } else {
    sigma.gamma.vect <- NULL
  }
  
  if (!(is.null(nu.gamma.vect))){
    nu.gamma.vect <- as.numeric(unlist(strsplit(nu.gamma.vect, split=",")))
  } else {
    nu.gamma.vect <- NULL
  }
  
  if (!(is.null(tau.gamma.vect))){
    tau.gamma.vect <- as.numeric(unlist(strsplit(tau.gamma.vect, split=",")))
  } else {
    tau.gamma.vect <- NULL
  }
  
  # get the beta and gammma vectors for the respective parameter
  beta.vect <- eval(parse(text=paste(parameter, ".beta.vect", sep="")), envir=environment())
  gamma.vect <- eval(parse(text=paste(parameter, ".gamma.vect", sep="")), envir=environment())
  
  # Convert knot boundaries from transmittable (character) format to numeric
  if (!(is.null(smoother.xl))){
    smoother.xl <- as.numeric(unlist(strsplit(smoother.xl, split=",")))
  } else {
    smoother.xl <- NULL
  }
  if (!(is.null(smoother.xr))){
    smoother.xr <- as.numeric(unlist(strsplit(smoother.xr, split=",")))
  } else {
    smoother.xr <- NULL
  }
  
  #**************************************************************************
  # II) Calculate matrix & vector to return to client ----  
  #**************************************************************************
  
  #*A) Get fitted model ----
  mod.gamlss.ds <- base::get("temp_mod.gamlss.ds", env=parent.frame())
  
  ## get design matrix for the parameter
  X.mat <- as.matrix(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".x", sep="")), envir=environment()))
  y <- as.vector(mod.gamlss.ds$y)
  
  ## get design matrix for the smoothers
  # get the control parameters for the smoothers
  coefficients <- eval(parse(text=paste("names(mod.gamlss.ds$", parameter, ".coefficients)", sep="")), envir=environment())
  smoother.coef <- coefficients[grep(pattern="pb(", x=tolower(coefficients), fixed=TRUE)]
  # only keep the arguments for the pb() function
  pb.args <- substr(smoother.coef, start=4, stop=nchar(smoother.coef)-1)
  pb.args <- strsplit(pb.args, split=",", fixed=TRUE)
  
  # create design matrices for all smoothers if possible
  for (i in 1:length(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), envir=environment()))){
    name <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", i, "]]$name", sep="")), envir=environment())
    if (length(grep(pattern="pb.control", x=pb.args[[i]], fixed=TRUE))>0) {
      # control parameters specified
      pb.control <- eval(parse(text=pb.args[[i]][grep(pattern="pb.control", x=pb.args[[i]], fixed=TRUE)]))
    } else {
      # no control parameters specified - use default
      pb.control <- eval(parse(text="pb.control()"))
    }
    x <- eval(parse(text=name), envir=parent.frame())
    basismatrix <- bbase(x=x, xl=smoother.xl[which(smoother.names==name)], xr=smoother.xr[which(smoother.names==name)],
                         ndx=pb.control$inter, deg=pb.control$degree)
    base::assign(paste("temp_Z", i, ".mat", sep=""), basismatrix, env=environment())
    if (i==smoother){
      nobs <- nrow(basismatrix)
      Z.mat <- basismatrix
    }
  }
  
  ## Get fitted values for all distribution parameters
  # necessary for all parameters to calculate deviance
  if("mu" %in% names(family$parameters)){
    mu <- base::get("temp_mu", env=parent.frame())
  }
  if("sigma" %in% names(family$parameters)){
    sigma <- base::get("temp_sigma", env=parent.frame())
  }
  if("nu" %in% names(family$parameters)){
    nu <- base::get("temp_nu", env=parent.frame())
  }
  if("tau" %in% names(family$parameters)){
    tau <- base::get("temp_tau", env=parent.frame())
  }

  ## Get smoothing fitted value matrix s
  coefSmo <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), envir=environment())
  if (!is.null(coefSmo)){
    s <- base::get(paste("temp_", parameter, ".s", sep=""), env=parent.frame())
  }
  
  ## Get gamma vector
  gamma.start <- 1
  for (i in 1:length(coefSmo)){
    gamma.length <- dim(coefSmo[[i]]$coef)[1]
    gamma.end <- gamma.start+gamma.length-1
    if (i==smoother){
      gamma <- gamma.vect[gamma.start:gamma.end]
    }
    gamma.start <- gamma.end+1
  }
  
  #*B) Update vectors----
  ## Calculate predictor vector eta for the parameter
  eta <- eval(parse(text=paste("family$", parameter, ".linkfun(", parameter, ")", sep="")), envir=environment())
  
  ## Calculate score and weights (for Fisher-scoring algorithm)
  dr <- eval(parse(text=paste("family$", parameter, ".dr(eta)", sep="")), envir=environment())  # dparameter/ deta
  dr <- 1/dr  # deta/ dparameter = 1/ (dparameter/ deta)
  
  # get the first and second derivatives of the log-likelihood
  if (parameter=="mu"){
    fv <- mu
    # first derivative of log-likelihood
    dldp.function <- family$dldm
    formals(dldp.function, envir=new.env()) <- alist(mu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldm2
    formals(d2ldp2.function, envir=new.env()) <- alist(mu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="sigma"){
    fv <- sigma
    # first derivative of log-likelihood
    dldp.function <- family$dldd
    formals(dldp.function, envir=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldd2
    formals(d2ldp2.function, envir=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="nu"){
    fv <- nu
    # first derivative of log-likelihood
    dldp.function <- family$dldv
    formals(dldp.function, envir=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldv2
    formals(d2ldp2.function, envir=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="tau"){
    fv <- tau
    # first derivative of log-likelihood
    dldp.function <- family$dldt
    formals(dldp.function, envir=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldt2
    formals(d2ldp2.function, envir=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  
  dldp <- dldp.function(fv)  # first derivative of log-likelihood with respect to parameter
  d2ldp2 <- d2ldp2.function(fv)  # expected  second derivative of log-Likelihood with respect to parameter
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2/(dr*dr)) # weights for Fisher-scoring algorithm =-(d2l/d2p)(dparameter/deta)^2
  # we need to stop the weights to go to Infty
  wt <- ifelse(wt>1e+10,1e+10,wt)
  wt <- ifelse(wt<1e-10,1e-10,wt)
  
  ## Update working variable vector wv
  wv <- eta+dldp/(dr*wt)
  if (family$type=="Mixed"){
    wv <-ifelse(is.nan(wv),0,wv)
  }
  
  #*C) Calculate matrix & vectors ----
  ## Calculate partial residuals
  partial.residuals <- wv - X.mat %*% beta.vect - base::rowSums(as.matrix(s[,-smoother]))
  fitted.partial.residuals <- Z.mat %*% gamma
  
  ## Calculate inner product to return to the client
  inner.product <- sum(wt*(partial.residuals - fitted.partial.residuals)^2)
  
  
  #**************************************************************************
  # III) Output ----
  #**************************************************************************
  
  return(list(nobs=nobs, inner.product=inner.product))
  
} 
# AGGREGATE FUNCTION
# gamlssDS4

