#'
#' @title gamlssDS6 called by ds.gamlss
#' @description This is the sixth server-side aggregate function called by \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @details It is an aggregation function that updates the distribution parameters and returns
#' the current deviance. The deviance can then be used on the client side to check whether the
#' inner iteration converged. This function is not intended for direct use by the user. For more details 
#' please see the extensive header of \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @param parameter A string specifing for which of the distribution parameters \code{c('mu', 'sigma', 'nu', 'tau')}
#' the model fitting should be performed.
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
#' @param autostep Logical, indicating whether the steps should be halved automatically 
#' if the new global deviance is greater than the old one. The default is \code{autostep=TRUE}.
#' @param inner.iteration.count An integer with the iteration count for the inner iteration.
#' @param autostep.count An integer with the iteration count for the autostep iteration.
#' @return A list with the following elements.
#' \describe{
#'  \item{\code{dv}}{Numeric value for the new deviance on the server.}
#'  \item{\code{pen}}{Numeric value with the sum of the quadratic penalties for the \code{parameter}.}
#'  \item{\code{errorMessage3}}{String for the validity of the fitted values in the inner iteration. \code{errorMessage3='Fitted values in the inner iteration out of range.'} indicates a
#'  problem, whereas \code{errorMessage3=NULL} means no problem with the fitted values.}
#' }
#' @author Annika Swenne
#' @import gamlss.dist
#' @export
#'

gamlssDS6 <- function(parameter = parameter, family = family, data = data, 
                      mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
                      nu.beta.vect = nu.beta.vect, tau.beta.vect = tau.beta.vect,
                      mu.gamma.vect = mu.gamma.vect, sigma.gamma.vect = sigma.gamma.vect,
                      nu.gamma.vect = nu.gamma.vect, tau.gamma.vect = tau.gamma.vect,
                      smoother.names = smoother.names,
                      smoother.xl = smoother.xl, smoother.xr = smoother.xr,
                      control = control, i.control = i.control, autostep = autostep,
                      inner.iteration.count = inner.iteration.count, 
                      autostep.count = autostep.count){
  
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
  
  dev.function <- family$G.dev.incr
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  mu.step <- c1[3]
  sigma.step <- c1[4]
  nu.step <- c1[5]
  tau.step <- c1[6]
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
  # II) Update distribution parameter vector ----  
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
  if (length(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), envir=environment()))>0){
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
      base::assign(paste("Z", i, ".mat", sep=""), basismatrix, env=environment())
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
  
  ## calculate smoothing fitted value matrix s
  # get the gamma vectors for the respective parameter & multiply them with the matrices
  gamma.start <- 1
  coefSmo <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), envir=environment())
  if (!is.null(coefSmo)){
    s.old <- base::get(paste("temp_", parameter, ".s", sep=""), env=parent.frame())
    s <- NULL
    for (i in 1:length(coefSmo)){
      gamma.length <- dim(coefSmo[[i]]$coef)[1]
      gamma.end <- gamma.start+gamma.length-1
      gamma <- gamma.vect[gamma.start:gamma.end]
      Z.mat <- eval(parse(text=paste("Z", i, ".mat", sep="")), envir=environment())
      s <- cbind(s, Z.mat %*% gamma)
      gamma.start <- gamma.end+1
    }
  } else{
    s <- rep(0, times=Ntotal)
  }
  s <- as.matrix(s)
  
  #*B) Update distribution parameter vector----
  # Calculate predictor vector eta for the parameter
  eta.old <- eval(parse(text=paste("family$", parameter, ".linkfun(", parameter, ")", sep="")), envir=environment())
  eta <- as.vector(X.mat %*% beta.vect + base::rowSums(as.matrix(s)))
  
  ## Weighting of the estimates
  # weight the new smoothing fitted value matrix and fitted eta with the old fitted values to avoid overjumping
  
  # fixed step size (method 1 as described in Stasinopolous et al. 2020, p.66)
  if (inner.iteration.count>1 & (autostep==FALSE | autostep.count==1)){
    # no weighting for the first inner iteration (old estimates not reasonable)
    step <- eval(parse(text=paste(parameter, ".step", sep="")), envir=environment())
    eta <- step*eta+(1-step)*eta.old
    if (!is.null(coefSmo)){
      s <- step*s+(1-step)*s.old
    }
  }
  
  # automatic halving of the step size (method 2 as described in Stasinopolous et al. 2020, p.66f)
  if (autostep==TRUE & autostep.count>0){
    # different calculation to avoid saving of the original lpold 
    # (does not work for 1st autostep iteration - division by zero)
    if (autostep.count==1){
      eta <- (eta.old + eta)/2
      if (!is.null(coefSmo)){
        s <- (s.old + s)/2
      }
    } else {
      eta <- (eta.old*(2**autostep.count - 1) - eta)/(2**autostep.count-2)
      if (!is.null(coefSmo)){
        s <- (s.old*(2**autostep.count - 1) - s)/(2**autostep.count-2)
      }
    }
  } 
  
  ## Save the new estimates
  fv <- eval(parse(text=paste("family$", parameter, ".linkinv(eta)", sep="")), envir=environment())
  
  if (autostep==FALSE | autostep.count>0){
    # Save the smoothing fitted values
    if(!is.null(coefSmo)){
      base::assign(paste("temp_", parameter, ".s", sep=""), s, env=parent.frame())
    }
    # Save the distribution parameter vector
    base::assign(paste("temp_", parameter, sep=""), fv, env=parent.frame())
  }

  ## Calculate derivatives & deviance
  dr <- eval(parse(text=paste("family$", parameter, ".dr(eta)", sep="")), envir=environment())  # dparameter/ deta
  dr <- 1/dr  # deta/ dparameter = 1/ (dparameter/ deta)
  
  if (parameter=="mu"){
    formals(dev.function, envir=new.env()) <- alist(mu = fv)
    # first derivative of log-likelihood
    dldp.function <- family$dldm
    formals(dldp.function, envir=new.env()) <- alist(mu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldm2
    formals(d2ldp2.function, envir=new.env()) <- alist(mu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="sigma"){
    formals(dev.function, envir=new.env()) <- alist(sigma = fv)
    # first derivative of log-likelihood
    dldp.function <- family$dldd
    formals(dldp.function, envir=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldd2
    formals(d2ldp2.function, envir=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="nu"){
    formals(dev.function, envir=new.env()) <- alist(nu = fv)
    # first derivative of log-likelihood
    dldp.function <- family$dldv
    formals(dldp.function, envir=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldv2
    formals(d2ldp2.function, envir=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="tau"){
    formals(dev.function, envir=new.env()) <- alist(tau = fv)
    # first derivative of log-likelihood
    dldp.function <- family$dldt
    formals(dldp.function, envir=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldt2
    formals(d2ldp2.function, envir=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  di <- dev.function(fv)  # deviance increment
  dv <- sum(di)  # the global deviance on the server
  dldp <- dldp.function(fv)  # first derivative of log-likelihood with respect to parameter
  d2ldp2 <- d2ldp2.function(fv)  # expected  second derivative of log-Likelihood with respect to parameter
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2/(dr*dr)) # weights for Fisher-scoring algorithm =-(d2l/d2p)(dparameter/deta)^2
  # we need to stop the weights to go to Infty
  wt <- ifelse(wt>1e+10,1e+10,wt)
  wt <- ifelse(wt<1e-10,1e-10,wt)
  wv <- eta+dldp/(dr*wt)
  if (family$type=="Mixed"){
    wv <-ifelse(is.nan(wv),0,wv)
  }
  
  ## Calculate penalty
  pen <- sum(eta*wt*(wv-eta))
  
  ## Check whether distribution parameter valid
  errorMessage <- NULL
  valid <- eval(parse(text=paste("family$", parameter, ".valid(fv)", sep="")), envir=environment())
  if (is.na(!valid)){
    errorMessage <- "Fitted values in the inner iteration out of range."
  } 
  
  
  #**************************************************************************
  # III) Output ----
  #**************************************************************************
  
  return(list(dv=dv, pen=pen, errorMessage3=errorMessage))
  
} 
# AGGREGATE FUNCTION
# gamlssDS6

