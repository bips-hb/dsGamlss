#'
#' @title gamlssDS2 called by ds.gamlss
#' @description This is the second server-side aggregate function called by \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @details It is an aggregation function that uses the model structure and starting
#' parameter vectors constructed by \code{gamlssDS1} to iteratively obtain the weighted least squares estimator for beta.
#' The function \code{gamlssDS2} also carries out a series of disclosure checks and if
#' the arguments or data fail any of those tests, model construction is blocked and an 
#' appropriate server-side error message is created and returned to \code{\link[dsGamlssClient]{ds.gamlss}} on the 
#' client-side. This function is not intended for direct use by the user. For more details 
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
#'  \item{\code{matrix}}{Numeric matrix that can be aggregated on the client-side to obtain the updated weighted least squares estimator for the distribution parameter.}
#'  \item{\code{vector}}{Numeric vector that can be aggregated on the client-side to obtain the updated weighted least squares estimator for the distribution parameter.}
#'  \item{\code{dv}}{Numeric value for the new deviance on the server.}
#'  \item{\code{disclosure.risk}}{Numeric value, either \code{0} or \code{1}, whereby \code{1} indicates a disclosure risk.}
#'  \item{\code{errorMessage2}}{String for the disclosure risk. \code{errorMessage='No errors'} means that no disclosure risk was identified.}
#' }
#' @author Annika Swenne
#' @import gamlss.dist
#' @export
#'

gamlssDS2 <- function(parameter = parameter, family = family, 
                      data = data, mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
                      nu.beta.vect = nu.beta.vect, tau.beta.vect = tau.beta.vect,
                      control = control, i.control = i.control){
  
  #**************************************************************************
  # I) Preparation ---- 
  # Reconvert the transfer strings into required variable types
  # Get function to calculate deviance increment
  #**************************************************************************
  
  ## Capture the nfilter settings
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  
  ## Get the value of the 'data' parameter provided as character on the client side
  dataname <- data
  if(is.null(dataname)){
    data <- NULL 
  }else{
    data <- eval(parse(text=dataname), envir=parent.frame())
  }
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family), envir=environment()))
  
  dev.function <- family$G.dev.incr
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  # Convert parameter vectors from transmittable (character) format to numeric 
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
  
  # get the beta and gammma vectors for the respective parameter
  beta.vect <- eval(parse(text=paste(parameter, ".beta.vect", sep="")), envir=environment())
  
  #**************************************************************************
  # II) Calculate matrix & vector to return to client ----  
  #**************************************************************************
  
  #*A) Get fitted model ----
  mod.gamlss.ds <- base::get("temp_mod.gamlss.ds", env=parent.frame())
  
  ## get design matrix for the parameter
  X.mat <- as.matrix(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".x", sep="")), envir=environment()))
  y <- as.vector(mod.gamlss.ds$y)
  
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
  } else{
    s <- rep(0, times=mod.gamlss.ds$N)
  }
  s <- as.matrix(s)
  
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
    # deviance
    formals(dev.function, envir=new.env()) <- alist(mu = fv)
  }
  if (parameter=="sigma"){
    fv <- sigma
    # first derivative of log-likelihood
    dldp.function <- family$dldd
    formals(dldp.function, envir=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldd2
    formals(d2ldp2.function, envir=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # deviance
    formals(dev.function, envir=new.env()) <- alist(sigma = fv)
  }
  if (parameter=="nu"){
    fv <- nu
    # first derivative of log-likelihood
    dldp.function <- family$dldv
    formals(dldp.function, envir=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldv2
    formals(d2ldp2.function, envir=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # deviance
    formals(dev.function, envir=new.env()) <- alist(nu = fv)
  }
  if (parameter=="tau"){
    fv <- tau
    # first derivative of log-likelihood
    dldp.function <- family$dldt
    formals(dldp.function, envir=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldt2
    formals(d2ldp2.function, envir=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # deviance
    formals(dev.function, envir=new.env()) <- alist(tau = fv)
  }
  
  dldp <- dldp.function(fv)  # first derivative of log-likelihood with respect to parameter
  d2ldp2 <- d2ldp2.function(fv)  # expected  second derivative of log-Likelihood with respect to parameter
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2/(dr*dr)) # weights for Fisher-scoring algorithm =-(d2l/d2p)(dparameter/deta)^2
  # we need to stop the weights to go to Infty
  wt <- ifelse(wt>1e+10,1e+10,wt)
  wt <- ifelse(wt<1e-10,1e-10,wt)
  
  ## Calculate deviance (this is the initial deviance)
  di <- dev.function(fv)  # deviance increment
  dv <- sum(di)  # the global deviance on the server
  
  ## Update working variable vector wv
  wv <- eta+dldp/(dr*wt)
  if (family$type=="Mixed"){
    wv <-ifelse(is.nan(wv),0,wv)
  }
  
  #*C) Calculate matrix & vectors ----
  ## Calculate partial residuals
  partial.residuals <- wv - base::rowSums(s)
  
  ## Calculate matrix and vector to return to the client
  vector <- t(X.mat) %*% (wt*partial.residuals)
  matrix <- t(X.mat) %*% diag(wt) %*% X.mat
  # remove the dimnames attributes
  attr(vector, "dimnames") <- NULL
  attr(matrix, "dimnames") <- NULL
  
  
  #**************************************************************************
  # III) Backup disclosure risk----
  # If y, X or w data are invalid but user has modified client-side
  # function (ds.gamlss) to circumvent trap, model will get to this point without
  # giving a controlled shut down with a warning about invalid data.
  # So as a safety measure, we will now use the same test that is used to
  # trigger a controlled trap in the client-side function to destroy the
  # vector and matrix in the study with the problem.
  # So this will make model fail without explanation
  
  # Disclosure code from gamlssDS1
  #**************************************************************************
  
  errorMessage.combined <- NULL
  disclosure.risk <- 0
  
  mu.x <- mod.gamlss.ds$mu.x
  sigma.x <- mod.gamlss.ds$sigma.x
  nu.x <- mod.gamlss.ds$nu.x
  tau.x <- mod.gamlss.ds$tau.x
  dim.mu.x <- dim(mu.x)
  dim.sigma.x <- dim(sigma.x)
  dim.nu.x <- dim(nu.x)
  dim.tau.x <- dim(tau.x)
  
  #*A) Oversaturated model----
  # (test against nfilter.glm)
  gamlss.saturation.invalid <- 0
  num.n <- mod.gamlss.ds$N
  num.mu.p <- dim.mu.x[2]
  num.sigma.p <- dim.sigma.x[2]
  num.nu.p <- dim.nu.x[2]
  num.tau.p <- dim.tau.x[2]
  
  if(mod.gamlss.ds$df.fit > nfilter.glm*num.n){
    gamlss.saturation.invalid <- 1
    errorMessage.combined <- c(errorMessage.combined,
                               "ERROR: Model has too many parameters, there is a possible risk of disclosure - please simplify model")
  }
  
  if(!is.null(num.mu.p)){
    if(num.mu.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage.combined <- c(errorMessage.combined,
                                 "ERROR: Model for mu has too many parameters, there is a possible risk of disclosure - please simplify model")
    }
  }
  
  if(!is.null(num.sigma.p)){
    if(num.sigma.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage.combined <- c(errorMessage.combined,
                                 "ERROR: Model for sigma has too many parameters, there is a possible risk of disclosure - please simplify model")
    }
  }
  
  if(!is.null(num.nu.p)){
    if(num.nu.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage.combined <- c(errorMessage.combined,
                                 "ERROR: Model for nu has too many parameters, there is a possible risk of disclosure - please simplify model")
    }
  }
  
  if(!is.null(num.tau.p)){
    if(num.tau.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage.combined <- c(errorMessage.combined,
                                 "ERROR: Model for tau has too many parameters, there is a possible risk of disclosure - please simplify model")
    }
  }
  
  #*B) Invalid y, mu.x, sigma.x, nu.x or tau.x ----
  # If y, X or w data are invalid but user has modified client-side
  # function (ds.gamlss) to circumvent trap, model will get to this point without
  # giving a controlled shut down with a warning about invalid data.
  # So as a safety measure, we will now use the same test that is used to
  # trigger a controlled trap in the client-side function to destroy the
  # vector and matrix in the study with the problem.
  
  ## check y vector validity
  y.invalid <- 0
  
  # count number of unique non-missing values (disclosure risk only arises with two levels)
  unique.values.noNA.y <- unique(y[stats::complete.cases(y)])
  
  # if two levels, check whether either level 0 < n < nfilter.tab
  if(length(unique.values.noNA.y)==2){
    tabvar <- table(y)[table(y)>=1]  # tabvar counts n in all categories with at least one observation
    min.category <- min(tabvar)
    if(min.category < nfilter.tab){
      y.invalid <- 1
      errorMessage.combined <- c(errorMessage.combined,
                                 "ERROR: y vector is binary with one category less than filter threshold for table cell size")
    }
  }
  
  ## check validity of design matrices
  # Check no dichotomous X vectors with between 1 and filter.threshold
  # observations at either level
  
  # mu
  mu.par.invalid <- NULL
  if(!is.null(num.mu.p)){
    mu.par.invalid <- rep(0, times=num.mu.p)
    for(pj in 1:num.mu.p){
      unique.values.noNA <- unique((mu.x[,pj])[stats::complete.cases(mu.x[,pj])])
      if(length(unique.values.noNA)==2){
        tabvar <- table(mu.x[,pj])[table(mu.x[,pj])>=1]  # tabvar counts n in all categories with at least one observation
        min.category <- min(tabvar)
        if(min.category < nfilter.tab){
          mu.par.invalid[pj] <- 1
          errorMessage.combined <- c(errorMessage.combined,
                                     "ERROR: at least one column in mu.x matrix is binary with one category less than filter threshold for table cell size")
        }
      }
    }
  }
  
  # sigma
  sigma.par.invalid <- NULL
  if(!is.null(num.sigma.p)){
    sigma.par.invalid <- rep(0, times=num.sigma.p)
    for(pj in 1:num.sigma.p){
      unique.values.noNA <- unique((sigma.x[,pj])[stats::complete.cases(sigma.x[,pj])])
      if(length(unique.values.noNA)==2){
        tabvar <- table(sigma.x[,pj])[table(sigma.x[,pj])>=1]  # tabvar counts n in all categories with at least one observation
        min.category <- min(tabvar)
        if(min.category < nfilter.tab){
          sigma.par.invalid[pj] <- 1
          errorMessage.combined <- c(errorMessage.combined,
                                     "ERROR: at least one column in sigma.x matrix is binary with one category less than filter threshold for table cell size")
        }
      }
    }
  }
  
  # nu
  nu.par.invalid <- NULL
  if(!is.null(num.nu.p)){
    nu.par.invalid <- rep(0, times=num.nu.p)
    for(pj in 1:num.nu.p){
      unique.values.noNA <- unique((nu.x[,pj])[stats::complete.cases(nu.x[,pj])])
      if(length(unique.values.noNA)==2){
        tabvar <- table(nu.x[,pj])[table(nu.x[,pj])>=1]  # tabvar counts n in all categories with at least one observation
        min.category <- min(tabvar)
        if(min.category < nfilter.tab){
          nu.par.invalid[pj] <- 1
          errorMessage.combined <- c(errorMessage.combined,
                                     "ERROR: at least one column in nu.x matrix is binary with one category less than filter threshold for table cell size")
        }
      }
    }
  }
  
  # tau
  tau.par.invalid <- NULL
  if(!is.null(num.tau.p)){
    tau.par.invalid <- rep(0, times=num.tau.p)
    for(pj in 1:num.tau.p){
      unique.values.noNA <- unique((tau.x[,pj])[stats::complete.cases(tau.x[,pj])])
      if(length(unique.values.noNA)==2){
        tabvar <- table(tau.x[,pj])[table(tau.x[,pj])>=1]  # tabvar counts n in all categories with at least one observation
        min.category <- min(tabvar)
        if(min.category < nfilter.tab){
          tau.par.invalid[pj] <- 1
          errorMessage.combined <- c(errorMessage.combined,
                                     "ERROR: at least one column in tau.x matrix is binary with one category less than filter threshold for table cell size")
        }
      }
    }
  }
  
  #*C) Combine disclosure risks----
  # If there is a disclosure risk destroy the matrix and vector (that should be returned to the client)
  if(!(y.invalid>0 || sum(mu.par.invalid)>0|| sum(sigma.par.invalid)>0 || sum(nu.par.invalid)>0 ||
       sum(tau.par.invalid)>0 || gamlss.saturation.invalid>0)){
    errorMessage.combined <- "No errors"
  }else{
    disclosure.risk <- 1
    matrix <- NA
    vector <- NA
    dv <- NA
    errorMessage.combined <- c(errorMessage.combined, "MODEL FAILED: model or data invalid, matrix and vector destroyed")
    rm(list=ls(pattern="^temp_", envir=parent.frame()), envir=parent.frame())
  }
  
  
  #**************************************************************************
  # V) Output ----
  #**************************************************************************
  
  return(list(matrix=matrix, vector=vector, dv=dv, disclosure.risk=disclosure.risk,
              errorMessage2=errorMessage.combined))
  
} 
# AGGREGATE FUNCTION
# gamlssDS2