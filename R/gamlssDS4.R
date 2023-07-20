#'
#' @title gamlssDS4 called by ds.gamlss
#' @description This is the fourth serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that uses the model structure and starting
#' parameter vectors constructed by gamlssDS1 to iteratively obtain the PWLSE for gamma.
#' For more details please see the extensive header of ds.gamlss and also the
#' gamlss function in native R gamlss package.
#' @param parameter a string specifing for which of the model parameters "mu", "sigma", "nu"
#' or "tau" the model fitting should be performed
#' @param smoother an integer indicating the number of the smoother that should be fitted
#' @param formula a formula object, with the response on the left of an ~ operator, 
#' and the terms, separated by + operators, on the right. Nonparametric smoothing
#' terms are indicated by pb() for penalised beta splines, cs for smoothing splines, 
#' lo for loess smooth terms and random or ra for random terms, 
#' e.g. y~cs(x,df=5)+x1+x2*x3. 
#' @param sigma.formula a formula object for fitting a model to the sigma parameter,
#' as in the formula above, e.g. sigma.formula=~cs(x,df=5).
#' @param nu.formula a formula object for fitting a model to the nu parameter, 
#' e.g. nu.formula=~x
#' @param tau.formula a formula object for fitting a model to the tau parameter, 
#' e.g. tau.formula=~cs(x,df=2)
#' @param family a gamlss.family object, which is used to define the distribution 
#' and the link functions of the various parameters. The distribution families 
#' supported by gamlss() can be found in gamlss.family. Functions such as BI() 
#' (binomial) produce a family object. Also can be given without the parentheses
#' i.e. BI. Family functions can take arguments, as in BI(mu.link=probit).
#' @param data an optional character string specifying a data.frame object holding 
#' the data to be analysed under the specified model 
#' @param mu.fix logical, indicate whether the mu parameter should be kept fixed
#' in the fitting processes.
#' @param sigma.fix logical, indicate whether the sigma parameter should be kept
#' fixed in the fitting processes.
#' @param nu.fix logical, indicate whether the nu parameter should be kept fixed 
#' in the fitting processes.
#' @param tau.fix logical, indicate whether the tau parameter should be kept fixed
#' in the fitting processes.
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
#' @param pb.xl a numeric vector created by the clientside function specifying the left
#' boundary for the knots for the pb-smoother.
#' @param pb.xr a numeric vector created by the clientside function specifying the right
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
#' @import gamlss
#' @import gamlss.dist
#' @export
#'

gamlssDS4 <- function(parameter = parameter, smoother = smoother, formula = formula, 
                      sigma.formula = sigma.formula, nu.formula = nu.formula, 
                      tau.formula = tau.formula, family = family, data = data, 
                      mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
                      nu.beta.vect = nu.beta.vect, tau.beta.vect = tau.beta.vect,
                      mu.gamma.vect = mu.gamma.vect, sigma.gamma.vect = sigma.gamma.vect,
                      nu.gamma.vect = nu.gamma.vect, tau.gamma.vect = tau.gamma.vect,
                      pb.xl = pb.xl, pb.xr = pb.xr,
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
    data <- eval(parse(text=dataname), env=parent.frame())
  }
  Ntotal <- dim(data)[1]
  
  ## Reconvert the special symbols to create the appropriate formula, gamlss.family objects, beta & gamma vectors
  formulatext <- gsub("left_parenthesis", "(", formula, fixed = TRUE)
  formulatext <- gsub("right_parenthesis", ")", formulatext, fixed = TRUE)
  formulatext <- gsub("tilde_symbol", "~", formulatext, fixed = TRUE)
  formulatext <- gsub("equal_symbol", "=", formulatext, fixed = TRUE)
  formulatext <- gsub("comma_symbol", ",", formulatext, fixed = TRUE)
  formulatext <- gsub("asterisk_symbol", "*", formulatext, fixed = TRUE)
  formula <- stats::as.formula(formulatext)
  formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  sigma.formulatext <- gsub("left_parenthesis", "(", sigma.formula, fixed = TRUE)
  sigma.formulatext <- gsub("right_parenthesis", ")", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("tilde_symbol", "~", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("equal_symbol", "=", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("comma_symbol", ",", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("asterisk_symbol", "*", sigma.formulatext, fixed = TRUE)
  sigma.formula <- stats::as.formula(sigma.formulatext)
  sigma.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  nu.formulatext <- gsub("left_parenthesis", "(", nu.formula, fixed = TRUE)
  nu.formulatext <- gsub("right_parenthesis", ")", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("tilde_symbol", "~", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("equal_symbol", "=", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("comma_symbol", ",", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("asterisk_symbol", "*", nu.formulatext, fixed = TRUE)
  nu.formula <- stats::as.formula(nu.formulatext)
  nu.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(nu.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  tau.formulatext <- gsub("left_parenthesis", "(", tau.formula, fixed = TRUE)
  tau.formulatext <- gsub("right_parenthesis", ")", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("tilde_symbol", "~", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("equal_symbol", "=", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("comma_symbol", ",", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("asterisk_symbol", "*", tau.formulatext, fixed = TRUE)
  tau.formula <- stats::as.formula(tau.formulatext)
  tau.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(tau.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family), env=environment()))
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  # Convert parameter vectors from transmittable (character) format to numeric 
  mu.beta.vect <- as.numeric(unlist(strsplit(mu.beta.vect, split=",")))
  sigma.beta.vect <- as.numeric(unlist(strsplit(sigma.beta.vect, split=",")))
  nu.beta.vect <- as.numeric(unlist(strsplit(nu.beta.vect, split=",")))
  tau.beta.vect <- as.numeric(unlist(strsplit(tau.beta.vect, split=",")))
  mu.gamma.vect <- as.numeric(unlist(strsplit(mu.gamma.vect, split=",")))
  sigma.gamma.vect <- as.numeric(unlist(strsplit(sigma.gamma.vect, split=",")))
  nu.gamma.vect <- as.numeric(unlist(strsplit(nu.gamma.vect, split=",")))
  tau.gamma.vect <- as.numeric(unlist(strsplit(tau.gamma.vect, split=",")))
  # get the beta and gammma vectors for the respective parameter
  beta.vect <- eval(parse(text=paste(parameter, ".beta.vect", sep="")), env=environment())
  gamma.vect <- eval(parse(text=paste(parameter, ".gamma.vect", sep="")), env=environment())
  
  # Convert knot boundaries from transmittable (character) format to numeric
  pb.xl <- as.numeric(unlist(strsplit(pb.xl, split=",")))
  pb.xr <- as.numeric(unlist(strsplit(pb.xr, split=",")))
  
  #**************************************************************************
  # II) Calculate matrix & vector to return to client ----  
  #**************************************************************************
  
  #*A) Fit the model ----
  # Now fit model specified in formula:
  # to increase computational speed the number of inner and backfitting iterations are set to 1
  # suppressWarnings to avoid the warning that the algorithm has not yet converged
  mod.gamlss.ds <- base::suppressWarnings(gamlss::gamlss(formula=formula2use, sigma.formula=sigma.formula2use, 
                                                        nu.formula=nu.formula2use, tau.formula=tau.formula2use,
                                                        family=family, data=data, method=RS(), mu.fix=mu.fix, 
                                                        sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                                                        control=gamlss.control(c.crit=c1[1], n.cyc=1, 
                                                                               mu.step=c1[3], sigma.step=c1[4], 
                                                                               nu.step=c1[5], tau.step=c1[6],
                                                                               gd.tol=c1[7], trace=FALSE),
                                                        i.control=glim.control(cc=c2[1], cyc=1, 
                                                                               bf.cyc=1, bf.tol=c2[4])))
  
  ## get design matrix for the parameter
  X.mat <- as.matrix(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".x", sep="")), env=environment()))
  y <- as.vector(mod.gamlss.ds$y)
  
  ## get design matrix for the smoothers
  # identify the variables for pb-smoothers (for all parameters)
  mu.coef.names <- names(mod.gamlss.ds$mu.coefficients)
  sigma.coef.names <- names(mod.gamlss.ds$sigma.coefficients)
  nu.coef.names <- names(mod.gamlss.ds$nu.coefficients)
  tau.coef.names <- names(mod.gamlss.ds$tau.coefficients)
  smoother.names <- c(mu.coef.names, sigma.coef.names, nu.coef.names, tau.coef.names)
  pb.names <- smoother.names[grep(pattern="pb(", x=smoother.names, fixed=TRUE)]
  pb.names <- unique(pb.names)
  pb.names <- gsub(pattern="pb(", replacement="", pb.names, fixed=TRUE)
  pb.names <- gsub(pattern=")", replacement="", pb.names, fixed=TRUE)
  # only extract the variables that are relevant for the current parameter
  coef.names <- eval(parse(text=paste(parameter, ".coef.names", sep="")), env=environment())
  pb.names.parameter <- coef.names[grep(pattern="pb(", x=coef.names, fixed=TRUE)]
  pb.names.parameter <- gsub(pattern="pb(", replacement="", pb.names.parameter, fixed=TRUE)
  pb.names.parameter <- gsub(pattern=")", replacement="", pb.names.parameter, fixed=TRUE)
  pb.xl.parameter <- pb.xl[which(pb.names %in% pb.names.parameter)]
  pb.xr.parameter <- pb.xr[which(pb.names %in% pb.names.parameter)]
  # create design matrices for current & previous smoother if possible
  if(length(pb.names.parameter[smoother])>0){
    name <- eval(parse(text=paste("pb.names.parameter[", smoother, "]", sep="")), env=environment())
    x <-  eval(parse(text=name), env=parent.frame())
    Z.mat <- bbase(x=x, xl=pb.xl.parameter[smoother], xr=pb.xr.parameter[smoother])
    # get design matrix for previous smoother if possible
    if (smoother>1){
      name <- eval(parse(text=paste("pb.names.parameter[", smoother, "]", sep="")), env=environment())
      x <-  eval(parse(text=name), env=parent.frame())
      Z.mat.old <- bbase(x=x, xl=pb.xl.parameter[smoother-1], xr=pb.xr.parameter[smoother-1])
    }
  }
  
  ## Get fitted values for all distribution parameters
  # necessary for all parameters to calculate deviance
  if("mu" %in% names(family$parameters)){
    mu <- base::get("mu", env=parent.frame())
  }
  if("sigma" %in% names(family$parameters)){
    sigma <- base::get("sigma", env=parent.frame())
  }
  if("nu" %in% names(family$parameters)){
    nu <- base::get("nu", env=parent.frame())
  }
  if("tau" %in% names(family$parameters)){
    tau <- base::get("tau", env=parent.frame())
  }
  
  ## Calculate smoothing fitted value matrix s
  gamma.start <- 1
  coefSmo <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), env=environment())
  s.old <- base::get(paste(parameter, ".s", sep=""), env=parent.frame())
  # get the gamma vectors for the respective parameter & multiply them with the matrices
  for (i in 1:length(coefSmo)){
    gamma.length <- dim(coefSmo[[i]]$coef)[1]
    gamma.end <- gamma.start+gamma.length-1
    gamma <- gamma.vect[gamma.start:gamma.end]
    # calculate new smoothing fitted values for previous smoother
    if (i == (smoother-1)){
      s.update <- as.vector(Z.mat.old %*% gamma)
    }
    gamma.start <- gamma.end+1
  }
  s <- s.old
  # update smoothing fitted value matrix for previous smoother
  if (smoother>1){
    s[,(smoother-1)] <- s.update
    base::assign(paste(parameter, ".s", sep=""), s, env=parent.frame())
  }
  
  #*B) Update vectors----
  ## Calculate predictor vector eta for the parameter
  eta <- eval(parse(text=paste("family$", parameter, ".linkfun(", parameter, ")", sep="")), env=environment())
  
  ## Calculate score and weights (for Fisher-scoring algorithm)
  dr <- eval(parse(text=paste("family$", parameter, ".dr(eta)", sep="")), env=environment())  # dparameter/ deta
  dr <- 1/dr  # deta/ dparameter = 1/ (dparameter/ deta)
  
  # get the first and second derivatives of the log-likelihood
  if (parameter=="mu"){
    fv <- mu
    # first derivative of log-likelihood
    dldp.function <- family$dldm
    formals(dldp.function, env=new.env()) <- alist(mu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldm2
    formals(d2ldp2.function, env=new.env()) <- alist(mu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="sigma"){
    fv <- sigma
    # first derivative of log-likelihood
    dldp.function <- family$dldd
    formals(dldp.function, env=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldd2
    formals(d2ldp2.function, env=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="nu"){
    fv <- nu
    # first derivative of log-likelihood
    dldp.function <- family$dldv
    formals(dldp.function, env=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldv2
    formals(d2ldp2.function, env=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  if (parameter=="tau"){
    fv <- tau
    # first derivative of log-likelihood
    dldp.function <- family$dldt
    formals(dldp.function, env=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldt2
    formals(d2ldp2.function, env=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
  }
  
  dldp <- dldp.function(fv)  # first derivative of log-likelihood with respect to parameter
  d2ldp2 <- d2ldp2.function(fv)  # expected  second derivative of log-Likelihood with respect to parameter
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2/(dr*dr)) # weights for Fisher-scoring algorithm =-(d2l/d2p)(dparameter/deta)^2
  # we need to stop the weights to go to Infty
  wt <- ifelse(wt>1e+10,1e+10,wt)
  wt <- ifelse(wt<1e-10,1e-10,wt)
  
  ## stopping criterion for backfitting
  # the sums are necessary to calculate deltaf on the server which is needed to determine
  # the stoppping criterion for backfitting
  if (smoother==1){
    sumofsquares <- 0
    sumofweights <- 1
  } else {
    sumofsquares <- sum((s[,(smoother-1)] - s.old[,(smoother-1)])^2*wt)
    sumofweights <- sum(wt)
  }
  
  ## Update working variable vector wv
  wv <- eta+dldp/(dr*wt)
  if (family$type=="Mixed"){
    wv <-ifelse(is.nan(wv),0,wv)
  }
  
  #*C) Calculate matrix & vectors ----
  
  ## Calculate partial residuals
  partial.residuals <- wv - X.mat %*% beta.vect - base::rowSums(as.matrix(s[,-smoother]))
  
  ## Calculate matrix and vector to return to the client
  vector <- t(Z.mat) %*% (wt*partial.residuals)
  matrix <- t(Z.mat) %*% diag(wt) %*% Z.mat
  # remove the dimnames attributes
  attr(vector, "dimnames") <- NULL
  attr(matrix, "dimnames") <- NULL
  
  
  #**************************************************************************
  # III) Output ----
  #**************************************************************************
  
  return(list(matrix=matrix, vector=vector, sumofsquares=sumofsquares,
              sumofweights=sumofweights))
  
} 
# AGGREGATE FUNCTION
# gamlssDS4

