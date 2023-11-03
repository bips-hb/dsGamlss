#'
#' @title gamlssDS6 called by ds.gamlss
#' @description This is the sixth serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that updates the distribution parameters and returns
#' the current deviance. The deviance can then be used on the client side to check whether the
#' innter iteration converged. For more details please see the extensive header of ds.gamlss and also the
#' gamlss function in native R gamlss package.
#' @param parameter a string specifing for which of the model parameters "mu", "sigma", "nu"
#' or "tau" the model fitting should be performed
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
#' @param autostep whether the steps should be halved automatically if the new global 
#' deviance is greater that the old one. 
#' @param inner.iteration.count an integer specifying the iteration count for the inner
#' iteration
#' @param autostep.count an integer specifying the iteration count for the autostep
#' iteration
#' @return a gamlss object with all components as in the native R gamlss function. 
#' Individual-level information like the components y (the response response) and 
#' residuals (the normalised quantile residuals of the model) are not disclosed to 
#' the client-side.
#' @author Annika Swenne
#' @import gamlss
#' @import gamlss.dist
#' @export
#'

gamlssDS6 <- function(parameter = parameter, formula = formula, 
                      sigma.formula = sigma.formula, nu.formula = nu.formula, 
                      tau.formula = tau.formula, family = family, data = data, 
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
  
  dev.function <- family$G.dev.incr
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  mu.step <- c1[3]
  sigma.step <- c1[4]
  nu.step <- c1[5]
  tau.step <- c1[6]
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
  smoother.xl <- as.numeric(unlist(strsplit(smoother.xl, split=",")))
  smoother.xr <- as.numeric(unlist(strsplit(smoother.xr, split=",")))
  
  #**************************************************************************
  # II) Update distribution parameter vector ----  
  #**************************************************************************
  
  #*A) Fit the model ----
  # Now fit model specified in formula:
  # to increase computational speed the number of inner and backfitting iterations are set to 1
  mod.gamlss.ds <- base::suppressWarnings(gamlss::gamlss(formula=formula2use, sigma.formula=sigma.formula2use, 
                                                         nu.formula=nu.formula2use, tau.formula=tau.formula2use,
                                                         family=family, data=data, method=RS(), mu.fix=mu.fix, 
                                                         sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                                                         control = gamlss.control(c.crit=c1[1], n.cyc=1, 
                                                                                  mu.step=c1[3], sigma.step=c1[4], 
                                                                                  nu.step=c1[5], tau.step=c1[6],
                                                                                  gd.tol=c1[7], trace=FALSE),
                                                         i.control = glim.control(cc=c2[1], cyc=1, 
                                                                                  bf.cyc=1, bf.tol=c2[4])))
  
  ## get design matrix for the parameter
  X.mat <- as.matrix(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".x", sep="")), env=environment()))
  y <- as.vector(mod.gamlss.ds$y)
  
  ## get design matrix for the smoothers
  # get the control parameters for the smoothers
  coefficients <- eval(parse(text=paste("names(mod.gamlss.ds$", parameter, ".coefficients)", sep="")), env=environment())
  smoother.coef <- coefficients[grep(pattern="pb(", x=tolower(coefficients), fixed=TRUE)]
  # only keep the arguments for the pb() function
  pb.args <- substr(smoother.coef, start=4, stop=nchar(smoother.coef)-1)
  pb.args <- strsplit(pb.args, split=",", fixed=TRUE)
  
  # create design matrices for all smoothers if possible
  if (length(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), env=environment()))>0){
    for (i in 1:length(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), env=environment()))){
      name <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", i, "]]$name", sep="")), env=environment())
      if (length(grep(pattern="pb.control", x=pb.args[[i]], fixed=TRUE))>0) {
        # control parameters specified
        pb.control <- eval(parse(text=pb.args[[i]][grep(pattern="pb.control", x=pb.args[[i]], fixed=TRUE)]))
      } else {
        # no control parameters specified - use default
        pb.control <- eval(parse(text="pb.control()"))
      }
      x <- eval(parse(text=name), env=parent.frame())
      basismatrix <- bbase(x=x, xl=smoother.xl[which(smoother.names==name)], xr=smoother.xr[which(smoother.names==name)],
                           ndx=pb.control$inter, deg=pb.control$degree)
      base::assign(paste("Z", i, ".mat", sep=""), basismatrix, env=environment())
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
  
  ## calculate smoothing fitted value matrix s
  # get the gamma vectors for the respective parameter & multiply them with the matrices
  gamma.start <- 1
  coefSmo <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo", sep="")), env=environment())
  if (!is.null(coefSmo)){
    s.old <- base::get(paste(parameter, ".s", sep=""), env=parent.frame())
    s <- NULL
    for (i in 1:length(coefSmo)){
      gamma.length <- dim(coefSmo[[i]]$coef)[1]
      gamma.end <- gamma.start+gamma.length-1
      gamma <- gamma.vect[gamma.start:gamma.end]
      Z.mat <- eval(parse(text=paste("Z", i, ".mat", sep="")), env=environment())
      s <- cbind(s, Z.mat %*% gamma)
      gamma.start <- gamma.end+1
    }
  } else{
    s <- rep(0, times=Ntotal)
  }
  s <- as.matrix(s)
  
  #*B) Update distribution parameter vector----
  # Calculate predictor vector eta for the parameter
  eta.old <- eval(parse(text=paste("family$", parameter, ".linkfun(", parameter, ")", sep="")), env=environment())
  eta <- as.vector(X.mat %*% beta.vect + base::rowSums(as.matrix(s)))
  
  ## Weighting of the estimates
  # weight the new smoothing fitted value matrix and fitted eta with the old fitted values to avoid overjumping
  
  # fixed step size (method 1 as described in Stasinopolous et al. 2020, p.66)
  if (inner.iteration.count>1 & (autostep==FALSE | autostep.count==1)){
    # no weighting for the first inner iteration (old estimates not reasonable)
    step <- eval(parse(text=paste(parameter, ".step", sep="")), env=environment())
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
  fv <- eval(parse(text=paste("family$", parameter, ".linkinv(eta)", sep="")), env=environment())
  
  if (autostep==FALSE | autostep.count>0){
    # Save the smoothing fitted values
    if(!is.null(coefSmo)){
      base::assign(paste(parameter, ".s", sep=""), s, env=parent.frame())
    }
    # Save the distribution parameter vector
    base::assign(parameter, fv, env=parent.frame())
  }

  ## Calculate deviance
  if (parameter=="mu"){
    formals(dev.function, env=new.env()) <- alist(mu = fv)
  }
  if (parameter=="sigma"){
    formals(dev.function, env=new.env()) <- alist(sigma = fv)
  }
  if (parameter=="nu"){
    formals(dev.function, env=new.env()) <- alist(nu = fv)
  }
  if (parameter=="tau"){
    formals(dev.function, env=new.env()) <- alist(tau = fv)
  }
  di <- dev.function(fv)  # deviance increment
  dv <- sum(di)  # the global deviance on the server
  
  ## Check whether distribution parameter valid
  errorMessage <- NULL
  valid <- eval(parse(text=paste("family$", parameter, ".valid(fv)", sep="")), env=environment())
  if (is.na(!valid)){
    errorMessage <- "Fitted values in the inner iteration out of range."
  } 
  
  
  #**************************************************************************
  # III) Output ----
  #**************************************************************************
  
  return(list(dv=dv, errorMessage3=errorMessage))
  
} 
# AGGREGATE FUNCTION
# gamlssDS6

