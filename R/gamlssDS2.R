#'
#' @title gamlssDS2 called by ds.gamlss
#' @description This is the second serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that that uses the model structure and starting
#' beta vectors constructed by gamlssDS1 to iteratively fit the gamlss model that has been
#' specified. The function gamlssDS2 also carries out a series of disclosure checks and if
#' the arguments or data fail any of those tests, model construction is blocked and an 
#' appropriate serverside error message is created and returned to ds.gamlss on the 
#' clientside. For more details please see the extensive header of ds.gamlss and also the
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
#' @param outer.iteration.count a numeric value indicating the number of the outer
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

gamlssDS2 <- function(parameter = parameter, formula = formula, sigma.formula = sigma.formula,
                      nu.formula = nu.formula, tau.formula = tau.formula, family = family, 
                      data=data, mu.fix = mu.fix, sigma.fix = sigma.fix, nu.fix = nu.fix, 
                      tau.fix = tau.fix, mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
                      nu.beta.vect = nu.beta.vect, tau.beta.vect = tau.beta.vect, control = control, 
                      i.control = i.control, outer.iteration.count = outer.iteration.count){
  
  #**************************************************************************
  # I) Preparation ---- 
  # Reconvert the transfer strings into required variable types
  # Extract the varnames
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
    data <- eval(parse(text=dataname), envir = parent.frame())
  }
  
  ## Reconvert the special symbols to create the appropriate formula, gamlss.family objects, beta & gamma vectors
  formulatext <- gsub("left_parenthesis", "(", formula, fixed = TRUE)
  formulatext <- gsub("right_parenthesis", ")", formulatext, fixed = TRUE)
  formulatext <- gsub("tilde_symbol", "~", formulatext, fixed = TRUE)
  formulatext <- gsub("equal_symbol", "=", formulatext, fixed = TRUE)
  formulatext <- gsub("comma_symbol", ",", formulatext, fixed = TRUE)
  formula <- stats::as.formula(formulatext)
  formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(formula))), env = parent.frame()) # here we need the formula as a 'call' object
  
  sigma.formulatext <- gsub("left_parenthesis", "(", sigma.formula, fixed = TRUE)
  sigma.formulatext <- gsub("right_parenthesis", ")", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("tilde_symbol", "~", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("equal_symbol", "=", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("comma_symbol", ",", sigma.formulatext, fixed = TRUE)
  sigma.formula <- stats::as.formula(sigma.formulatext)
  sigma.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula))), env = parent.frame()) # here we need the formula as a 'call' object
  
  nu.formulatext <- gsub("left_parenthesis", "(", nu.formula, fixed = TRUE)
  nu.formulatext <- gsub("right_parenthesis", ")", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("tilde_symbol", "~", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("equal_symbol", "=", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("comma_symbol", ",", nu.formulatext, fixed = TRUE)
  nu.formula <- stats::as.formula(nu.formulatext)
  nu.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(nu.formula))), env = parent.frame()) # here we need the formula as a 'call' object
  
  tau.formulatext <- gsub("left_parenthesis", "(", tau.formula, fixed = TRUE)
  tau.formulatext <- gsub("right_parenthesis", ")", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("tilde_symbol", "~", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("equal_symbol", "=", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("comma_symbol", ",", tau.formulatext, fixed = TRUE)
  tau.formula <- stats::as.formula(tau.formulatext)
  tau.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(tau.formula))), env = parent.frame()) # here we need the formula as a 'call' object
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family)))
  
  dev.function <- family$G.dev.incr
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  ## Get the variable names
  # Rewrite formulas extracting variables nested in structures like data frame or list
  # (e.g. D$A~D$B will be re-written A~B)
  # Note final product is a list of the variables in the model (yvector and covariates)
  # it is NOT a list of model terms - these are derived later
  
  # Convert formula strings into separate variable names split by |
  formulas <- paste(formulatext, sigma.formulatext, nu.formulatext, tau.formulatext, sep="|")
  formulas <- gsub("pb(", "", formulas, fixed=TRUE) 
  formulas <- gsub(")", "", formulas, fixed=TRUE) 
  formulas <- gsub(" ", "", formulas, fixed=TRUE)
  formulas <- gsub("~", "|", formulas, fixed=TRUE)
  formulas <- gsub("+", "|", formulas, fixed=TRUE)
  formulas <- gsub("*", "|", formulas, fixed=TRUE)
  formulas <- gsub("||", "|", formulas, fixed=TRUE)
  
  # Remember model.variables and then varnames include both yvect and linear predictor components 
  model.variables <- unlist(strsplit(formulas, split="|", fixed=TRUE))

  varnames <- c()
  for(i in 1:length(model.variables)){
    elt <- unlist(strsplit(model.variables[i], split="$", fixed=TRUE))
    if(length(elt) > 1){
      assign(elt[length(elt)], eval(parse(text=model.variables[i]), envir = parent.frame()), envir = parent.frame())
      varnames <- append(varnames, elt[length(elt)])
    }else{
      varnames <- append(varnames, elt)
    }
  }
  varnames <- unique(varnames)
  
  if(!is.null(dataname)){
    for(v in 1:length(varnames)){
      varnames[v] <- paste0(dataname,"$",varnames[v])
      test.string.0 <- paste0(dataname,"$","0")
      test.string.1 <- paste0(dataname,"$","1")
      if(varnames[v]==test.string.0) varnames[v] <- "0"
      if(varnames[v]==test.string.1) varnames[v] <- "1"
    }
    cbindraw.text <- paste0("cbind(", paste(varnames, collapse=","), ")")
  }else{
    cbindraw.text <- paste0("cbind(", paste(varnames, collapse=","), ")")
  }
  
  #**************************************************************************
  # II) Identify missings----  
  #**************************************************************************
  
  # Identify and use variable names to count missings
  all.data <- eval(parse(text=cbindraw.text), envir = parent.frame())
  
  Ntotal <- dim(all.data)[1]
  
  nomiss.any <- stats::complete.cases(all.data)
  nomiss.any.data <- all.data[nomiss.any,]
  N.nomiss.any <- dim(nomiss.any.data)[1]
  
  Nvalid <- N.nomiss.any
  Nmissing <- Ntotal-Nvalid
  
  #**************************************************************************
  # III) Calculate matrix & vector to return to client ----  
  #**************************************************************************
  
  #*i) Fit the model ----
  # Now fit model specified in formula:
  # to increase computational speed the number of inner and backfitting iterations are set to 1
  mod.gamlss.ds <- gamlss::gamlss(formula=formula2use, sigma.formula=sigma.formula2use, 
                                  nu.formula=nu.formula2use, tau.formula=tau.formula2use,
                                  family=family, data=data, method=RS(),
                                  mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix,
                                  tau.fix=tau.fix,
                                  control = gamlss.control(c.crit=c1[1], n.cyc=1, 
                                                           mu.step=c1[3], sigma.step=c1[4], 
                                                           nu.step=c1[5], tau.step=c1[6],
                                                           gd.tol=c1[7]),
                                  i.control = glim.control(cc=c2[1], cyc=1, 
                                                           bf.cyc=1, bf.tol=c2[4]))
  
  # get design matrix for the parameter
  X.mat.orig <- as.matrix(eval(parse(text=paste("mod.gamlss.ds$", parameter, ".x", sep=""))))
  y <- as.vector(mod.gamlss.ds$y)
  
  #*ii) Update vectors----
  ## Calculate predictor vector eta for the parameter
  if("mu" %in% names(family$parameters)){
    mu <- base::get("mu", envir = parent.frame())
  }
  if("sigma" %in% names(family$parameters)){
    sigma <- base::get("sigma", envir = parent.frame())
  }
  if("nu" %in% names(family$parameters)){
    nu <- base::get("nu", envir = parent.frame())
  }
  if("tau" %in% names(family$parameters)){
    tau <- base::get("tau", envir = parent.frame())
  }
  eta <- eval(parse(text=paste("family$", parameter, ".linkfun(", parameter, ")", sep="")))
  
  ## Calculate score and weights (for Fisher-scoring algorithm)
  dr <- eval(parse(text=paste("family$", parameter, ".dr(eta)", sep="")))  # dparameter/ deta
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
  
  ## Update working variable vector wv
  wv <- eta+dldp/(dr*wt)
  if (family$type=="Mixed"){
    wv <-ifelse(is.nan(wv),0,wv)
  }
  
  ## Calculate deviance
  di <- dev.function(fv)  # deviance increment
  dv <- sum(di)  # the global deviance on the server
  
  
  #**************************************************************************
  # IV) Save working variable on the server  ----  
  #**************************************************************************
  
  # We need to save the working variable on the server for the next iteration
  # Since the working variable could be disclosive we keep it on the server and do not
  # reveal it to the client
  
  
  
  
  
} 
# AGGREGATE FUNCTION
# gamlssDS2


