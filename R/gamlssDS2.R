#'
#' @title gamlssDS2 called by ds.gamlss
#' @description This is the second serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that uses the model structure and starting
#' parameter vectors constructed by gamlssDS1 to iteratively obtain the WLSE for beta.
#' The function gamlssDS2 also carries out a series of disclosure checks and if
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

gamlssDS2 <- function(parameter = parameter, formula = formula, sigma.formula = sigma.formula,
                      nu.formula = nu.formula, tau.formula = tau.formula, family = family, 
                      data = data, mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
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
  
  ## Capture the nfilter settings
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  
  ## Get the value of the 'data' parameter provided as character on the client side
  dataname <- data
  if(is.null(dataname)){
    data <- NULL 
  }else{
    data <- eval(parse(text=dataname), env=parent.frame())
  }
  
  ## Reconvert the special symbols to create the appropriate formula, gamlss.family objects, beta & gamma vectors
  formulatext <- gsub("left_parenthesis", "(", formula, fixed = TRUE)
  formulatext <- gsub("right_parenthesis", ")", formulatext, fixed = TRUE)
  formulatext <- gsub("tilde_symbol", "~", formulatext, fixed = TRUE)
  formulatext <- gsub("equal_symbol", "=", formulatext, fixed = TRUE)
  formulatext <- gsub("comma_symbol", ",", formulatext, fixed = TRUE)
  formula <- stats::as.formula(formulatext)
  formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  sigma.formulatext <- gsub("left_parenthesis", "(", sigma.formula, fixed = TRUE)
  sigma.formulatext <- gsub("right_parenthesis", ")", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("tilde_symbol", "~", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("equal_symbol", "=", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("comma_symbol", ",", sigma.formulatext, fixed = TRUE)
  sigma.formula <- stats::as.formula(sigma.formulatext)
  sigma.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  nu.formulatext <- gsub("left_parenthesis", "(", nu.formula, fixed = TRUE)
  nu.formulatext <- gsub("right_parenthesis", ")", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("tilde_symbol", "~", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("equal_symbol", "=", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("comma_symbol", ",", nu.formulatext, fixed = TRUE)
  nu.formula <- stats::as.formula(nu.formulatext)
  nu.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(nu.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  tau.formulatext <- gsub("left_parenthesis", "(", tau.formula, fixed = TRUE)
  tau.formulatext <- gsub("right_parenthesis", ")", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("tilde_symbol", "~", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("equal_symbol", "=", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("comma_symbol", ",", tau.formulatext, fixed = TRUE)
  tau.formula <- stats::as.formula(tau.formulatext)
  tau.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(tau.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family), env=environment()))
  
  dev.function <- family$G.dev.incr
  
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
      assign(elt[length(elt)], eval(parse(text=model.variables[i]), env=parent.frame()), env=parent.frame())
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
  all.data <- eval(parse(text=cbindraw.text), env=parent.frame())
  
  Ntotal <- dim(all.data)[1]
  
  nomiss.any <- stats::complete.cases(all.data)
  nomiss.any.data <- all.data[nomiss.any,]
  N.nomiss.any <- dim(nomiss.any.data)[1]
  
  Nvalid <- N.nomiss.any
  Nmissing <- Ntotal-Nvalid
  
  #**************************************************************************
  # III) Calculate matrix & vector to return to client ----  
  #**************************************************************************
  
  #*A) Fit the model ----
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
  # create design matrices for them
  if(length(pb.names.parameter)>0){
    for (i in 1:length(pb.names.parameter)){
      name <- eval(parse(text=paste("pb.names.parameter[", i, "]", sep="")), env=environment())
      x <- eval(parse(text=name), env=parent.frame())
      basismatrix <- bbase(x=x, xl=pb.xl.parameter[i], xr=pb.xr.parameter[i])
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
    # deviance
    formals(dev.function, env=new.env()) <- alist(mu = fv)
  }
  if (parameter=="sigma"){
    fv <- sigma
    # first derivative of log-likelihood
    dldp.function <- family$dldd
    formals(dldp.function, env=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldd2
    formals(d2ldp2.function, env=new.env()) <- alist(sigma = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # deviance
    formals(dev.function, env=new.env()) <- alist(sigma = fv)
  }
  if (parameter=="nu"){
    fv <- nu
    # first derivative of log-likelihood
    dldp.function <- family$dldv
    formals(dldp.function, env=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldv2
    formals(d2ldp2.function, env=new.env()) <- alist(nu = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # deviance
    formals(dev.function, env=new.env()) <- alist(nu = fv)
  }
  if (parameter=="tau"){
    fv <- tau
    # first derivative of log-likelihood
    dldp.function <- family$dldt
    formals(dldp.function, env=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # second derivative of log-likelihood
    d2ldp2.function <- family$d2ldt2
    formals(d2ldp2.function, env=new.env()) <- alist(tau = fv)  # replace function parameters ($y, $mu, $sigma, $nu, $tau)
    # deviance
    formals(dev.function, env=new.env()) <- alist(tau = fv)
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
  # IV) Backup disclosure risk----
  # If y, X or w data are invalid but user has modified clientside
  # function (ds.gamlss) to circumvent trap, model will get to this point without
  # giving a controlled shut down with a warning about invalid data.
  # So as a safety measure, we will now use the same test that is used to
  # trigger a controlled trap in the clientside function to destroy the
  # score.vector and information.matrix in the study with the problem.
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
  # If y, X or w data are invalid but user has modified clientside
  # function (ds.gamlss) to circumvent trap, model will get to this point without
  # giving a controlled shut down with a warning about invalid data.
  # So as a safety measure, we will now use the same test that is used to
  # trigger a controlled trap in the clientside function to destroy the
  # score.vector and information.matrix in the study with the problem.
  
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
    info.matrix <- NA
    score.vector <- NA
    disclosure.risk <- 1
    errorMessage.combined <- c(errorMessage.combined, "MODEL FAILED: model or data invalid, matrix and vector destroyed")
  }
  
  
  #**************************************************************************
  # V) Output ----
  #**************************************************************************
  
  return(list(matrix=matrix, vector=vector, dv=dv, disclosure.risk=disclosure.risk,
              errorMessage2=errorMessage.combined))
  
} 
# AGGREGATE FUNCTION
# gamlssDS2