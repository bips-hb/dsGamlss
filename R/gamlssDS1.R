#'
#' @title gamlssDS1 called by ds.gamlss
#' @description This is the first server-side aggregate function called by \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @details It is an aggregation function that sets up the appropriate model structure
#' and dimensions to fit a \code{ds.gamlss} model. This function is not intended for direct
#' use by the user. For more details please see the extensive header of \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @param formula A formula string in the legal transmission format for DataSHIELD, 
#' specifying the model for the mu distribution parameter. The DataSHIELD legal transmission
#' format means that special characters, like '(' are replaced with the corresponding verbal
#' descriptions, e.g. 'left_parenthesis'.
#' @param sigma.formula A formula string in the legal transmission format for DataSHIELD, 
#' specifying the model for the sigma distribution parameter. The DataSHIELD legal transmission
#' format means that special characters, like '(' are replaced with the corresponding verbal
#' descriptions, e.g. 'left_parenthesis'.
#' @param nu.formula A formula string in the legal transmission format for DataSHIELD, 
#' specifying the model for the nu distribution parameter. The DataSHIELD legal transmission
#' format means that special characters, like '(' are replaced with the corresponding verbal
#' descriptions, e.g. 'left_parenthesis'.
#' @param tau.formula A formula string in the legal transmission format for DataSHIELD, 
#' specifying the model for the tau distribution parameter. The DataSHIELD legal transmission
#' format means that special characters, like '(' are replaced with the corresponding verbal
#' descriptions, e.g. 'left_parenthesis'.
#' @param family A family string in the legal transmission format for DataSHIELD, which
#' is used to define the distribution of the response variable. The DataSHIELD legal transmission
#' format means that special characters, like '(' are replaced with the corresponding verbal
#' descriptions, e.g. 'left_parenthesis'. Currently, only the following families are supported:
#' \code{family=c('NOleft_parenthesisright_parenthesis', 'NO2left_parenthesisright_parenthesis', 
#' 'BCCGleft_parenthesisright_parenthesis', 'BCPEleft_parenthesisright_parenthesis')}.
#' @param data A character string specifying a data.frame object holding 
#' the data to be analysed under the specified model. 
#' @param mu.fix Logical, indicating whether the mu parameter should be kept fixed
#' in the fitting processes. Default \code{mu.fix=FALSE}.
#' @param sigma.fix Logical, indicating whether the sigma parameter should be kept
#' fixed in the fitting processes. Default \code{sigma.fix=FALSE}.
#' @param nu.fix Logical, indicating whether the nu parameter should be kept fixed 
#' in the fitting processes. Default \code{nu.fix=FALSE}.
#' @param tau.fix Logical, indicating whether the tau parameter should be kept fixed
#' in the fitting processes. Default \code{tau.fix=FALSE}.
#' @param global.mean Numeric value, giving the global mean of the outcome variable,
#' which is necessary to initialize the distribution parameters for some families.
#' Otherwise \code{global.mean=NULL}.
#' @param global.sd Numeric value, giving the global sd of the outcome variable
#' which is necessary to initialize the distribution parameters for some families.
#' Otherwise \code{global.sd=NULL}.
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
#' @return A list with the following elements.
#' \describe{
#'  \item{\code{mod.gamlss.ds}}{A \code{gamlss} object with all components as in the native R \code{\link[gamlss]{gamlss}}
#'  function. Individual-level information like the components \code{y} (the response) and \code{residuals} (the normalised
#'  quantile residuals of the model) are not disclosed to the client-side.}
#'  \item{\code{G.dev}}{Numeric value for the initial deviance on the server.}
#'  \item{\code{dim.mu.x}}{Numeric vector with two elements, specifying the dimension of the design matrix for mu.}
#'  \item{\code{dim.sigma.x}}{Numeric vector with two elements, specifying the dimension of the design matrix for sigma.}
#'  \item{\code{dim.nu.x}}{Numeric vector with two elements, specifying the dimension of the design matrix for nu.}
#'  \item{\code{dim.tau.x}}{Numeric vector with two elements, specifying the dimension of the design matrix for tau.}
#'  \item{\code{smoother.names}}{String vector with the unique variable names that are used for smoothing.}
#'  \item{\code{smoother.xmin}}{Numeric vector with anononymized minima for the variables in \code{smoother.names}.}
#'  \item{\code{smoother.xmax}}{Numeric vector with anononymized maxima for the variables in \code{smoother.names}.}
#'  \item{\code{y.invalid}}{Numeric value, either \code{0} or \code{1}, whereby \code{1} indicates a disclosure risk in the response variable.}
#'  \item{\code{mu.par.invalid}}{Numeric vector with elements \code{0} or \code{1}, whereby \code{1} indicates a disclosure risk in the corresponding
#'  explanatory variable for mu.}
#'  \item{\code{sigma.par.invalid}}{Numeric vector with elements \code{0} or \code{1}, whereby \code{1} indicates a disclosure risk in the corresponding
#'  explanatory variable for sigma.}
#'  \item{\code{nu.par.invalid}}{Numeric vector with elements \code{0} or \code{1}, whereby \code{1} indicates a disclosure risk in the corresponding
#'  explanatory variable for nu.}
#'  \item{\code{tau.par.invalid}}{Numeric vector with elements \code{0} or \code{1}, whereby \code{1} indicates a disclosure risk in the corresponding
#'  explanatory variable for tau.}
#'  \item{\code{gamlss.saturation.invalid}}{Numeric value, either \code{0} or \code{1}, whereby \code{1} indicates a disclosure risk from an oversaturated model.}
#'  \item{\code{errorMessage}}{String for the disclosure risk. \code{errorMessage='Study data or applied model invalid for this source'} indicates a
#'  disclosure risk, whereas \code{errorMessage='No errors'} means that no disclosure risk was identified.}
#' }
#' @author Annika Swenne
#' @import gamlss
#' @import gamlss.dist
#' @export
#'

gamlssDS1 <- function(formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula,
                     tau.formula = tau.formula, family = family, data = data, mu.fix=mu.fix, sigma.fix = sigma.fix, 
                     nu.fix = nu.fix, tau.fix = tau.fix, global.mean = global.mean, 
                     global.sd = global.sd, control = control, i.control = i.control,
                     autostep = autostep){
  
  
  #**************************************************************************
  # I) Preparation ---- 
  # Reconvert the transfer strings into required variable types
  #**************************************************************************
  
  ## Capture the nfilter settings
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  nfilter.noise <- as.numeric(thr$nfilter.noise)
  
  errorMessage <- "No errors"
  
  ## Get the value of the 'data' parameter provided as character on the client side
  if(is.null(data)){
    dataTable <- NULL 
  }else{
    dataTable <- eval(parse(text=data), envir=parent.frame())
  }
  
  # Create an environment with elements from parent.frame() and gamlss namespace 
  combined_env <- new.env()
  # Populate with objects from gamlss namespace
  list2env(as.list(asNamespace("gamlss")), envir = combined_env)
  # Populate with objects from parent.frame()  (overwriting duplicates if any)
  list2env(as.list(parent.frame()), envir = combined_env)
  
  ## Reconvert the special symbols to create the appropriate formula & gamlss.family objects
  formulatext <- gsub("left_parenthesis", "(", formula, fixed = TRUE)
  formulatext <- gsub("right_parenthesis", ")", formulatext, fixed = TRUE)
  formulatext <- gsub("tilde_symbol", "~", formulatext, fixed = TRUE)
  formulatext <- gsub("equal_symbol", "=", formulatext, fixed = TRUE)
  formulatext <- gsub("comma_symbol", ",", formulatext, fixed = TRUE)
  formulatext <- gsub("asterisk_symbol", "*", formulatext, fixed = TRUE)
  formulatext <- gsub("caret_symbol", "^", formulatext, fixed = TRUE)
  formula <- stats::as.formula(formulatext, env=combined_env)
  formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(formula))), env=combined_env) # here we need the formula as a 'call' object
  
  sigma.formulatext <- gsub("left_parenthesis", "(", sigma.formula, fixed = TRUE)
  sigma.formulatext <- gsub("right_parenthesis", ")", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("tilde_symbol", "~", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("equal_symbol", "=", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("comma_symbol", ",", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("asterisk_symbol", "*", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("caret_symbol", "^", sigma.formulatext, fixed = TRUE)
  sigma.formula <- stats::as.formula(sigma.formulatext, env=combined_env)
  sigma.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula))), env=combined_env) # here we need the formula as a 'call' object
  
  nu.formulatext <- gsub("left_parenthesis", "(", nu.formula, fixed = TRUE)
  nu.formulatext <- gsub("right_parenthesis", ")", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("tilde_symbol", "~", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("equal_symbol", "=", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("comma_symbol", ",", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("asterisk_symbol", "*", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("caret_symbol", "^", nu.formulatext, fixed = TRUE)
  nu.formula <- stats::as.formula(nu.formulatext, env=combined_env)
  nu.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(nu.formula))), env=combined_env) # here we need the formula as a 'call' object
  
  tau.formulatext <- gsub("left_parenthesis", "(", tau.formula, fixed = TRUE)
  tau.formulatext <- gsub("right_parenthesis", ")", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("tilde_symbol", "~", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("equal_symbol", "=", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("comma_symbol", ",", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("asterisk_symbol", "*", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("caret_symbol", "^", tau.formulatext, fixed = TRUE)
  tau.formula <- stats::as.formula(tau.formulatext, env=combined_env)
  tau.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(tau.formula))), env=combined_env) # here we need the formula as a 'call' object
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  familytext <- family
  family <- gamlss.dist::as.family(eval(parse(text=family), envir=environment()))
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  #**************************************************************************
  # II) First outer iteration of gamlss ---- 
  # To get the desired matrices and vectors (and their dimensions)
  #**************************************************************************
  
  # to increase computational speed the number of inner and backfitting iterations are set to 1
  # suppressWarnings to avoid the warning that the algorithm has not yet converged
  mod.gamlss.ds <- base::suppressWarnings(gamlss::gamlss(formula=formula2use, sigma.formula=sigma.formula2use, 
                                                         nu.formula=nu.formula2use, tau.formula=tau.formula2use,
                                                         family=family, data=dataTable,
                                                         mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix,
                                                         tau.fix=tau.fix,
                                                         control = gamlss::gamlss.control(c.crit=c1[1], n.cyc=1, 
                                                                                  mu.step=c1[3], sigma.step=c1[4], 
                                                                                  nu.step=c1[5], tau.step=c1[6],
                                                                                  gd.tol=c1[7], trace=FALSE, 
                                                                                  autostep=autostep),
                                                         i.control = gamlss::glim.control(cc=c2[1], cyc=1, 
                                                                                  bf.cyc=1, bf.tol=c2[4])))
  
  
  base::assign("temp_mod.gamlss.ds", mod.gamlss.ds, env=parent.frame())
  
  mu.x <- mod.gamlss.ds$mu.x
  sigma.x <- mod.gamlss.ds$sigma.x
  nu.x <- mod.gamlss.ds$nu.x
  tau.x <- mod.gamlss.ds$tau.x
  
  dim.mu.x <- dim(mu.x)
  dim.sigma.x <- dim(sigma.x)
  dim.nu.x <- dim(nu.x)
  dim.tau.x <- dim(tau.x)
  
  mu.coef.names <- names(mod.gamlss.ds$mu.coefficients)
  sigma.coef.names <- names(mod.gamlss.ds$sigma.coefficients)
  nu.coef.names <- names(mod.gamlss.ds$nu.coefficients)
  tau.coef.names <- names(mod.gamlss.ds$tau.coefficients)
  smoother.names <- NULL
  
  y <- as.vector(mod.gamlss.ds$y)
  
  ## Block individual level information & get variables for the smoothers
  mod.gamlss.ds$y <- "The response variable is not disclosed!"
  mod.gamlss.ds$residuals <- "The residuals of the model are not disclosed!"
  mod.gamlss.ds$control$n.cyc <- c1[2] #to get the correct maximum number of iterations
  mod.gamlss.ds$control$trace <- TRUE  #this is not relevant for ds.gamlss
  if("mu" %in% names(family$parameters)){
    mod.gamlss.ds$mu.fv <- "The fitted values of the mu model are not disclosed!"
    mod.gamlss.ds$mu.lp <- "The linear predictors of the mu model are not disclosed!"
    mod.gamlss.ds$mu.wv <- "The working variable of the mu model are not disclosed!"
    mod.gamlss.ds$mu.wt <- "The working weights of the mu model are not disclosed!"
    mod.gamlss.ds$mu.x <- "The design matrix of the mu model is not disclosed!"
    mod.gamlss.ds$mu.qr <- "The QR decomposition of the mu model is not disclosed!"
    mod.gamlss.ds$mu.formula <- ~1
    if(length(mod.gamlss.ds$mu.coefSmo)>0){
      mod.gamlss.ds$mu.s <- "The smoothing fitted values of the mu model are not disclosed!"
      mod.gamlss.ds$mu.var <- "The variances for the smoothing fitted values of the mu model are not disclosed!"
      for(i in 1:length(mod.gamlss.ds$mu.coefSmo)){
        mod.gamlss.ds$mu.coefSmo[[i]]$fv <- "The smoothing fitted values of the mu model are not disclosed!"
        mod.gamlss.ds$mu.coefSmo[[i]]$fun <- "The function for the knots of the mu model is not disclosed!"
        smoother.names <- c(smoother.names, mod.gamlss.ds$mu.coefSmo[[i]]$name)
      }
    }
  }
  if("sigma" %in% names(family$parameters)){
    mod.gamlss.ds$sigma.fv <- "The fitted values of the sigma model are not disclosed!"
    mod.gamlss.ds$sigma.lp <- "The linear predictors of the sigma model are not disclosed!"
    mod.gamlss.ds$sigma.wv <- "The working variable of the sigma model are not disclosed!"
    mod.gamlss.ds$sigma.wt <- "The working weights of the sigma model are not disclosed!"
    mod.gamlss.ds$sigma.x <- "The design matrix of the sigma model is not disclosed!"
    mod.gamlss.ds$sigma.qr <- "The QR decomposition of the sigma model is not disclosed!"
    mod.gamlss.ds$sigma.formula <- ~1
    if(length(mod.gamlss.ds$sigma.coefSmo)>0){
      mod.gamlss.ds$sigma.s <- "The smoothing fitted values of the sigma model are not disclosed!"
      mod.gamlss.ds$sigma.var <- "The variances for the smoothing fitted values of the sigma model are not disclosed!"
      for(i in 1:length(mod.gamlss.ds$sigma.coefSmo)){
        mod.gamlss.ds$sigma.coefSmo[[i]]$fv <- "The smoothing fitted values of the sigma model are not disclosed!"
        mod.gamlss.ds$sigma.coefSmo[[i]]$fun <- "The function for the knots of the sigma model is not disclosed!"
        smoother.names <- c(smoother.names, mod.gamlss.ds$sigma.coefSmo[[i]]$name)
      }
    }
  }
  if("nu" %in% names(family$parameters)){
    mod.gamlss.ds$nu.fv <- "The fitted values of the nu model are not disclosed!"
    mod.gamlss.ds$nu.lp <- "The linear predictors of the nu model are not disclosed!"
    mod.gamlss.ds$nu.wv <- "The working variable of the nu model are not disclosed!"
    mod.gamlss.ds$nu.wt <- "The working weights of the nu model are not disclosed!"
    mod.gamlss.ds$nu.x <- "The design matrix of the nu model is not disclosed!"
    mod.gamlss.ds$nu.qr <- "The QR decomposition of the nu model is not disclosed!"
    mod.gamlss.ds$nu.formula <- ~1
    if(length(mod.gamlss.ds$nu.coefSmo)>0){
      for(i in 1:length(mod.gamlss.ds$nu.coefSmo)){
        mod.gamlss.ds$nu.s <- "The smoothing fitted values of the nu model are not disclosed!"
        mod.gamlss.ds$nu.var <- "The variances for the smoothing fitted values of the nu model are not disclosed!"
        mod.gamlss.ds$nu.coefSmo[[i]]$fv <- "The smoothing fitted values of the nu model are not disclosed!"
        mod.gamlss.ds$nu.coefSmo[[i]]$fun <- "The function for the knots of the nu model is not disclosed!"
        smoother.names <- c(smoother.names, mod.gamlss.ds$nu.coefSmo[[i]]$name)
      }
    }
  }
  if("tau" %in% names(family$parameters)){
    mod.gamlss.ds$tau.fv <- "The fitted values of the tau model are not disclosed!"
    mod.gamlss.ds$tau.lp <- "The linear predictors of the tau model are not disclosed!"
    mod.gamlss.ds$tau.wv <- "The working variable of the tau model are not disclosed!"
    mod.gamlss.ds$tau.wt <- "The working weights of the tau model are not disclosed!"
    mod.gamlss.ds$tau.x <- "The design matrix of the tau model is not disclosed!"
    mod.gamlss.ds$tau.qr <- "The QR decomposition of the tau model is not disclosed!"
    mod.gamlss.ds$tau.formula <- ~1
    if(length(mod.gamlss.ds$tau.coefSmo)>0){
      mod.gamlss.ds$tau.s <- "The smoothing fitted values of the tau model are not disclosed!"
      mod.gamlss.ds$tau.var <- "The variances for the smoothing fitted values of the tau model are not disclosed!"
      for(i in 1:length(mod.gamlss.ds$tau.coefSmo)){
        mod.gamlss.ds$tau.coefSmo[[i]]$fv <- "The smoothing fitted values of the tau model are not disclosed!"
        mod.gamlss.ds$tau.coefSmo[[i]]$fun <- "The function for the knots of the tau model is not disclosed!"
        smoother.names <- c(smoother.names, mod.gamlss.ds$tau.coefSmo[[i]]$name)
      }
    }
  }
  
  # If a smoothing variable occurs several times for different distribution parameters it should only
  # occur ones in the smoother.names and minimum and maximum
  smoother.names <- unique(smoother.names)
  
  # Get the anonymized minimum and maximum for each variable in smoother.names (similar to scatterPlotDS)
  # The minimum and maximum are needed to use same knots on all servers during the fitting of the 
  # smoothing terms
  # Note that for simplicity at the moment only probabilistic anonymization is implemented
  if(length(smoother.names)>0){  # the mode includes smoothers
    smoother.xmin <- rep(NA, length(smoother.names))
    smoother.xmax <- rep(NA, length(smoother.names))
    for (i in 1:length(smoother.names)){
      x <- eval(parse(text=smoother.names[i]), envir=parent.frame())
      # the study-specific seed for random number generation
      seed <- getOption("datashield.seed")
      if (is.null(seed)){
        stop("gamlssDS1 with smoothers requires 'datashield.seed' R option to operate", call.=FALSE)
      }else{
        set.seed(seed)
        smoother.xmin[i] <- min(x) - abs(stats::rnorm(n=1, mean=0, sd=sqrt(nfilter.noise*stats::var(x))))
        smoother.xmax[i] <- max(x) + abs(stats::rnorm(n=1, mean=0, sd=sqrt(nfilter.noise*stats::var(x))))
      }
    }
  }else{  # the model does not include smoothers
    smoother.xmin <- NULL
    smoother.xmax <- NULL
  }
  
  #**************************************************************************
  # III) Initialization ----
  # Initialize the parameter vectors, smoothing fitted values & the deviance
  #**************************************************************************

  ### Initialize & save the parameter vectors on the server-side
  # since they might be disclosive they cannot be returned to the client
  mu <- NULL
  sigma <- NULL
  nu <- NULL
  tau <- NULL
  
  ## compute start values for parameters
  if("mu" %in% names(family$parameters)){
    if (familytext=="NO2()"){
      mu <- eval(family$mu.initial, envir=environment())
    } else {
      mu <- (y + global.mean)/2
    }
    base::assign("temp_mu", mu, env=parent.frame())
  }
  
  if("sigma" %in% names(family$parameters)){
    if (familytext=="NO()"){
      sigma <- rep(global.sd, length(y))
    } else if (familytext=="NO2()"){
      sigma <- rep(global.sd^2, length(y))
    } else {
      eval(family$sigma.initial, envir=environment())
    }
    base::assign("temp_sigma", sigma, env=parent.frame())
  }
  
  if("nu" %in% names(family$parameters)){
    eval(family$nu.initial, envir=environment())
    base::assign("temp_nu", nu, env=parent.frame())
  }
  
  if("tau" %in% names(family$parameters)){
    eval(family$tau.initial, envir=environment())
    base::assign("temp_tau", tau, env=parent.frame())
  }
  
  ## Initialize & save the matrix with the smoothing fitted values on the server-side
  # since they might be disclosive they cannot be returned to the client
  if (!is.null(mod.gamlss.ds$mu.coefSmo)){
    base::assign("temp_mu.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$mu.coefSmo)), env=parent.frame())
  }
  
  if (!is.null(mod.gamlss.ds$sigma.coefSmo)){
    base::assign("temp_sigma.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$sigma.coefSmo)), env=parent.frame())
  }
  
  if (!is.null(mod.gamlss.ds$nu.coefSmo)){
    base::assign("temp_nu.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$nu.coefSmo)), env=parent.frame())
  }
  
  if (!is.null(mod.gamlss.ds$tau.coefSmo)){
    base::assign("temp_tau.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$tau.coefSmo)), env=parent.frame())
  }

  ## Initialize deviance
  G.dev.incr  <- eval(body(family$G.dev.incr), envir=environment())  # deviance increment (function provided by gamlss.family object)
  G.dev <- sum(G.dev.incr)  # the weighted global deviance

  #**************************************************************************
  # IV) Disclosure risk----
  #**************************************************************************

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
    errorMessage <- "ERROR: Model has too many parameters, there is a possible risk of disclosure - please simplify model"
  }

  if(!is.null(num.mu.p)){
    if(num.mu.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage <- "ERROR: Modelfor mu has too many parameters, there is a possible risk of disclosure - please simplify model"
    }
  }

  if(!is.null(num.sigma.p)){
    if(num.sigma.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage <- "ERROR: Model for sigma has too many parameters, there is a possible risk of disclosure - please simplify model"
    }
  }

  if(!is.null(num.nu.p)){
    if(num.nu.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage <- "ERROR: Model for nu has too many parameters, there is a possible risk of disclosure - please simplify model"
    }
  }

  if(!is.null(num.tau.p)){
    if(num.tau.p > nfilter.glm*num.n){
      gamlss.saturation.invalid <- 1
      errorMessage <- "ERROR: Model for tau has too many parameters, there is a possible risk of disclosure - please simplify model"
    }
  }

  #*B) Invalid y, mu.x, sigma.x, nu.x or tau.x ----
  # If y, X or w data are invalid but user has modified client-side
  # function (ds.gamlss) to circumvent trap, model will get to this point without
  # giving a controlled shut down with a warning about invalid data.
  # So as a safety measure, we will now use the same test that is used to
  # trigger a controlled trap in the client-side function to destroy the
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
      errorMessage <- "ERROR: y vector is binary with one category less than filter threshold for table cell size"
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
          errorMessage <- "ERROR: at least one column in mu.x matrix is binary with one category less than filter threshold for table cell size"
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
          errorMessage <- "ERROR: at least one column in sigma.x matrix is binary with one category less than filter threshold for table cell size"
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
          errorMessage <- "ERROR: at least one column in nu.x matrix is binary with one category less than filter threshold for table cell size"
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
          errorMessage <- "ERROR: at least one column in tau.x matrix is binary with one category less than filter threshold for table cell size"
        }
      }
    }
  }

  #*C) Combine disclosure risks----
  # If y, mu.x, sigma.x, nu.x or tau.x are invalid, or the model is overparameterized, this will be detected by gamlssDS1
  # and passed to ds.gamlss resulting in a warning and a controlled shut down of the function.
  # But in case someone modifies the client side function to circumvent the trap, so the
  # error is only apparent once the main main iterations have started via gamlssDS2
  # the equivalent tests in gamlssDS2 will destroy the info.matrix and score.vector in the affected study so
  # the model fitting will simply terminate.
  if(!(y.invalid>0 || sum(mu.par.invalid)>0|| sum(sigma.par.invalid)>0 || sum(nu.par.invalid)>0 ||
         sum(tau.par.invalid)>0 || gamlss.saturation.invalid>0)){
    errorMessage <- "No errors"
  }else{
    errorMessage <- "Study data or applied model invalid for this source"
    mod.gamlss.ds <- NA
    G.dev <- NA
    dim.mu.x <- NA
    dim.sigma.x <- NA
    dim.nu.x <- NA
    dim.tau.x <- NA
    smoother.names <- NA 
    smoother.xmin  <- NA
    smoother.xmax <- NA
    rm(list=ls(pattern="^temp_", envir=parent.frame()), envir=parent.frame())
  }

  #**************************************************************************
  # V) Output ----
  #**************************************************************************
  return(list(mod.gamlss.ds=mod.gamlss.ds, G.dev=G.dev,
              dim.mu.x=dim.mu.x, dim.sigma.x=dim.sigma.x, dim.nu.x=dim.nu.x, dim.tau.x=dim.tau.x,
              smoother.names=smoother.names, smoother.xmin=smoother.xmin, smoother.xmax=smoother.xmax,
              y.invalid=y.invalid, mu.par.invalid=mu.par.invalid, sigma.par.invalid=sigma.par.invalid,
              nu.par.invalid=nu.par.invalid, tau.par.invalid=tau.par.invalid,
              gamlss.saturation.invalid=gamlss.saturation.invalid, errorMessage=errorMessage))

}
# AGGREGATE FUNCTION
# gamlssDS1
