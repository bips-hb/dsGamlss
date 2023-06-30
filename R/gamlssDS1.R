#'
#' @title gamlssDS1 called by ds.gamlss
#' @description This is the first serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that sets up the model structure and creates
#' the starting beta vectors that feeds, via ds.gamlss into gamlssDS2 to enable iterative
#' fitting of the gamlss model that has been specified. For more details please see the 
#' extensive header of ds.gamlss and also the gamlss function in native R gamlss package.
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

gamlssDS1 <- function(formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula,
                     tau.formula = tau.formula, family = family, data=data,
                     control = control, i.control = i.control){
  
  
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
    dataTable <- eval(parse(text=data), env=parent.frame())
  }
  
  ## Reconvert the special symbols to create the appropriate formula & gamlss.family objects
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
  
  #**************************************************************************
  # II) First outer iteration of gamlss ---- 
  # To get the desired matrices and vectors (and their dimensions)
  #**************************************************************************

  # to increase computational speed the number of inner and backfitting iterations are set to 1
  # suppressWarnings to avoid the warning that the algorithm has not yet converged
  mod.gamlss.ds <- base::suppressWarnings(gamlss::gamlss(formula=formula2use, sigma.formula=sigma.formula2use, 
                                                         nu.formula=nu.formula2use, tau.formula=tau.formula2use,
                                                         family=family, data=dataTable, method=RS(),
                                                         mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix,
                                                         tau.fix=tau.fix,
                                                         control = gamlss.control(c.crit=c1[1], n.cyc=1, 
                                                                                  mu.step=c1[3], sigma.step=c1[4], 
                                                                                  nu.step=c1[5], tau.step=c1[6],
                                                                                  gd.tol=c1[7], trace=FALSE),
                                                         i.control = glim.control(cc=c2[1], cyc=1, 
                                                                                  bf.cyc=1, bf.tol=c2[4])))
  
  parameters <- mod.gamlss.ds$parameters
  
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
  
  y <- as.vector(mod.gamlss.ds$y)

  ## Smoothers
  smoother.names <- c(mu.coef.names, sigma.coef.names, nu.coef.names, tau.coef.names)
  pb.names <- smoother.names[grep(pattern="pb(", x=smoother.names, fixed=TRUE)]
  pb.names <- unique(pb.names)
  pb.names <- gsub(pattern="pb(", replacement="", pb.names, fixed=TRUE)
  pb.names <- gsub(pattern=")", replacement="", pb.names, fixed=TRUE)
  
  # Get the anonymized minimum and maximum for each variable in pb.names (similar to scatterPlotDS)
  # The minimum and maximum are needed to use same knots on all servers during the fitting of the 
  # smoothing terms
  # Note that for simplicity at the moment only probabilistic anonymization is implemented
  if(length(pb.names)>0){  # the mode includes pb-smoothers
    pb.xmin <- rep(NA, length(pb.names))
    pb.xmax <- rep(NA, length(pb.names))
    for (i in 1:length(pb.names)){
      x <- eval(parse(text=pb.names[i]), env=parent.frame())
      # the study-specific seed for random number generation
      seed <- getOption("datashield.seed")
      if (is.null(seed)){
        stop("gamlssDS1 with pb-smoothers requires 'datashield.seed' R option to operate", call.=FALSE)
      }else{
        set.seed(seed)
        x.new <- x + stats::rnorm(n=length(x), mean=0, sd=sqrt(nfilter.noise*stats::var(x)))
        pb.xmin[i] <- min(x.new)
        pb.xmax[i] <- max(x.new)
      }
    }
  }else{  # the model does not include pb-smoothers
    pb.xmin <- NULL
    pb.xmax <- NULL
  }

  # get the length of the gamma vectors
  # mu
  mu.gamma.length <- NULL
  if (!is.null(mod.gamlss.ds$mu.coefSmo)){
    for (s in 1:length(mod.gamlss.ds$mu.coefSmo)){
      mu.gamma.length <- c(mu.gamma.length, dim(mod.gamlss.ds$mu.coefSmo[[s]]$coef)[1])
    }
  }
  
  # sigma
  sigma.gamma.length <- NULL
  if (!is.null(mod.gamlss.ds$sigma.coefSmo)){
    for (s in 1:length(mod.gamlss.ds$sigma.coefSmo)){
      sigma.gamma.length <- c(sigma.gamma.length, dim(mod.gamlss.ds$sigma.coefSmo[[s]]$coef)[1])
    }
  }
  
  # nu
  nu.gamma.length <- NULL
  if (!is.null(mod.gamlss.ds$nu.coefSmo)){
    for (s in 1:length(mod.gamlss.ds$nu.coefSmo)){
      nu.gamma.length <- c(nu.gamma.length, dim(mod.gamlss.ds$nu.coefSmo[[s]]$coef)[1])
    }
  }
  
  # tau
  tau.gamma.length <- NULL
  if (!is.null(mod.gamlss.ds$tau.coefSmo)){
    for (s in 1:length(mod.gamlss.ds$tau.coefSmo)){
      tau.gamma.length <- c(tau.gamma.length, dim(mod.gamlss.ds$tau.coefSmo[[s]]$coef)[1])
    }
  }
  
  #**************************************************************************
  # III) Initialization ----
  # Initialize the parameter vectors, smoothing fitted values & the deviance
  #**************************************************************************

  ## Initialize & save the parameter vectors on the server-side
  # since they might be disclosive they cannot be returned to the client
  if("mu" %in% names(family$parameters)){
    eval(family$mu.initial, env=environment())
    base::assign("mu", mu, env=parent.frame())
  }
  if("sigma" %in% names(family$parameters)){
    eval(family$sigma.initial, env=environment())
    base::assign("sigma", sigma, env=parent.frame())
  }
  if("nu" %in% names(family$parameters)){
    eval(family$nu.initial, env=environment())
    base::assign("nu", nu, env=parent.frame())
  }
  if("tau" %in% names(family$parameters)){
    eval(family$tau.initial, env=environment())
    base::assign("tau", tau, env=parent.frame())
  }
  
  ## Initialize & save the matrix with the smoothing fitted values on the server-side
  # since they might be disclosive they cannot be returned to the client
  if (!is.null(mod.gamlss.ds$mu.coefSmo)){
    base::assign("mu.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$mu.coefSmo)), env=parent.frame())
  }
  
  if (!is.null(mod.gamlss.ds$sigma.coefSmo)){
    base::assign("sigma.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$sigma.coefSmo)), env=parent.frame())
  }
  
  if (!is.null(mod.gamlss.ds$nu.coefSmo)){
    base::assign("nu.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$nu.coefSmo)), env=parent.frame())
  }
  
  if (!is.null(mod.gamlss.ds$tau.coefSmo)){
    base::assign("tau.s", matrix(0, mod.gamlss.ds$N, length(mod.gamlss.ds$tau.coefSmo)), env=parent.frame())
  }

  ## Initialize deviance
  G.dev.incr  <- eval(body(family$G.dev.incr), env=environment())  # deviance increment (function provided by gamlss.family object)
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
  }

  #**************************************************************************
  # V) Output ----
  #**************************************************************************
  return(list(G.dev=G.dev, parameters=parameters,
              dim.mu.x=dim.mu.x, dim.sigma.x=dim.sigma.x, dim.nu.x=dim.nu.x, dim.tau.x=dim.tau.x,
              mu.coef.names=mu.coef.names, sigma.coef.names=sigma.coef.names, nu.coef.names=nu.coef.names, tau.coef.names=tau.coef.names,
              pb.names=pb.names, pb.xmin=pb.xmin, pb.xmax=pb.xmax,
              mu.gamma.length=mu.gamma.length, sigma.gamma.length=sigma.gamma.length,
              nu.gamma.length=nu.gamma.length, tau.gamma.length=tau.gamma.length,
              y.invalid=y.invalid, mu.par.invalid=mu.par.invalid, sigma.par.invalid=sigma.par.invalid,
              nu.par.invalid=nu.par.invalid, tau.par.invalid=tau.par.invalid,
              gamlss.saturation.invalid=gamlss.saturation.invalid, errorMessage=errorMessage))

}
# AGGREGATE FUNCTION
# gamlssDS1
