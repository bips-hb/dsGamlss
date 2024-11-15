#'
#' @title gamlssDS1 called by ds.gamlss
#' @description This is the first serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that sets up the model structure and creates
#' the starting beta and gamma vectors that feeds, via ds.gamlss into gamlssDS2 to enable iterative
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
#' @param mu.coef.start optional vector of regression coefficients to compute improved 
#' start values for mu. Default NULL.
#' @param sigma.coef.start optional vector of regression coefficients to compute improved
#' start values for sigma. Default NULL.
#' @param nu.coef.start optional vector of regression coefficients to compute improved
#' start values for nu. Default NULL.
#' @param tau.coef.start optional vector of regression coefficients to compute improved
#' start values for tau. Default NULL.
#' @param mu.coef.start.names vector with names for the regression coefficients in 
#' \code{mu.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for mu. If values are given in \code{mu.coef.start} but 
#' \code{mu.coef.start.names} is NULL then the same formula as in \code{formula} is used to
#' obtain the design matrix. Default NULL.
#' @param sigma.coef.start.names vector with names for the regression coefficients in 
#' \code{sigma.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for sigma. If values are given in \code{sigma.coef.start} but 
#' \code{sigma.coef.start.names} is NULL then the same formula as in \code{sigma.formula} 
#' is used to obtain the design matrix. Default NULL.
#' @param nu.coef.start.names vector with names for the regression coefficients in 
#' \code{nu.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for nu. If values are given in \code{nu.coef.start} but 
#' \code{nu.coef.start.names} is NULL then the same formula as in \code{nu.formula} 
#' is used to obtain the design matrix. Default NULL.
#' @param tau.coef.start.names vector with names for the regression coefficients in 
#' \code{tau.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for tau. If values are given in \code{tau.coef.start} but 
#' \code{tau.coef.start.names} is NULL then the same formula as in \code{tau.formula} 
#' is used to obtain the design matrix. Default NULL.
#' @param mu.fix logical, indicating whether the mu parameter should be kept fixed
#' in the fitting processes.
#' @param sigma.fix logical, indicating whether the sigma parameter should be kept
#' fixed in the fitting processes.
#' @param nu.fix logical, indicating whether the nu parameter should be kept fixed 
#' in the fitting processes.
#' @param tau.fix logical, indicating whether the tau parameter should be kept fixed
#' in the fitting processes.
#' @param global.mean numeric value, giving the global mean of the outcome variable
#' (necessary to initialize the distribution parameter for some families, otherwise NULL)
#' @param global.sd numeric value, giving the global sd of the outcome variable
#' (necessary to initialize the distribution parameter for some families, otherwise NULL)
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
#' @export
#'

gamlssDS1 <- function(formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula,
                     tau.formula = tau.formula, family = family, data = data, mu.coef.start = mu.coef.start, 
                     sigma.coef.start = sigma.coef.start, nu.coef.start = nu.coef.start,
                     tau.coef.start = tau.coef.start, mu.coef.start.names = mu.coef.start.names, 
                     sigma.coef.start.names = sigma.coef.start.names, nu.coef.start.names = nu.coef.start.names,
                     tau.coef.start.names = tau.coef.start.names, mu.fix=mu.fix, sigma.fix = sigma.fix, 
                     nu.fix = nu.fix, tau.fix = tau.fix, global.mean = global.mean, 
                     global.sd = global.sd, control = control, i.control = i.control){
  
  
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
  
  ## Reconvert the special symbols to create the appropriate formula & gamlss.family objects
  formulatext <- gsub("left_parenthesis", "(", formula, fixed = TRUE)
  formulatext <- gsub("right_parenthesis", ")", formulatext, fixed = TRUE)
  formulatext <- gsub("tilde_symbol", "~", formulatext, fixed = TRUE)
  formulatext <- gsub("equal_symbol", "=", formulatext, fixed = TRUE)
  formulatext <- gsub("comma_symbol", ",", formulatext, fixed = TRUE)
  formulatext <- gsub("asterisk_symbol", "*", formulatext, fixed = TRUE)
  formulatext <- gsub("caret_symbol", "^", formulatext, fixed = TRUE)
  formula <- stats::as.formula(formulatext)
  formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  sigma.formulatext <- gsub("left_parenthesis", "(", sigma.formula, fixed = TRUE)
  sigma.formulatext <- gsub("right_parenthesis", ")", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("tilde_symbol", "~", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("equal_symbol", "=", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("comma_symbol", ",", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("asterisk_symbol", "*", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("caret_symbol", "^", sigma.formulatext, fixed = TRUE)
  sigma.formula <- stats::as.formula(sigma.formulatext)
  sigma.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  nu.formulatext <- gsub("left_parenthesis", "(", nu.formula, fixed = TRUE)
  nu.formulatext <- gsub("right_parenthesis", ")", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("tilde_symbol", "~", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("equal_symbol", "=", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("comma_symbol", ",", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("asterisk_symbol", "*", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("caret_symbol", "^", nu.formulatext, fixed = TRUE)
  nu.formula <- stats::as.formula(nu.formulatext)
  nu.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(nu.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  tau.formulatext <- gsub("left_parenthesis", "(", tau.formula, fixed = TRUE)
  tau.formulatext <- gsub("right_parenthesis", ")", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("tilde_symbol", "~", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("equal_symbol", "=", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("comma_symbol", ",", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("asterisk_symbol", "*", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("caret_symbol", "^", tau.formulatext, fixed = TRUE)
  tau.formula <- stats::as.formula(tau.formulatext)
  tau.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(tau.formula))), env=parent.frame()) # here we need the formula as a 'call' object
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  familytext <- family
  family <- gamlss.dist::as.family(eval(parse(text=family), envir=environment()))
  
  if(!is.null(mu.coef.start.names)){
    mu.coef.start.names <- gsub("left_parenthesis", "(", mu.coef.start.names, fixed = TRUE)
    mu.coef.start.names <- gsub("right_parenthesis", ")", mu.coef.start.names, fixed = TRUE)
    mu.coef.start.names <- gsub("asterisk_symbol", "*", mu.coef.start.names, fixed = TRUE)
    mu.coef.start.names <- gsub("caret_symbol", "^", mu.coef.start.names, fixed = TRUE)
    mu.coef.start.names <- unlist(strsplit(mu.coef.start.names, split=","))
  }
  
  if(!is.null(sigma.coef.start.names)){
    sigma.coef.start.names <- gsub("left_parenthesis", "(", sigma.coef.start.names, fixed = TRUE)
    sigma.coef.start.names <- gsub("right_parenthesis", ")", sigma.coef.start.names, fixed = TRUE)
    sigma.coef.start.names <- gsub("asterisk_symbol", "*", sigma.coef.start.names, fixed = TRUE)
    sigma.coef.start.names <- gsub("caret_symbol", "^", sigma.coef.start.names, fixed = TRUE)
    sigma.coef.start.names <- unlist(strsplit(sigma.coef.start.names, split=","))
  }
  
  if(!is.null(nu.coef.start.names)){
    nu.coef.start.names <- gsub("left_parenthesis", "(", nu.coef.start.names, fixed = TRUE)
    nu.coef.start.names <- gsub("right_parenthesis", ")", nu.coef.start.names, fixed = TRUE)
    nu.coef.start.names <- gsub("asterisk_symbol", "*", nu.coef.start.names, fixed = TRUE)
    nu.coef.start.names <- gsub("caret_symbol", "^", nu.coef.start.names, fixed = TRUE)
    nu.coef.start.names <- unlist(strsplit(nu.coef.start.names, split=","))
  }
  
  if(!is.null(tau.coef.start.names)){
    tau.coef.start.names <- gsub("left_parenthesis", "(", tau.coef.start.names, fixed = TRUE)
    tau.coef.start.names <- gsub("right_parenthesis", ")", tau.coef.start.names, fixed = TRUE)
    tau.coef.start.names <- gsub("asterisk_symbol", "*", tau.coef.start.names, fixed = TRUE)
    tau.coef.start.names <- gsub("caret_symbol", "^", tau.coef.start.names, fixed = TRUE)
    tau.coef.start.names <- unlist(strsplit(tau.coef.start.names, split=","))
  }
  
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
                                                                                  gd.tol=c1[7], trace=FALSE),
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
  if("mu" %in% names(family$parameters)){
    mod.gamlss.ds$mu.fv <- "The fitted values of the mu model are not disclosed!"
    mod.gamlss.ds$mu.lp <- "The linear predictors of the mu model are not disclosed!"
    mod.gamlss.ds$mu.wv <- "The working variable of the mu model are not disclosed!"
    mod.gamlss.ds$mu.wt <- "The working weights of the mu model are not disclosed!"
    mod.gamlss.ds$mu.terms <- NULL
    mod.gamlss.ds$mu.x <- "The design matrix of the mu model is not disclosed!"
    mod.gamlss.ds$mu.qr <- "The QR decomposition of the mu model is not disclosed!"
    mod.gamlss.ds$mu.formula <- ~1
    mod.gamlss.ds$mu.s <- "The smoothing fitted values of the mu model are not disclosed!"
    mod.gamlss.ds$mu.var <- "The variances for the smoothing fitted values of the mu model are not disclosed!"
    if(length(mod.gamlss.ds$mu.coefSmo)>0){
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
    mod.gamlss.ds$sigma.terms <- NULL
    mod.gamlss.ds$sigma.x <- "The design matrix of the sigma model is not disclosed!"
    mod.gamlss.ds$sigma.qr <- "The QR decomposition of the sigma model is not disclosed!"
    mod.gamlss.ds$sigma.formula <- ~1
    mod.gamlss.ds$sigma.s <- "The smoothing fitted values of the sigma model are not disclosed!"
    mod.gamlss.ds$sigma.var <- "The variances for the smoothing fitted values of the sigma model are not disclosed!"
    if(length(mod.gamlss.ds$sigma.coefSmo)>0){
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
    mod.gamlss.ds$nu.terms <- NULL
    mod.gamlss.ds$nu.x <- "The design matrix of the nu model is not disclosed!"
    mod.gamlss.ds$nu.qr <- "The QR decomposition of the nu model is not disclosed!"
    mod.gamlss.ds$nu.formula <- ~1
    mod.gamlss.ds$nu.s <- "The smoothing fitted values of the nu model are not disclosed!"
    mod.gamlss.ds$nu.var <- "The variances for the smoothing fitted values of the nu model are not disclosed!"
    if(length(mod.gamlss.ds$nu.coefSmo)>0){
      for(i in 1:length(mod.gamlss.ds$nu.coefSmo)){
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
    mod.gamlss.ds$tau.terms <- NULL
    mod.gamlss.ds$tau.x <- "The design matrix of the tau model is not disclosed!"
    mod.gamlss.ds$tau.qr <- "The QR decomposition of the tau model is not disclosed!"
    mod.gamlss.ds$tau.formula <- ~1
    mod.gamlss.ds$tau.s <- "The smoothing fitted values of the tau model are not disclosed!"
    mod.gamlss.ds$tau.var <- "The variances for the smoothing fitted values of the tau model are not disclosed!"
    if(length(mod.gamlss.ds$tau.coefSmo)>0){
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
  
  ## Create the design matrix to compute the start values
  if("mu" %in% names(family$parameters)){
    if(!is.null(mu.coef.start.names)){
      # Check for presence of any form of intercept in mu.coef.start.names
      mu.include.intercept <- any(grepl("\\bintercept\\b", mu.coef.start.names, ignore.case = TRUE))
      # Remove any form of "intercept" from the terms if it exists
      mu.terms <- mu.coef.start.names[!grepl("\\bintercept\\b", mu.coef.start.names, ignore.case = TRUE)]
      # Create the formula string based on the presence of intercept
      if (mu.include.intercept==TRUE) {
        if (length(mu.terms) > 0){
          mu.formula.start <- stats::as.formula(paste("~", paste(mu.terms, collapse = " + ")))  # includes intercept by default
        } else {
          # only intercept
          mu.formula.start <- stats::as.formula("~ 1")
        }
      } else {
        mu.formula.start <- stats::as.formula(paste("~ 0 +", paste(mu.terms, collapse = " + ")))  # no intercept
      }
      # Here we need the formula as a 'call' object
      mu.formula.start <- stats::as.formula(paste0(Reduce(paste, deparse(mu.formula.start))), env=parent.frame())
      # Generate the design matrix
      mu.x.start <- stats::model.matrix(mu.formula.start, data = dataTable)
      mu.x.start.source <- "mu.coef.start.names"
      
    } else {
      mu.x.start <- mu.x
      mu.x.start.source <- "the formula for mu"
    }
  }
  
  if("sigma" %in% names(family$parameters)){
    if(!is.null(sigma.coef.start.names)){
      # Check for presence of any form of intercept in sigma.coef.start.names
      sigma.include.intercept <- any(grepl("\\bintercept\\b", sigma.coef.start.names, ignore.case = TRUE))
      # Remove any form of "intercept" from the terms if it exists
      sigma.terms <- sigma.coef.start.names[!grepl("\\bintercept\\b", sigma.coef.start.names, ignore.case = TRUE)]
      # Create the formula string based on the presence of intercept
      if (sigma.include.intercept==TRUE) {
        if (length(sigma.terms) > 0){
          sigma.formula.start <- stats::as.formula(paste("~", paste(sigma.terms, collapse = " + ")))  # includes intercept by default
        } else {
          # only intercept
          sigma.formula.start <- stats::as.formula("~ 1")
        }
      } else {
        sigma.formula.start <- stats::as.formula(paste("~ 0 +", paste(sigma.terms, collapse = " + ")))  # no intercept
      }
      # Here we need the formula as a 'call' object
      sigma.formula.start <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula.start))), env=parent.frame())
      # Generate the design matrix
      sigma.x.start <- stats::model.matrix(sigma.formula.start, data = dataTable)
      sigma.x.start.source <- "sigma.coef.start.names"
    } else {
      sigma.x.start <- sigma.x
      sigma.x.start.source <- "the sigma.formula"
    }
  }
  
  if("nu" %in% names(family$parameters)){
    if(!is.null(nu.coef.start.names)){
      # Check for presence of any form of intercept in nu.coef.start.names
      nu.include.intercept <- any(grepl("\\bintercept\\b", nu.coef.start.names, ignore.case = TRUE))
      # Remove any form of "intercept" from the terms if it exists
      nu.terms <- nu.coef.start.names[!grepl("\\bintercept\\b", nu.coef.start.names, ignore.case = TRUE)]
      # Create the formula string based on the presence of intercept
      if (nu.include.intercept==TRUE) {
        if (length(nu.terms) > 0){
          nu.formula.start <- stats::as.formula(paste("~", paste(nu.terms, collapse = " + ")))  # includes intercept by default
        } else {
          # only intercept
          nu.formula.start <- stats::as.formula("~ 1")
        }
      } else {
        nu.formula.start <- stats::as.formula(paste("~ 0 +", paste(nu.terms, collapse = " + ")))  # no intercept
      }
      # Here we need the formula as a 'call' object
      nu.formula.start <- stats::as.formula(paste0(Reduce(paste, deparse(nu.formula.start))), env=parent.frame())
      # Generate the design matrix
      nu.x.start <- stats::model.matrix(nu.formula.start, data = dataTable)
      nu.x.start.source <- "nu.coef.start.names"
    } else {
      nu.x.start <- nu.x
      nu.x.start.source <- "the nu.formula"
    }
  }
  
  if("tau" %in% names(family$parameters)){
    if(!is.null(tau.coef.start.names)){
      # Check for presence of any form of intercept in tau.coef.start.names
      tau.include.intercept <- any(grepl("\\bintercept\\b", tau.coef.start.names, ignore.case = TRUE))
      # Remove any form of "intercept" from the terms if it exists
      tau.terms <- tau.coef.start.names[!grepl("\\bintercept\\b", tau.coef.start.names, ignore.case = TRUE)]
      # Create the formula string based on the presence of intercept
      if (tau.include.intercept==TRUE) {
        if (length(tau.terms) > 0){
          tau.formula.start <- stats::as.formula(paste("~", paste(tau.terms, collapse = " + ")))  # includes intercept by default
        } else {
          # only intercept
          tau.formula.start <- stats::as.formula("~ 1")
        }
      } else {
        tau.formula.start <- stats::as.formula(paste("~ 0 +", paste(tau.terms, collapse = " + ")))  # no intercept
      }
      # Here we need the formula as a 'call' object
      sigma.formula.start <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula.start))), env=parent.frame())
      # Generate the design matrix
      tau.x.start <- stats::model.matrix(tau.formula.start, data = dataTable)
      tau.x.start.source <- "tau.coef.start.names"
    } else {
      tau.x.start <- tau.x
      tau.x.start.source <- "the tau.formula"
    }
  }
  
  ## compute start values for parameters
  if("mu" %in% names(family$parameters)){
    if(!is.null(mu.coef.start)){ 
      if(length(mu.coef.start)==dim(mu.x.start)[2]){
        eta <- mu.x.start %*% mu.coef.start
        mu <- as.vector(eval(parse(text="family$mu.linkinv(eta)"), envir=environment()))
      } else{
        warning(paste("The length of mu.coef.start does not match with the length of", dim(mu.x.start)[2], "that is implied by ", mu.x.start.source, ". 
                      Therefore, mu.coef.start is ignored and instead the default values are used for the initialization of mu.", sep=" "))
        if (familytext=="NO2()"){
          mu <- eval(family$mu.initial, envir=environment())
        } else {
          mu <- (y + global.mean)/2
        }
      }
    } else{
      if (familytext=="NO2()"){
        mu <- eval(family$mu.initial, envir=environment())
      } else {
        mu <- (y + global.mean)/2
      }
    }
    base::assign("temp_mu", mu, env=parent.frame())
  }
  
  if("sigma" %in% names(family$parameters)){
    if(!is.null(sigma.coef.start)){
      if(length(sigma.coef.start)==dim(sigma.x.start)[2]){
        eta <- sigma.x.start %*% sigma.coef.start
        sigma <- as.vector(eval(parse(text="family$sigma.linkinv(eta)"), envir=environment()))
      } else{
        warning(paste("The length of sigma.coef.start does not match with the length of", dim(sigma.x.start)[2], "that is implied by ", sigma.x.start.source, ". 
                      Therefore, sigma.coef.start is ignored and instead the default values are used for the initialization of sigma.", sep=" "))
        if (familytext=="NO()"){
          sigma <- rep(global.sd, length(y))
        } else if (familytext=="NO2()"){
          sigma <- rep(global.sd^2, length(y))
        } else {
          eval(family$sigma.initial, envir=environment())
        }
      }
    } else {
      if (familytext=="NO()"){
        sigma <- rep(global.sd, length(y))
      } else if (familytext=="NO2()"){
        sigma <- rep(global.sd^2, length(y))
      } else {
        eval(family$sigma.initial, envir=environment())
      }
    }
    base::assign("temp_sigma", sigma, env=parent.frame())
  }
  
  if("nu" %in% names(family$parameters)){
    if(!is.null(nu.coef.start)){
      if(length(nu.coef.start)==dim(nu.x.start)[2]){
        eta <- nu.x.start %*% nu.coef.start
        nu <- as.vector(eval(parse(text="family$nu.linkinv(eta)"), envir=environment()))
      } else{
        warning(paste("The length of nu.coef.start does not match with the length of", dim(nu.x.start)[2], "that is implied by ", nu.x.start.source, ". 
                      Therefore, nu.coef.start is ignored and instead the default values are used for the initialization of nu.", sep=" "))
        eval(family$nu.initial, envir=environment())
      }
    } else{
      eval(family$nu.initial, envir=environment())
    }
    base::assign("temp_nu", nu, env=parent.frame())
  }
  
  if("tau" %in% names(family$parameters)){
    if(!is.null(tau.coef.start)){
      if(length(tau.coef.start)==dim(tau.x.start)[2]){
        eta <- tau.x.start %*% tau.coef.start
        tau <- as.vector(eval(parse(text="family$tau.linkinv(eta)"), envir=environment()))
      } else{
        warning(paste("The length of tau.coef.start does not match with the length of", dim(tau.x.start)[2], "that is implied by ", tau.x.start.source, ".
                      Therefore, tau.coef.start is ignored and instead the default values are used for the initialization of tau.", sep=" "))
        eval(family$tau.initial, envir=environment())
      }
    } else{
      eval(family$tau.initial, envir=environment())
    }
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
  return(list(mod.gamlss.ds=mod.gamlss.ds, G.dev=G.dev,
              dim.mu.x=dim.mu.x, dim.sigma.x=dim.sigma.x, dim.nu.x=dim.nu.x, dim.tau.x=dim.tau.x,
              smoother.names=smoother.names, smoother.xmin=smoother.xmin, smoother.xmax=smoother.xmax,
              y.invalid=y.invalid, mu.par.invalid=mu.par.invalid, sigma.par.invalid=sigma.par.invalid,
              nu.par.invalid=nu.par.invalid, tau.par.invalid=tau.par.invalid,
              gamlss.saturation.invalid=gamlss.saturation.invalid, errorMessage=errorMessage))

}
# AGGREGATE FUNCTION
# gamlssDS1
