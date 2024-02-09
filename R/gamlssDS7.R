#'
#' @title gamlssDS7 called by ds.gamlss
#' @description This is the seventh serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that calculates the normalized quantile
#' residuals. For more details please see the extensive header of ds.gamlss and also the 
#' gamlss function in native R gamlss package.
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
#' @return anonymized normalized quantile residuals for the gamlss model fitted with 
#' ds.gamlss. 
#' @author Annika Swenne
#' @import gamlss
#' @import gamlss.dist
#' @export

gamlssDS7 <- function(formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula,
                      tau.formula = tau.formula, family = family, data=data, mu.fix=mu.fix,
                      sigma.fix = sigma.fix, nu.fix = nu.fix, tau.fix = tau.fix,
                      control = control, i.control = i.control){
  
  #**************************************************************************
  # I) Preparation ---- 
  #**************************************************************************
  
  ## Define internal function as in gamlss that is used to estimate normalized quantile residuals
  rqres <- function (pfun = "pNO", type = c("Continuous", "Discrete", "Mixed"), 
                     censored = NULL, ymin = NULL, mass.p = NULL, prob.mp = NULL, 
                     y = y, ...) 
  {
    type <- match.arg(type)
    cdf <- eval(parse(text = pfun))
    switch(type, Continuous = {
      rqres <- qnorm(cdf(q = y, ...))
    }, Discrete = {
      if (is.null(censored)) {
        aval <- cdf(ifelse(y == ymin, y, y - 1), ...)
        aval <- ifelse(y == ymin, 0, aval)
        bval <- cdf(q = y, ...)
        uval <- runif(length(y), aval, bval)
        uval <- ifelse(uval > 0.999999, uval - 1e-16, uval)
        uval <- ifelse(uval < 1e-06, uval + 1e-16, uval)
        rqres <- qnorm(uval)
      } else {
        qq <- ifelse(y[, 1] == ymin, y[, 1], y[, 1] - 1)
        aval <- cdf(Surv(qq, y[, 2]), ...)
        aval <- ifelse(y[, 1] == ymin, 0, aval)
        bval <- cdf(q = y, ...)
        uval <- runif(length(y[, 1]), min = aval, max = bval)
        uval <- ifelse(uval > 0.999999, uval - 1e-16, uval)
        uval <- ifelse(uval < 1e-06, uval + 1e-16, uval)
        rqres <- qnorm(ifelse(y[, "status"] == 1, uval, bval))
      }
    }, Mixed = {
      if (is.null(mass.p) && is.null(prob.mp)) stop("For mixed distributions mass.p and prob.mp arguments have to be specified")
      length.mass.p <- length(mass.p)
      switch(length.mass.p, {
        if (mass.p == 0) {
          uval <- ifelse(y == mass.p, runif(length(y), 
                                            0, prob.mp), cdf(q = y, ...))
        } else if (mass.p == 1) {
          uval <- ifelse(y == mass.p, runif(length(y), 
                                            1 - prob.mp, 1), cdf(q = y, ...))
        } else {
          stop("mass point is not at zero or one")
        }
      }, {
        uval <- ifelse(y == mass.p[1], runif(length(y), 0, 
                                             prob.mp[, 1]), 0)
        uval <- ifelse(y > mass.p[1] & y < mass.p[2], cdf(q = y, 
                                                          ...), uval)
        uval <- ifelse(y == 1, runif(length(y), 1 - prob.mp[, 
                                                            2], 1), uval)
      })
      rqres <- qnorm(uval)
    })
    rqres
  }
  
  body(rqres) <-  eval(quote(body(rqres)), envir = parent.frame())
  
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
  mu.step <- c1[3]
  sigma.step <- c1[4]
  nu.step <- c1[5]
  tau.step <- c1[6]
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  #**************************************************************************
  # II) Get the outcome ----  
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
  
  ## get outcome
  y <- as.vector(mod.gamlss.ds$y)
  
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
  
  #**************************************************************************
  # III) Calculate residuals ----  
  #**************************************************************************
  
  residuals <- eval(family$rqres, env=parent.frame())
  return(residuals)
}