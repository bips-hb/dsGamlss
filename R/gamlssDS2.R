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
#' @return a gamlss object with all components as in the native R gamlss function. 
#' Individual-level information like the components y (the response response) and 
#' residuals (the normalised quantile residuals of the model) are not disclosed to 
#' the client-side.
#' @author Annika Swenne
#' @import gamlss
#' @import gamlss.dist
#' @export
#'

gamlssDS2 <- function(formula = formula, sigma.formula = sigma.formula, nu.formula = nu.formula,
                      tau.formula = tau.formula, family = family, data=data, mu.fix = mu.fix,  
                      sigma.fix = sigma.fix, nu.fix = nu.fix, tau.fix = tau.fix,
                      mu.beta.vect = mu.beta.vect, sigma.beta.vect = sigma.beta.vect,
                      nu.beta.vect = nu.beta.vect, tau.beta.vect = tau.beta.vect,
                      control = control, i.control = i.control){
  #**************************************************************************
  # I) Preparation ---- 
  # Reconvert the transfer strings into required variable types
  #**************************************************************************
  
  ## Capture the nfilter settings
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  
  ## Get the value of the 'data' parameter provided as character on the client side
  if(is.null(data)){
    dataTable <- NULL 
  }else{
    dataTable <- eval(parse(text=data), envir = parent.frame())
  }
  
  ## Reconvert the special symbols to create the appropriate formula & gamlss.family objects
  formulatext <- gsub("left_parenthesis", "(", formula, fixed = TRUE)
  formulatext <- gsub("right_parenthesis", ")", formulatext, fixed = TRUE)
  formulatext <- gsub("tilde_symbol", "~", formulatext, fixed = TRUE)
  formulatext <- gsub("equal_symbol", "=", formulatext, fixed = TRUE)
  formulatext <- gsub("comma_symbol", ",", formulatext, fixed = TRUE)
  formula <- stats::as.formula(formulatext)
  formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(formula))), env = parent.frame()) # here we need the formula as a 'call' object
  original.formulatext <- formulatext
  
  sigma.formulatext <- gsub("left_parenthesis", "(", sigma.formula, fixed = TRUE)
  sigma.formulatext <- gsub("right_parenthesis", ")", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("tilde_symbol", "~", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("equal_symbol", "=", sigma.formulatext, fixed = TRUE)
  sigma.formulatext <- gsub("comma_symbol", ",", sigma.formulatext, fixed = TRUE)
  sigma.formula <- stats::as.formula(sigma.formulatext)
  sigma.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(sigma.formula))), env = parent.frame()) # here we need the formula as a 'call' object
  sigma.original.formulatext <- sigma.formulatext
  
  nu.formulatext <- gsub("left_parenthesis", "(", nu.formula, fixed = TRUE)
  nu.formulatext <- gsub("right_parenthesis", ")", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("tilde_symbol", "~", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("equal_symbol", "=", nu.formulatext, fixed = TRUE)
  nu.formulatext <- gsub("comma_symbol", ",", nu.formulatext, fixed = TRUE)
  nu.formula <- stats::as.formula(nu.formulatext)
  nu.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(nu.formula))), env = parent.frame()) # here we need the formula as a 'call' object
  nu.original.formulatext <- nu.formulatext
  
  tau.formulatext <- gsub("left_parenthesis", "(", tau.formula, fixed = TRUE)
  tau.formulatext <- gsub("right_parenthesis", ")", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("tilde_symbol", "~", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("equal_symbol", "=", tau.formulatext, fixed = TRUE)
  tau.formulatext <- gsub("comma_symbol", ",", tau.formulatext, fixed = TRUE)
  tau.formula <- stats::as.formula(tau.formulatext)
  tau.formula2use <- stats::as.formula(paste0(Reduce(paste, deparse(tau.formula))), env = parent.frame()) # here we need the formula as a 'call' object
  tau.original.formulatext <- tau.formulatext
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family)))
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  ## Get the variable names
  # Rewrite formula extracting variables nested in strutures like data frame or list
  # (e.g. D$A~D$B will be re-written A~B)
  # Note final product is a list of the variables in the model (yvector and covariates)
  # it is NOT a list of model terms - these are derived later
  
  # Convert formula string into separate variable names split by |
  formulatext <- gsub("pb(", "", formulatext, fixed=TRUE) 
  formulatext <- gsub(")", "", formulatext, fixed=TRUE) 
  formulatext <- gsub(" ", "", formulatext, fixed=TRUE)
  formulatext <- gsub("~", "|", formulatext, fixed=TRUE)
  formulatext <- gsub("+", "|", formulatext, fixed=TRUE)
  formulatext <- gsub("*", "|", formulatext, fixed=TRUE)
  formulatext <- gsub("||", "|", formulatext, fixed=TRUE)
  
  sigma.formulatext <- gsub("pb(", "", sigma.formulatext, fixed=TRUE) 
  sigma.formulatext <- gsub(")", "", sigma.formulatext, fixed=TRUE) 
  sigma.formulatext <- gsub(" ", "", sigma.formulatext, fixed=TRUE)
  sigma.formulatext <- gsub("~", "|", sigma.formulatext, fixed=TRUE)
  sigma.formulatext <- gsub("+", "|", sigma.formulatext, fixed=TRUE)
  sigma.formulatext <- gsub("*", "|", sigma.formulatext, fixed=TRUE)
  sigma.formulatext <- gsub("||", "|", sigma.formulatext, fixed=TRUE)
  
  nu.formulatext <- gsub("pb(", "", nu.formulatext, fixed=TRUE) 
  nu.formulatext <- gsub(")", "", nu.formulatext, fixed=TRUE) 
  nu.formulatext <- gsub(" ", "", nu.formulatext, fixed=TRUE)
  nu.formulatext <- gsub("~", "|", nu.formulatext, fixed=TRUE)
  nu.formulatext <- gsub("+", "|", nu.formulatext, fixed=TRUE)
  nu.formulatext <- gsub("*", "|", nu.formulatext, fixed=TRUE)
  nu.formulatext <- gsub("||", "|", nu.formulatext, fixed=TRUE)
  
  tau.formulatext <- gsub("pb(", "", tau.formulatext, fixed=TRUE) 
  tau.formulatext <- gsub(")", "", tau.formulatext, fixed=TRUE) 
  tau.formulatext <- gsub(" ", "", tau.formulatext, fixed=TRUE)
  tau.formulatext <- gsub("~", "|", tau.formulatext, fixed=TRUE)
  tau.formulatext <- gsub("+", "|", tau.formulatext, fixed=TRUE)
  tau.formulatext <- gsub("*", "|", tau.formulatext, fixed=TRUE)
  tau.formulatext <- gsub("||", "|", tau.formulatext, fixed=TRUE)
  
  # Remember model.variables and then varnames include both yvect and linear predictor components 
  mu.model.variables <- unlist(strsplit(formulatext, split="|", fixed=TRUE))
  sigma.model.variables <- unlist(strsplit(sigma.formulatext, split="|", fixed=TRUE))
  nu.model.variables <- unlist(strsplit(nu.formulatext, split="|", fixed=TRUE))
  tau.model.variables <- unlist(strsplit(tau.formulatext, split="|", fixed=TRUE))
  
  # mu
  mu.varnames <- c()
  for(i in 1:length(mu.model.variables)){
    elt <- unlist(strsplit(mu.model.variables[i], split="$", fixed=TRUE))
    if(length(elt) > 1){
      assign(elt[length(elt)], eval(parse(text=mu.model.variables[i]), envir = parent.frame()), envir = parent.frame())
      modified.formulatext <- gsub(mu.model.variables[i], elt[length(elt)], original.formulatext, fixed=TRUE)
      mu.varnames <- append(mu.varnames, elt[length(elt)])
    }else{
      mu.varnames <- append(mu.varnames, elt)
    }
  }
  mu.varnames <- unique(mu.varnames)
  
  # sigma
  sigma.varnames <- c()
  for(i in 1:length(sigma.model.variables)){
    elt <- unlist(strsplit(sigma.model.variables[i], split="$", fixed=TRUE))
    if(length(elt) > 1){
      assign(elt[length(elt)], eval(parse(text=sigma.model.variables[i]), envir = parent.frame()), envir = parent.frame())
      sigma.modified.formulatext <- gsub(sigma.model.variables[i], elt[length(elt)], sigma.original.formulatext, fixed=TRUE)
      sigma.varnames <- append(sigma.varnames, elt[length(elt)])
    }else{
      sigma.varnames <- append(sigma.varnames, elt)
    }
  }
  sigma.varnames <- unique(sigma.varnames)
  
  # nu
  nu.varnames <- c()
  for(i in 1:length(nu.model.variables)){
    elt <- unlist(strsplit(nu.model.variables[i], split="$", fixed=TRUE))
    if(length(elt) > 1){
      assign(elt[length(elt)], eval(parse(text=nu.model.variables[i]), envir = parent.frame()), envir = parent.frame())
      nu.modified.formulatext <- gsub(nu.model.variables[i], elt[length(elt)], nu.original.formulatext, fixed=TRUE)
      nu.varnames <- append(nu.varnames, elt[length(elt)])
    }else{
      nu.varnames <- append(nu.varnames, elt)
    }
  }
  nu.varnames <- unique(nu.varnames)
  
  # tau
  tau.varnames <- c()
  for(i in 1:length(tau.model.variables)){
    elt <- unlist(strsplit(tau.model.variables[i], split="$", fixed=TRUE))
    if(length(elt) > 1){
      assign(elt[length(elt)], eval(parse(text=tau.model.variables[i]), envir = parent.frame()), envir = parent.frame())
      tau.modified.formulatext <- gsub(tau.model.variables[i], elt[length(elt)], tau.original.formulatext, fixed=TRUE)
      tau.varnames <- append(tau.varnames, elt[length(elt)])
    }else{
      tau.varnames <- append(tau.varnames, elt)
    }
  }
  tau.varnames <- unique(tau.varnames)
  
} 
# AGGREGATE FUNCTION
# gamlssDS2


