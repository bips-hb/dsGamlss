#'
#' @title gamlssDS1 an assign function called by ds.galmss
#' @description This function calls the gamlssDS1 that sets up the model frame 
#' on the server side and performs some checks on the data.
#' @details For additional details please see the extensive header of ds.gamlss and also
#' the gamlss function in native R gamlss package.
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
#' @param data a data frame containing the variables occurring in the formula. 
#' If this is missing, the variables should be on the parent environment.
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
                     tau.formula = tau.formula, family = family, data=data, mu.fix = mu.fix,  
                     sigma.fix = sigma.fix, nu.fix = nu.fix, tau.fix = tau.fix){
  
  
  #**************************************************************************
  #I) Preparation ----  
  #**************************************************************************
  
  ## Capture the nfilter settings
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  
  data <- eval(parse(text = data), envir = parent.frame())
  
  ## Reconvert the special symbols to create the appropriate formula & gamlss.family objects
  formula <- gsub("left_parenthesis", "(", formula, fixed = TRUE)
  formula <- gsub("right_parenthesis", ")", formula, fixed = TRUE)
  formula <- gsub("tilde_symbol", "~", formula, fixed = TRUE)
  formula <- gsub("equal_symbol", "=", formula, fixed = TRUE)
  formula <- gsub("comma_symbol", ",", formula, fixed = TRUE)
  formula <- stats::as.formula(formula)
  
  sigma.formula <- gsub("left_parenthesis", "(", sigma.formula, fixed = TRUE)
  sigma.formula <- gsub("right_parenthesis", ")", sigma.formula, fixed = TRUE)
  sigma.formula <- gsub("tilde_symbol", "~", sigma.formula, fixed = TRUE)
  sigma.formula <- gsub("equal_symbol", "=", sigma.formula, fixed = TRUE)
  sigma.formula <- gsub("comma_symbol", ",", sigma.formula, fixed = TRUE)
  sigma.formula <- stats::as.formula(sigma.formula)
  
  nu.formula <- gsub("left_parenthesis", "(", nu.formula, fixed = TRUE)
  nu.formula <- gsub("right_parenthesis", ")", nu.formula, fixed = TRUE)
  nu.formula <- gsub("tilde_symbol", "~", nu.formula, fixed = TRUE)
  nu.formula <- gsub("equal_symbol", "=", nu.formula, fixed = TRUE)
  nu.formula <- gsub("comma_symbol", ",", nu.formula, fixed = TRUE)
  nu.formula <- stats::as.formula(nu.formula)
  
  tau.formula <- gsub("left_parenthesis", "(", tau.formula, fixed = TRUE)
  tau.formula <- gsub("right_parenthesis", ")", tau.formula, fixed = TRUE)
  tau.formula <- gsub("tilde_symbol", "~", tau.formula, fixed = TRUE)
  tau.formula <- gsub("equal_symbol", "=", tau.formula, fixed = TRUE)
  tau.formula <- gsub("comma_symbol", ",", tau.formula, fixed = TRUE)
  tau.formula <- stats::as.formula(tau.formula)
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family)))
  
  #**************************************************************************
  #I) Mu: Create modelframe----  
  #**************************************************************************
  
  ## Save call to transfer parameters to model.frame function
  gamlsscall <- match.call()  # the function call with arguments specified by full names
  
  ## Check for NA in the data 
  if(!missing(data)){
   if (any(is.na(data))){
     stop("The data contains NA's, use data = na.omit(mydata)")
   }  
  }  
  
  ## Evaluate the model frame for mu
  mnames <- c("", "formula", "data")  # relevant names that should be extracted from the calls
  cnames <- names(gamlsscall)  # get the names of the arguments of the call (first element "")
  cnames <- cnames[match(mnames,cnames,0)]  # keep only the ones that match with mnames
  mcall <- gamlsscall[cnames]  # get in mcall all the relevant information but remember that the first element
  # (function name) will be NULL
  mcall[[1]] <- as.name("model.frame")  # replace NULL with model.frame (to be able to execute model.frame 
  # function later)
  
  ## Specials for smoothing
  # add specials attribute for smoothing to formula object
  mcall$formula <- terms(formula, specials=.gamlss.sm.list, data=data)
  
  mu.frame <- eval(mcall, sys.parent())  # calls the model.frame function inside mcall to create the 
  # modelframe with the variables needed to use formula
  # also uses pb() function to create model frame for smoothing
  
  ## This part deals with the family 
  family <- as.gamlss.family(family)  # bring first the gamlss family
  G.dev.expr <- body(family$G.dev.inc)  # expression to calculate deviance increment for family 
  
}
