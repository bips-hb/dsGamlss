#'
#' @title gamlssDS7 called by ds.gamlss
#' @description This is the seventh serverside aggregate function called by ds.gamlss.
#' @details It is an aggregation function that calculates the normalized quantile
#' residuals. Furthermore, variables that were saved during the computation are deleted 
#' from the server. For more details please see the extensive header of ds.gamlss and 
#' also the gamlss function in native R gamlss package.
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
#' @param k the number of the nearest neighbours for which the mean is calculated
#' in the anonymization procedure.
#' @return anonymized normalized quantile residuals for the gamlss model fitted with 
#' ds.gamlss. 
#' @author Annika Swenne
#' @import gamlss.dist
#' @export

gamlssDS7 <- function(family = family, data=data, mu.fix=mu.fix,
                      sigma.fix = sigma.fix, nu.fix = nu.fix, tau.fix = tau.fix,
                      control = control, i.control = i.control, k = k){
  
  #**************************************************************************
  # I) Preparation ---- 
  #**************************************************************************
  
  ## Capture the nfilter settings
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.kNN <- as.numeric(thr$nfilter.kNN)   
  
  ## Replicate rqres within gamlss enviroment 
  rqres <- function (pfun="pNO", type=c("Continuous", "Discrete", "Mixed"), censored=NULL,  
                     ymin=NULL, mass.p=NULL, prob.mp=NULL, y=y, ... ){ 
    # function to calculate the normalized (randomized) quantile residuals of the gamlss object
  }
  body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
  
  ## Get the value of the 'data' parameter provided as character on the client side
  dataname <- data
  if(is.null(dataname)){
    data <- NULL 
  }else{
    data <- eval(parse(text=dataname), envir=parent.frame())
  }
  
  family <- gsub("left_parenthesis", "(", family, fixed = TRUE)
  family <- gsub("right_parenthesis", ")", family, fixed = TRUE)
  family <- gsub("equal_symbol", "=", family, fixed = TRUE)
  family <- gsub("comma_symbol", ",", family, fixed = TRUE)
  family <- gamlss.dist::as.family(eval(parse(text=family), envir=environment()))
  
  c1 <- as.numeric(unlist(strsplit(control, split=",")))
  mu.step <- c1[3]
  sigma.step <- c1[4]
  nu.step <- c1[5]
  tau.step <- c1[6]
  c2 <- as.numeric(unlist(strsplit(i.control, split=",")))
  
  #**************************************************************************
  # II) Get the outcome ----  
  #**************************************************************************
  
  #*A) Get fitted model ----
  mod.gamlss.ds <- base::get("temp_mod.gamlss.ds", env=parent.frame())
  
  ## get outcome
  y <- as.vector(mod.gamlss.ds$y)
  
  ## Get fitted values for all distribution parameters
  # necessary for all parameters to calculate deviance
  if("mu" %in% names(family$parameters)){
    mu <- base::get("temp_mu", env=parent.frame())
  }
  if("sigma" %in% names(family$parameters)){
    sigma <- base::get("temp_sigma", env=parent.frame())
  }
  if("nu" %in% names(family$parameters)){
    nu <- base::get("temp_nu", env=parent.frame())
  }
  if("tau" %in% names(family$parameters)){
    tau <- base::get("temp_tau", env=parent.frame())
  }
  
  #**************************************************************************
  # III) Calculate residuals ----  
  # apply k-nearest neighbour anonymisation to residuals as used in
  # scatterPlotDS (method.indicator=1)
  #**************************************************************************
  
  residuals <- eval(family$rqres)
  residuals <- stats::na.omit(residuals) 
  nresid <- length(residuals)
  
  # standardise the residuals (maybe not necessary)
  # residuals approximate standard normal distribution if model correct
  residuals.standardised <- (residuals-mean(residuals))/stats::sd(residuals)
  
  # Check if k is integer and has a value greater than or equal to the pre-specified threshold
  # and less than or equal to the length of rows of data.complete minus the pre-specified threshold
  if(k < nfilter.kNN | k > (nresid - nfilter.kNN)){
    stop(paste0("k must be greater than or equal to ", nfilter.kNN, "and less than or equal to ", (nresid-nfilter.kNN), "."), call.=FALSE)
  }else{
    neighbours <- k
  }
  
  # Find the k-1 nearest neighbours of each data point
  nearest <- RANN::nn2(residuals.standardised, k=neighbours)
  
  # Calculate the centroid of each n nearest data points
  residuals.centroid <- matrix()
  for (i in 1:nresid){
    residuals.centroid[i] <- mean(residuals.standardised[nearest$nn.idx[i,1:neighbours]])
  }
  
  # Calculate the scaling factor
  residuals.scalingFactor <- stats::sd(residuals.standardised)/stats::sd(residuals.centroid)
  
  # Apply the scaling factor to the centroids
  residuals.masked <- residuals.centroid * residuals.scalingFactor
  
  # Shift the centroids back to the actual position and scale of the original data
  residuals.new <- (residuals.masked * stats::sd(residuals)) + mean(residuals)
  
  #**************************************************************************
  # IV) Delete variables ----  
  # delete all variables that were used during the computation and stored
  # in the parent.frame() environment
  #**************************************************************************
  
  # delete temporary variables (prefix temp_) will be deleted
  rm(list=ls(pattern="^temp_", envir=parent.frame()), envir=parent.frame())
  
  return(residuals.new)
}