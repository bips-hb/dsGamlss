#' 
#' @title Local function that creates the basis matrix for p-splines
#' @description This is an internal function called by the functions 'gamlssDS3',
#' 'gamlssDS4', 'gamlssDS5' and 'gamlssDS6'.
#' @details The function is based on the internal bbase-function defined in gamlss::pb (version 5.4-12).
#'          Note that unlike in the original function the use of quantiles is currently not implemented in
#'          bbase
#' @param x A vector of x-values.
#' @param xl Left limit for the knots.
#' @param xr Right limit for the knots.
#' @param ndx Number of equal space intervals in x.
#' @param deg Degree of the polynomial.
#' @return A B-basis matrix for p-splines.
#' @keywords internal
#' @author Annika Swenne

bbase <- function(x, xl, xr, ndx=20, deg=3){
  
  #**************************************************************************
  #I) Define functions----  
  #**************************************************************************
  tpower <- function(x, t, p){
    # Truncated p-th power function (defined for single values)
    # used to construct B-splines
    # *************************************************************************
    # x : single x-value
    # t : single knot
    # p : power
    (x-t)^p*(x>t)  # x>t indicator function TRUE=1, FALSE=0
  }
  
  #**************************************************************************
  #II) Construct B-spline basis----  
  #************************************************************************** 
  dx <- (xr-xl)/ndx # increment to ensure the desired number of intervals ndx 
  
  ## Calculate B-splines for equally spaced knots via truncated power function as described in Eilers & Marx (2010)
  # add deg extra knots below xl & above xr for truncated power function
  knots <- seq(xl-deg*dx, xr+deg*dx, by=dx)
  # calculate the power in the knots (evaluate function tpower for all combinations of x-values and knots)
  P <- outer(x, knots, tpower, deg)  # matrix of dimension length(x)*length(knots)
  n <- dim(P)[2]
  D <- diff(diag(n), diff=deg+1)/(gamma(deg+1)*dx^deg) 
    # diff(diag(n), diff=deg+1) difference matrix of order deg+1
    # gamma(deg+1) returns gamma function of deg+1 (single value)
  B <- (-1)^(deg+1)*P %*% t(D)
  attr(B, "knots") <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
  B 
}



