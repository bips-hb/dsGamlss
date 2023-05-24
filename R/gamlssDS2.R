gamlssDS <- function(formula = formula,
                   sigma.formula = sigma.formula,
                   nu.formula = nu.formula,
                   tau.formula = tau.formual,
                   family = family,
                   data=data,
                   mu.fix = mu.fix,  
                   sigma.fix = sigma.fix,
                   nu.fix = nu.fix,
                   tau.fix = tau.fix,
                   control = control, 
                   i.control = i.control){
  ##***************************************************************************************
  ## require(stats) Thursday, June 10, 2004 at 09:58 MS
  #require(splines) # this will be removed with namespaces
  #****************************************************************************************
  #gamlss.rc.list<-c("EX.rc","Exponential.rc") # the right censoring distribution list 
  #gamlss.bi.list<-c("BI", "Binomial", "BB", "Beta Binomial") # binomial denominators
  #.gamlss.multin.list<-c("MULTIN", "MN3", "MN4", "MN5")
  #****************************************************************************************
  # to fit a gamlss model
  #****************************************************************************************
  #       formula : formula for the mu parameter
  # sigma.formula : formula for the sigma parameter (does not require response variable)
  #    nu.formula : formula for the sigma parameter (does not require response variable)
  #   tau.formula : formula for the tau parameter (does not require response variable)
  #        family : distribution of the response variable (defines link functions, derivatives, etc.)
  #          data : data frame containing variables in the formulas
  #        mu.fix : whether the parameter mu is fixed 
  #     sigma.fix : whether the parameter sigma is fixed
  #        nu.fix : whether the parameter nu is fixed 
  #       tau.fix : whether the parameter tau is fixed 
  #       control : control parameters for outer iteration 
  #     i.control : control parameters for inner iteration & backfitting
  
  #**************************************************************************
  #I) Define functions----  
  #**************************************************************************
  
  # this is to replicate rqres within gamlss enviroment DS Friday, March 31, 2006 at 10:30
  rqres <- function (pfun="pNO", type=c("Continuous", "Discrete", "Mixed"), censored=NULL,  
                     ymin=NULL, mass.p=NULL, prob.mp=NULL, y=y, ... ){ 
    # function to calculate the normalized (randomized) quantile residuals of the gamlss object
  }
  body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
  
  RS <- function(n.cyc=control$n.cyc, no.warn=TRUE){
    # executes the outer iteration of the gamlss algorithm
    #**************************************************************************
    #   n.cyc : maximum number of cycles for the outer iteration
    # no.warn : whether a warning should be printed if the algorithm has not 
    #           converged after n.cyc iterations (set to FALSE inside mixed())
    
    #**************************************************************************
    #1) Define function----  
    #**************************************************************************
    glim.fit <- function(f, X, y, w, fv, os, step = 1, control = glim.control(), 
                         auto, gd.tol){ 
      # emulates the GLIM iterative algorithm created by Mikis Stasinopoulos 
      # on Monday, February 18, 2002 at 09:23
      # last change:
      # Wednesday, February 23, 2005 at 20:23 where the step has changed to apply 
      # to the fitted values 
      #**************************************************************************
      #       f : parameter.object (e.g. mu.object)
      #       X : design matrix for parameter (e.g. mu.X)
      #       y : response vector
      #       w : vector of weights for weighted likelihood analysis (not the same as in GLM's) 
      #      fv : current parameter estimate (e.g. mu)
      #      os : offset for parameter in linear predictor
      #    step : step length for the parameter (e.g. mu.step)
      # control : control parameters for inner iteration & backfitting
      #  gd.tol : global deviance tolerance level by which the deviance is allowed to increase in 
      #           each iteration 
      #    auto : whether the steps should be halved automatically if the new global deviance is 
      #           greater than the old one
      
      #**************************************************************************
      #i) Initialization ----  
      #**************************************************************************
      cc <- control$cc  # convergence criterion for inner iteration
      cyc <- control$cyc  # maximum number of cycles for inner iteration
      trace <- control$glm.trace  # whether to print at each inner iteration (TRUE) or not (FALSE)
      bf.cyc <- control$bf.cyc  # maximum number of cycles for backfitting algorithm
      bf.tol <- control$bf.tol  # convergence criterion for backfitting algorithm
      bf.trace <- control$bf.trace  # whether to print at each backfitting iteration (TRUE) or not (FALSE)
      itn <- 0  # number of inner iteration
      
      ## Calculation of weights (for Fisher-scoring algorithm)
      lp <- eta <- f$linkfun(fv)  # predictor (fv=mu)
      dr <- f$dr(eta)  # dmu/deta
      dr <- 1/dr  # deta/dmu = 1 / (dmu/deta)
      di <- f$G.di(fv)  # deviance increment
      dv <- sum(w*di)  # the global deviance
      olddv <- dv+1  # the old global deviance (only for the first step when old global deviance does not
      # exist to ensure that algorithm continues, i.e. abs(olddv-dv) > cc)
      dldp <- f$dldp(fv)  # first derivative of log-likelihood with respect to mu
      d2ldp2 <- f$d2ldp2(fv)  # expected  second derivative of log-Likelihood with respect to mu
      d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)  # added 26-10-07  
      wt <- -(d2ldp2/(dr*dr)) # weights for Fisher-scoring algorithm =-(d2l/d2p)(dmu/deta)^2
      # we need to stop the weights to go to Infty
      wt <- ifelse(wt>1e+10,1e+10,wt) # Mikis 9-10-14 
      wt <- ifelse(wt<1e-10,1e-10,wt) 
      wv <- (eta-os)+dldp/(dr*wt) # working variable z 
      if (family$type=="Mixed"){
        wv <-ifelse(is.nan(wv),0,wv) ## TEST
      }
      iterw <- FALSE  # indicator whether deviance has increased in an inner iteration
      # for smoothing
      who <- f$who  # formula for smoothing
      smooth.frame <- f$smooth.frame  # data frame for smoother
      s <- f$smooth  # matrix for fitted smoothers
      
      #**************************************************************************
      #ii) Inner Iteration  ----     
      #**************************************************************************
      while (abs(olddv-dv) > cc && itn < cyc){  # MS Wednesday, June 26, 2002 
        itn <- itn+1  # the glim inner iteration number
        lpold <- lp  # old eta
        sold <- s  # old matrix with fitted smoothers
        
        if (any(is.na(wt))||any(is.na(wv))){
          stop("NA's in the working vector or weights for parameter ", names(formals(f$valid)), "\n")
        }
        if (any(!is.finite(wt))||any(!is.finite(wv))){
          stop("Inf values in the working vector or weights for parameter ", names(formals(f$valid)), "\n")
        }
        
        #*a) Backfitting  ----  
        if(length(who) > 0){  
          #*a.1) with smoothing terms
          fit <- additive.fit(x=X, y=wv, w=wt*w, s=s, who=who, smooth.frame, maxit=bf.cyc, 
                              tol=bf.tol, trace=bf.trace)
          # avoid overjumping via step parameter (method 1 as described in Stasinopolous et al. 2020, p.66f)
          # for method 2 (autostep=TRUE) it holds 1-step=0 and the autostep method is implemented in part b)
          if (itn==1){
            lp <- fit$fitted.values  # new eta = fitted values for working variable (linear & nonlinear part)
            s <- fit$smooth  # new s = fitted values for smoothers (nonlinear part)
          } else {
            lp <- step*fit$fitted.values+(1-step)*lpold
            s <- step*fit$smooth+(1-step)*sold
          }
        } else {  
          #*a.2) without smoothing terms
          fit <- lm.wfit(x=X, y=wv, w=wt*w, method="qr") 
          if (itn==1){
            lp <- fit$fitted.values  # new eta = fitted values for working variable (linear & nonlinear part)
          } else {
            lp <- step*fit$fitted.values+(1-step)*lpold
          }
        }
        
        eta <- lp+os  # add offset to new eta
        fv <- f$linkinv(eta)  # new mu
        olddv <- dv
        di <- f$G.di(fv)  # deviance increment
        dv <- sum(w*di)  # new global deviance 
        
        #*b) Autostep ----
        # to avoid overjumping (method 2 as described in Stasinopolous et al. 2020, p.66f)
        if (dv > olddv && itn >= 2 && auto==TRUE){
          for(i in 1:5){ # MS Thursday, September 22, 2005 
            lp <- (lp+lpold)/2  # weighted eta
            eta <- lp+os  # add offset to new eta
            fv <- f$linkinv(eta)  # mu
            di <- f$G.di(fv)   # deviance increment
            dv <- sum(w*di)  # new global deviance
            #  cat("try",i,"\n")
            if(length(who) > 0){
              s <- (s+sold)/2  # uses the same weights for the smoothing matrix
            }
            if ((olddv-dv) > cc){
              break # new deviance smaller than old one
            }
          }
        }
        if ((dv > olddv+gd.tol ) && itn >= 2 && iterw==FALSE){
          warning("The deviance has increased in an inner iteration for ",
                  names(formals(f$valid)), "\n","Increase gd.tol and if persist, try different steps",  "\n", "or model maybe inappropriate")
          iterw <-TRUE
        }                      
        if (is.na(!f$valid(fv))){ # MS Saturday, April 6, 2002 at 18:06 
          stop( "fitted values in the inner iteration out of range")
        } 
        
        #*c) Calculate weights for new eta ----
        dr <- f$dr(eta)   # dmu/deta
        dr <- 1/dr  # deta/dmu = 1 / (dmu/deta)
        dldp <- f$dldp(fv)  # first derivative of log-likelihood with respect to mu  
        d2ldp2 <- f$d2ldp2(fv)  # expected  second derivative of log-Likelihood with respect to mu
        d2ldp2 <-  ifelse(d2ldp2 < -1e-15, d2ldp2,-1e-15) # added 26-10-07   
        wt <- -(d2ldp2/(dr*dr))  # weights for Fisher-scoring algorithm =-(d2l/d2p)(dmu/deta)^2
        wt <- ifelse(wt>1e+10,1e+10,wt) # Mikis 9-10-14    
        wt <- ifelse(wt<1e-10,1e-10,wt) 
        wv <- (eta-os)+dldp/(dr*wt)  # working variable z 
        if (family$type=="Mixed"){
          wv <-ifelse(is.nan(wv),0,wv) ## TEST  
        }
        if (trace){
          cat("GLIM iteration ", itn, " for ",  names(formals(f$valid)), ": Global Deviance = ",
              format(round(dv, 4)), " \n", sep = "")
        }
      } #end of inner iteration
      
      #**************************************************************************
      #iii) Return estimates  ----
      #**************************************************************************
      pen <- 0  # sum of the quadratic penalties for parameter
      if(length(who) > 0){
        pen <- sum(eta*wt*(wv-eta))
      }
      c(fit, list(fv=fv, wv=wv, wt=wt, eta=eta, os=os, pen=pen)) #ms Saturday, December 4, 2004 
    } # end of glim.fit function   
    
    #**************************************************************************
    #2) Initialization----  
    #**************************************************************************
    ## getting the control parameters from control object
    c.crit <- control$c.crit  # convergence criterion for outer iteration
    n.cyc <- control$n.cyc  # maximum number of cycles for outer iteration
    trace <- control$trace  # whether to print at each outer iteration (TRUE) or not (FALSE)
    autostep <- control$autostep  # default autostep=TRUE, i.e. steps are halved automatically 
    # up to 5 times if the new global deviance is greater than the 
    # old one
    mu.step <- control$mu.step  # step length for parameter mu (default mu.step=1)
    sigma.step <- control$sigma.step  # step length for parameter sigma (default sigma.step=1)
    nu.step <- control$nu.step  # step length for parameter nu (default nu.step=1)
    tau.step <- control$tau.step  # step length for parameter tau (default tau.step=1)
    gd.tol <- control$gd.tol  # global deviance tolerance level by which the deviance is allowed to increase in 
    # each iteration (default gd.tol=Inf)
    iter <- control$iter  # starting value for the number of outer iterations (typically 0 unless refit is used)
    conv <- FALSE  # the algorithm has not converged yet
    
    ## initial Gloval deviance
    G.dev.incr  <- eval(G.dev.expr)  # deviance increment (function provided by gamlss.family object)
    G.dev <- sum(w*G.dev.incr)  # the weighted global deviance
    G.dev.old <- G.dev+1  # the old global deviance (only for the first step when old global deviance does not
    # exist to ensure that algorithm continues, i.e. abs(G.dev.old-G.dev) > c.crit)
    
    #**************************************************************************
    #3) Outer Iteration----  
    #**************************************************************************
    while (abs(G.dev.old-G.dev) > c.crit && iter < n.cyc){
      
      #*a) Mu: Inner Iteration ----
      if ("mu" %in% names(family$parameters)){
        if (family$parameter$mu==TRUE & mu.fix==FALSE){
          mu.fit <<- glim.fit(f=mu.object, X=mu.X, y=y, w=w,
                              fv=mu, os=mu.offset, step=mu.step,
                              control=i.control, gd.tol=gd.tol,
                              auto=autostep)
          mu <<- mu.fit$fv  # update parameter estimate
          mu.object$smooth <- mu.fit$smooth  # fitted values for smoothers (nonlinear part)
        }
      }
      
      #*b) Sigma: Inner Iteration ----
      if ("sigma" %in% names(family$parameters)){
        if (family$parameter$sigma==TRUE & sigma.fix==FALSE){  
          sigma.fit  <<- glim.fit(f=sigma.object, X=sigma.X, y=y, 
                                  w=w, fv=sigma, os=sigma.offset, 
                                  step=sigma.step, control=i.control,
                                  gd.tol=gd.tol, auto=autostep)
          sigma <<- sigma.fit$fv  # update parameter estimate
          sigma.object$smooth <- sigma.fit$smooth  # fitted values for smoothers (nonlinear part)
        }
      }
      
      #*c) Nu: Inner Iteration ----
      if ("nu" %in% names(family$parameters)){
        if (family$parameter$nu==TRUE & nu.fix==FALSE){
          nu.fit <<- glim.fit(f=nu.object, X=nu.X, y=y,
                              w=w, fv=nu, os=nu.offset, 
                              step=nu.step, control=i.control, 
                              gd.tol=gd.tol, auto=autostep)
          nu <<- nu.fit$fv  # update parameter estimate
          nu.object$smooth <- nu.fit$smooth  # fitted values for smoothers (nonlinear part)
        }
      }
      
      #*d) Tau: Inner Iteration ----
      if ("tau" %in% names(family$parameters)){
        if (family$parameter$tau ==TRUE & tau.fix==FALSE){ 
          tau.fit <<- glim.fit(f=tau.object, X=tau.X, y=y,
                               w=w, fv=tau, os=tau.offset,
                               step=tau.step, control=i.control, 
                               gd.tol=gd.tol, auto=autostep)
          tau <<- tau.fit$fv  # update parameter estimate
          tau.object$smooth <- tau.fit$smooth  # fitted values for smoothers (nonlinear part)
        }
      }
      
      #*e) Check deviance ----
      G.dev.old <- G.dev
      G.dev.incr  <- eval(G.dev.expr)  # deviance increment (function provided by gamlss.family object)
      G.dev <- sum(w*G.dev.incr)  # the weighted global deviance
      
      iter <- iter+1  # number of outer iterations
      fiter <<- iter  
      if(trace){
        cat("GAMLSS-RS iteration ", iter, ": Global Deviance = ",
            format(round(G.dev, 4)), " \n", sep = "")
      }
      if (G.dev > (G.dev.old+gd.tol) && iter >1 ){
        stop(paste("The global deviance is increasing", "\n", 
                   "Try different steps for the parameters or the model maybe inappropriate"))
      }
    } # end of while loop
    
    # check whether the outer iteration algorithm converged
    if (abs(G.dev.old-G.dev) < c.crit){ # MS Wednesday, June 11, 2003 at 11:58 
      #taken out (abs((G.dev-G.dev.old)/(0.1+abs(G.dev.old)))<c.crit&&iter<=n.cyc) 
      conv <- TRUE
    } else {
      conv <- FALSE
    }
    if (!conv && no.warn ){
      warning("Algorithm RS has not yet converged")
    }
    conv
    
  } # end of GAMLSS RS fitting algorithm
  
  get.smoothers <- function(term){
    # gets the smoothers from a terms object
    # *************************************************************************
    # terms : terms object from which the smoothers should be extracted
    
    a <- attributes(term)
    smoothers <- a$specials  # extract specials attribute that includes the smoothers
    if (length(smoothers) > 0){
      # identify smoothers that are non-zero
      smoothers <- smoothers[sapply(smoothers, length) > 0]
      # iterate over all elements in smoothers
      for (i in seq(along = smoothers)){ 
        tt <- smoothers[[i]]
        # apply over columns of a$factors 
        ff <- apply(a$factors[tt, ,drop=FALSE], 2, any)  # drop=FALSE to maintain array structure
        if(any(ff)){
          smoothers[[i]] <- seq(along = ff)[a$order==1 & ff]
        } else {
          smoothers[[i]] <- NULL
        }
      }
    }
    smoothers
  }
  
  get.object <- function(what){
    # create parameter object, i.e. list of certain attributes (inherits components from family object)
    # *************************************************************************
    # what : string with the name of the parameter ("mu", "sigma", "nu", "tau")
    
    # extract attributes from family object
    link <- eval(parse(text=(paste("family$", what, ".link", sep=""))))  # name of link function
    linkfun <- eval(parse(text=(paste("family$", what, ".linkfun", sep=""))))  # link function 
    linkinv <- eval(parse(text=(paste("family$", what, ".linkinv", sep=""))))  # inverse link function   
    dr <- eval(parse(text=(paste("family$",what,".dr", sep=""))))  # dparameter/deta
    dldp <- switch(what,
                   "mu" = family$dldm,
                   "sigma" = family$dldd,
                   "nu" = family$dldv,
                   "tau" = family$dldt)  # first derivative of log-likelihood 
    d2ldp2 <- switch(what,
                     "mu" = family$d2ldm2,
                     "sigma" = family$d2ldd2,
                     "nu" = family$d2ldv2,
                     "tau" = family$d2ldt2)  # second derivative of log-likelihood   
    G.di <- family$G.dev.incr  # deviance function for single observation           
    valid <- eval(parse(text=(paste("family$",what,".valid", sep=""))))  # function specifying valid 
    # values for the parameter
    
    object <- list(link = link, linkfun = linkfun, linkinv = linkinv, dr = dr,
                   dldp = dldp,  d2ldp2 = d2ldp2, G.di = G.di, valid = valid)  
    
    # additional attributes for smoothing
    if(length(eval(parse(text=(paste(what,".smoothers", sep=""))))) > 0){  # if smoothing 
      parAttrTermlevels <- eval(parse(text=(paste(what,".a$term.labels", sep=""))))  # formulas for predictors
      boo <- unlist(eval(parse(text=(paste(what,".smoothers", sep="")))))  # which kind of smoothers (e.g. pb)
      who <- parAttrTermlevels[boo[order(boo)]]  # formulas for smoothing
      smooth.frame <- eval(parse(text=(paste(what,".frame", sep="")))) 
      s <- matrix(0, N, length(who))
      dimnames(s) <- list(names(y), who)
      object$smooth <- s  # smoothing matrix (dimension n* #smoothers)
      object$who <- who  # formulas for smoothing (as string)
      object$smooth.frame <- smooth.frame  # model frame for smoothing
    } 
    object   
  }
  
  other.formula <- function(form){
    # to add the response to the formula for other parameters than mu 
    # (they do not require a specification of the response variable)
    # *************************************************************************
    # form : formula as specified in the gamlss call
    dform <- formula(form)
    if (length(dform)==2){  # formula does not include response variable
      dform[3] <- dform[2]      # taking predictors in position [3]
      if (is(formula,"terms")){
        dform[2] <- formula[[2]]
      } else {
        dform[2] <- formula[2] # ms 31-12-08   # put y in position 2 
      }
    }
    dform 
  }
  
  parameterOut <- function(what="mu", save){
    # this function is used for outputting the results for a parameter
    # *************************************************************************
    # what : string with the name of the parameter ("mu", "sigma", "nu", "tau")
    # save : whether all information should be saved on exist. If save=FALSE only a reduced
    #        output is saved
    out <- list()
    if(save==TRUE){  
      ## all information should be saved
      if(family$parameter[[what]]==TRUE && eval(parse(text=paste(what,".fix",sep="")))==FALSE){
        ## parameter exists for the family & is not fixed
        out$fv <- eval(parse(text=what))  # fitted values for the parameter
        out$lp <- eval(parse(text=(paste(what,".fit$eta", sep=""))))  # estimate for eta
        out$wv <- eval(parse(text=(paste(what,".fit$wv", sep=""))))  # working variable z 
        out$wt <- eval(parse(text=(paste(what,".fit$wt", sep=""))))  # iterative weights
        out$link <- eval(parse(text=(paste(what,".object$link", sep=""))))  # link function (as character)
        out$terms <- eval(parse(text=(paste(what,".terms", sep=""))))  # terms object
        out$x <- eval(parse(text=(paste(what,".X", sep=""))))  # design matrix
        out$qr <- eval(parse(text=(paste(what,".fit$qr", sep=""))))  # qr decomposition for sqrt(w)*X
        out$coefficients <- eval(parse(text=(paste(what,".fit$coefficients", sep=""))))  # linear coefficients
        out$offset <- eval(parse(text=(paste(what,".fit$os", sep=""))))  # offset in the linear predictor
        out$xlevels  <- .getXlevels(eval(parse(text=paste(what,".terms",sep=""))), 
                                    eval(parse(text=paste(what,".frame",sep=""))))  # record of the factor levels
        out$formula <- eval(parse(text=paste(what,".formula",sep="")))  # formula for the parameter (as character)
        
        if(length(eval(parse(text=paste(what,".smoothers",sep="")))) > 0){  
          ## the model includes smoothers
          out$df <- eval(parse(text=paste(what,".fit$nl.df",sep="")))+  
            eval(parse(text=paste(what,".fit$rank",sep="")))  # degrees of freedom (linear + nonlinear)
          out$nl.df <- eval(parse(text=paste(what,".fit$nl.df",sep="")))  # nonlinear (e.g. smoothing) df
          # (does not include 2 df for fitted 
          # constant & linear part)
          out$s  <- eval(parse(text=paste(what,".fit$smooth",sep="")))  # fitted values for smoothers (nonlinear part)
          out$var <- eval(parse(text=paste(what,".fit$var",sep="")))  # variance for the smoothing fitted values
          out$coefSmo <- eval(parse(text=paste(what,".fit$coefSmo",sep="")))  # GAM object fitted within backfitting algorithm
          out$lambda <- eval(parse(text=paste(what,".fit$lambda",sep="")))  # smoothing parameter
          out$pen <- eval(parse(text=paste(what,".fit$pen",sep="")))  # sum of quadratic penalties         
        } else {  
          ## the model does not include smoothers
          out$df <- eval(parse(text=paste(what,".fit$rank",sep="")))  # degrees of freedom (linear + nonlinear)
          out$nl.df <- 0  # 0 nonlinear df (model does not include smoothers)
          out$pen  <- 0  # sum of quadratic penalties 
        }
      } else { 
        ## parameter does not exist for the family or is fixed
        out$fix <- eval(parse(text=paste(what,".fix",sep="")))  # whether the parameter is fixed
        out$df <- 0  # degrees of freedom for the parameter
        out$fv <- eval(parse(text=what))  # values for the parameter
      }
    } else {  
      ## if(save==FALSE)
      if(family$parameter[[what]]==TRUE && eval(parse(text=paste(what,".fix",sep="")))==FALSE){
        ## parameter exists for the family & is not fixed
        out$terms <- eval(parse(text=(paste(what,".terms", sep=""))))  # terms object
        out$formula <- eval(parse(text=paste(what,".formula",sep="")))  # formula for the parameter (as character)
        
        if(length(eval(parse(text=paste(what,".smoothers",sep="")))) > 0){
          ## the model includes smoothers
          out$df <- eval(parse(text=paste(what,".fit$nl.df",sep="")))+
            eval(parse(text=paste(what,".fit$rank",sep="")))  # degrees of freedom (linear + nonlinear)
          out$nl.df <- eval(parse(text=paste(what,".fit$nl.df",sep="")))  # nonlinear (e.g. smoothing) df
          # (does not include 2 df for fitted 
          # constant & linear part)
        } else {  
          ## the model does not include smoothers
          out$df <- eval(parse(text=paste(what,".fit$rank",sep="")))  # degrees of freedom (linear + nonlinear)
          out$nl.df <- 0  # 0 nonlinear df (model does not include smoothers)
        }
      } else {   
        ## parameter does not exist for the family or is fixed
        out$df <- 0  # degrees of freedom for the parameter
      }
    }
    out
  } # end of parameterOut function
  
  
  #**************************************************************************
  #III) Checks ----
  # does the data meet the requirements?
  #**************************************************************************
  
  #*a) Response variable ----
  ##*************************************************************************
  ## Extract the model components using model.extra and model.matrix
  Y <- model.extract(mu.frame, "response")  # extracting the y variable from the modelframe
  if(is.null(dim(Y))){  # if y not matrix
    N <- length(Y) 
  } else {
    N <- dim(Y)[1]   # number of observations in Y
  }
  
  ## Checking whether response matches specified gamlss.family
  if(any(family$family%in%.gamlss.bi.list)){  ## binomial checking
    # extracting now the y and the binomial denominator in case we use BI or BB
    if (NCOL(Y) == 1) {
      if (is.factor(Y)){  # convert factor to numeric using tthat boolean can be represented
        # as 0-1
        y <- Y != levels(Y)[1]
      } else {
        y <- Y
      }
      bd <- rep(1, N)
      if (any(y < 0 | y > 1)){
        stop("y values must be 0 <= y <= 1")
      } 
    } else if (NCOL(Y) == 2){
      if (any(abs(Y - round(Y)) > 0.001)) {
        warning("non-integer counts in a binomial GAMLSS!")
      }
      bd <- Y[,1] + Y[,2]  # total number of games
      y <-  Y[,1]
      if (any(y < 0 | y > bd)) {
        stop("y values must be 0 <= y <= N") # MS Monday, October 17, 2005
      }
    } else {
      stop(paste("For the binomial family, Y must be", 
                 "a vector of 0 and 1's or a 2 column", 
                 "matrix where col 1 is no. successes", 
                 "and col 2 is no. failures"))
    }
    
  } else if(any(family$family%in%.gamlss.multin.list)){  ## multinomial checking
    if(is.factor(Y)){  # convert factor to numbers
      y <- unclass(Y)
    } else {
      y <- Y
    }
    
  } else if(is.Surv(Y)){  ## survival checking
    # checking that the family is censored
    if (length(grep("censored",family$family[[2]]))==0) {
      stop(paste("the family in not a censored distribution, use cens()"))
    }
    # checking compatability of Surv object and censored distribution
    if (length(grep(attr(Y,"type"),family$family[[2]]))==0){
      stop(paste("the Surv object and the censored distribution are not of the same type"))
    }
    y <- Y
  } else {  ## other families
    y <- Y
  }
  
  ## checking the permissible y values      
  if (!family$y.valid(y)){ # MS Thursday, June 20, 2002 at 16:30
    stop( "response variable out of range")
  }  
  
  #*b) Other parameters ----
  ##*************************************************************************
  ## this part is used if start.from is used as argument
  if(!is.null(start.from)){
    if (!is.gamlss(start.from)){
      stop(paste("The object in start.from is not a gamlss object", "\n", ""))
    }
    mu.start <- NULL
    sigma.start <- NULL
    nu.start <- NULL
    tau.start <- NULL
    ## location model
    if ("mu"%in%start.from$parameters){
      mu.start <- start.from$mu.fv
    }
    ## scale-dispersion submodel
    if ("sigma"%in%start.from$parameters) {
      sigma.start <- start.from$sigma.fv
    }                
    ##  nu submodel
    if ("nu"%in%start.from$parameters){
      nu.start <- start.from$nu.fv
    } 
    ##  tau submodel
    if ("tau"%in%start.from$parameters){
      tau.start <- start.from$tau.fv 
    } 
  }
  
  ## checking that parameter.fix is boolean
  if (!is.logical(mu.fix)){
    stop("mu.fix should be logical TRUE or FALSE")
  } 
  if (!is.logical(sigma.fix)){
    stop("sigma.fix should be logical TRUE or FALSE")
  }
  if (!is.logical(nu.fix)){
    stop("nu.fix should be logical TRUE or FALSE")
  }
  if (!is.logical(tau.fix)){
    stop("tau.fix should be logical TRUE or FALSE")
  } 
  
  ## checking the weights
  w <- model.extract(mu.frame, weights) # weights for the likelihood
  if (is.null(w)){
    w <- rep(1, N)  # if no weights are specified they are set to 1
  } else if (any(w < 0)){
    stop("negative weights not allowed")
  } 
  
  ## checking whether proper algorithm  (RS, CG or mixed) 
  name.method <- substitute(method)
  name.method <- deparse(name.method[1])
  list.methods <- c("RS()","CG()","mixed()")
  i.method <- pmatch(name.method,list.methods,nomatch=0)
  if(!i.method){
    stop("Method must be RS(), CG() or mixed()")
  }
  
  #******************************************************************************************
  #IV) Set up submodels ----
  #******************************************************************************************
  ##  Set up location-mean submodel:
  ##             mu.X   design matrix
  ##        mu.offset   offset in linear predictor
  ##         mu.start   starting values for mu (optional)
  ##*****************************************************************************************
  mu.fit <- list()  # MS Thursday, January 23, 2003 at 14:46
  mu.formula <- formula  # ms Wednesday, December 29, 2004 
  mu.terms <- attr(mu.frame, "terms")  # get terms attribute from modelframe
  mu.smoothers <- get.smoothers(term=mu.terms)  
  mu.a <- attributes(mu.terms)  # get attributes from terms object
  mu.X <- model.matrix(mu.terms, mu.frame, contrasts)  # creates the design matrix for the mean
  mu.offset <- model.extract(mu.frame, offset)  # the mean-location offset
  if(is.null(mu.offset)){
    mu.offset <- rep(0,N)
  }
  mu.object <- get.object(what="mu")  # large list with different attributes (some inherited from family object)
  # the following is to get the right GLIM arguments 
  # replace function parameters ($y, $mu, $sigma, $nu, $tau) with (mu=fv)
  formals(mu.object$dldp, envir=new.env()) <- alist(mu = fv)
  formals(mu.object$d2ldp2, envir=new.env()) <- alist(mu = fv) 
  formals(mu.object$G.di, envir=new.env()) <- alist(mu = fv)
  formals(mu.object$valid, envir=new.env()) <- alist(mu = fv)
  # initial values for mu
  if(!is.null(mu.start)){ 
    if(length(mu.start)>1){
      mu <- mu.start
    } else {
      mu <- rep(mu.start,N)
    }
  } else {
    eval(family$mu.initial)  # use function in family object to calculate initial values
  }
  
  ##*****************************************************************************************
  ##  Set up dispersion-scale submodel:
  ##           sigma.X   design matrix
  ##      sigma.offset   offset in linear predictor
  ##       sigma.start   starting values for sigma (optional)
  ##*****************************************************************************************
  if ("sigma"%in%names(family$parameters)){
    orig.Envir  <- attr(mcall$formula, ".Environment")    # saves environment in which formula was defined
    sigma.fit <- list()    
    form.sigma <- other.formula(form=sigma.formula)  # to add response to formula
    
    ## Specials for smoothing
    # add specials attribute for smoothing to formula object
    if(missing(data)){
      sigma.terms <- terms(form.sigma, specials = .gamlss.sm.list)
    } else {
      sigma.terms <- terms(form.sigma, specials = .gamlss.sm.list, data = data)
    }
    mcall$formula <- sigma.terms
    attr(mcall$formula, ".Environment") <- orig.Envir # DS fix for Willem Thursday, March 18, 2010
    
    sigma.frame <- eval(mcall,sys.parent())  # calls the model.frame function inside mcall to create the 
    # modelframe with the variables needed to use formula
    # also uses pb() function to create model frame for smoothing
    sigma.terms <- attr(sigma.frame, "terms")  # get terms attribute from modelframe
    sigma.smoothers <- get.smoothers(sigma.terms)  
    sigma.a <- attributes(sigma.terms)  # get attributes from terms object     
    sigma.X <- model.matrix(sigma.terms, sigma.frame, contrasts)  # creates the design matrix for sigma
    sigma.offset <- model.extract(sigma.frame, offset)  # the sigma-location offset
    if(is.null(sigma.offset)){
      sigma.offset <- rep(0,N)
    }
    sigma.object <- get.object("sigma")  # large list with different attributes (some inherited from family object)
    # the following is to get the right GLIM arguments 
    # replace function parameters ($y, $mu, $sigma, $nu, $tau) with (sigma=fv)    
    formals(sigma.object$dldp, envir=new.env()) <- alist(sigma=fv)
    formals(sigma.object$d2ldp2, envir=new.env()) <- alist(sigma=fv)
    formals(sigma.object$G.di, envir=new.env()) <- alist(sigma=fv)
    formals(sigma.object$valid, envir=new.env()) <- alist(sigma=fv)
    formals(family$d2ldmdd, envir=new.env()) <- alist(sigma=sigma)  # ?? I do not think is needed
    
    ## initial values for sigma
    if(!is.null(sigma.start)) {
      if(length(sigma.start)>1){
        sigma <- sigma.start 
      } else {
        sigma <- rep(sigma.start,N)
      }
    } else {
      eval(family$sigma.initial)  # use function in family object to calculate initial values
    } 
  }
  
  ##*****************************************************************************************
  ##  Set up for the 3rd parameter submodel:
  ##            nu.X   design matrix
  ##       nu.offset   offset in linear predictor
  ##        nu.start   starting values for nu (optional)
  ##*****************************************************************************************
  if ("nu"%in%names(family$parameters)){
    nu.fit <- list() # MS Thursday, January 23, 2003 at 14:48    
    form.nu <- other.formula(form = nu.formula)  # add response to formula
    
    ## Specials for smoothing
    # add specials attribute for smoothing to formula object
    if(missing(data)){
      nu.terms <- terms(form.nu, specials = .gamlss.sm.list) 
    } else {
      nu.terms <- terms(form.nu, specials = .gamlss.sm.list, data = data)
    }
    mcall$formula <- nu.terms
    attr(mcall$formula, ".Environment") <- orig.Envir # DS fix for Willem Thursday, March 18, 2010
    
    nu.frame <- eval(mcall,sys.parent())  # calls the model.frame function inside mcall to create the 
    # modelframe with the variables needed to use formula
    # also uses pb() function to create model frame for smoothing   
    nu.terms <- attr(nu.frame, "terms")  # get terms attribute from modelframe 
    nu.a <- attributes(nu.terms)  # get attributes from terms object 
    nu.smoothers <- get.smoothers(nu.terms)  
    nu.X <- model.matrix(form.nu, nu.frame, contrasts)  # creates the design matrix for nu 
    nu.offset <- model.extract(nu.frame, offset)  # the nu-location offset
    if(is.null(nu.offset)){
      nu.offset <- rep(0,N)
    }
    nu.object <- get.object("nu")  # large list with different attributes (some inherited from family object)
    # the following is to get the right GLIM arguments 
    # replace function parameters ($y, $mu, $sigma, $nu, $tau) with (nu=fv)     
    formals(nu.object$dldp, envir=new.env()) <- alist(nu=fv)
    formals(nu.object$d2ldp2, envir=new.env()) <- alist(nu=fv)    
    formals(nu.object$G.di, envir=new.env()) <- alist(nu=fv)
    formals(nu.object$valid, envir=new.env()) <- alist(nu=fv)
    formals(family$d2ldmdv, envir=new.env()) <- alist(nu=nu)
    formals(family$d2ldddv, envir=new.env()) <- alist(nu=nu) 
    
    ## initial values for nu 
    if(!is.null(nu.start)){
      if(length(nu.start)>1){
        nu <- nu.start  
      } else {
        nu <- rep(nu.start,N)
      }
    } else {
      eval(family$nu.initial)  # use function in family object to calculate initial values
    }
  }    
  
  ##*****************************************************************************************
  ##  Set up for the 4rd parameter submodel:
  ##            tau.X   design matrix
  ##       tau.offset   offset in linear predictor
  ##        tau.start   starting values for tau (optional)
  ##*****************************************************************************************
  if ("tau"%in%names(family$parameters)){   
    tau.fit <- list() # MS Thursday, January 23, 2003 at 14:48
    form.tau <- other.formula(form = tau.formula)  # add response to formula
    
    ## Specials for smoothing
    # add specials attribute for smoothing to formula object
    if(missing(data)){
      tau.terms <- terms(form.tau, specials = .gamlss.sm.list)
    } else {
      tau.terms <- terms(form.tau, specials = .gamlss.sm.list, data = data)
    }
    mcall$formula <- tau.terms
    attr(mcall$formula, ".Environment") <- orig.Envir # DS fix for Willem Thursday, March 18, 2010
    
    tau.frame <- eval(mcall,sys.parent())  # calls the model.frame function inside mcall to create the 
    # modelframe with the variables needed to use formula
    # also uses pb() function to create model frame for smoothing
    tau.terms <- attr(tau.frame, "terms")  # get terms attribute from modelframe 
    tau.a <- attributes(tau.terms)  # get attributes from terms object
    tau.smoothers <- get.smoothers(tau.terms)  
    tau.X <- model.matrix(form.tau, tau.frame, contrasts)  # creates the design matrix for tau 
    tau.offset <- model.extract(tau.frame, offset)  # the tau-location offset
    if(is.null(tau.offset)){
      tau.offset <- rep(0,N)
    }
    tau.object <- get.object("tau")  # large list with different attributes (some inherited from family object)
    # the following is to get the right GLIM arguments 
    # replace function parameters ($y, $mu, $sigma, $nu, $tau) with (tau=fv)     
    formals(tau.object$dldp, envir=new.env()) <- alist(tau=fv) 
    formals(tau.object$d2ldp2, envir=new.env()) <- alist(tau=fv) 
    formals(tau.object$G.di, envir=new.env()) <- alist(tau=fv)
    formals(tau.object$valid, envir=new.env()) <- alist(tau=fv)
    formals(family$d2ldmdt, envir=new.env()) <- alist(tau=tau)
    formals(family$d2ldddt, envir=new.env()) <- alist(tau=tau)  
    formals(family$d2ldvdt, envir=new.env()) <- alist(tau=tau) 
    
    ## initial values for tau
    if(!is.null(tau.start)) {
      if(length(tau.start)>1){
        tau <- tau.start
      } else {
        tau <- rep(tau.start,N)
      }
    } else {
      eval(family$tau.initial)  # use function in family object to calculate initial values
    }  
  } 
  
  #******************************************************************************************
  #V) Fit Model ----
  #******************************************************************************************
  fiter <- 0  # number of outer iteration
  conv <- eval(substitute(method))  # execute method (returns TRUE if algorithm converged)
  method <- substitute(method)  # return parse tree substituting variables bound in current environment
  
  #******************************************************************************************
  #VI) Create output ----
  #******************************************************************************************
  ## first the general output 
  # calculate the Global deviance again
  G.dev.incr  <- eval(G.dev.expr)  # deviance increment (function provided by gamlss.family object)
  G.dev <- sum(w*G.dev.incr)  # the weighted global deviance
  
  out <- list(family=family$family, parameters=names(family$parameters), 
              call=gamlsscall, y=y, control=control, weights=w, 
              G.deviance=G.dev, N=N, rqres=family$rqres, iter=fiter, 
              type=family$type, method=method, contrasts=contrasts) 
  # family : name of the family for the response variable
  # parameters : names of the family parameters
  # call : gamlss call object with parameters & data set
  # y : response variable
  # control : control parameters for outer iteration
  # weights : weights for weighted likelihood analysis (not the same as in GLM's)
  # G.deviance : weighted global deviance
  # N : length of the response variable (= number of observations if no weights are used)
  # reqres : function to calculate the normalized (randomized) quantile residuals of the object
  # iter : number of outer iterations
  # type : type of the distribution of the response variable (continuous, discrete or mixed)
  # method : Which algorithm is used for the fit, RS(), CG() or mixed()
  # contrasts : contrasts that were used in the fit
  
  out$converged <- conv  # Whether the model has converged
  out$residuals <- eval(family$rqres)  # normalized (randomized) quantile residuals of the model
  if(all(trunc(w)==w)){  # frequency count weights are used
    noObs <- sum(w)
  } else {
    noObs <- N 
  }
  out$noObs <- noObs  # actual number of observations
  
  ## binomial denominator 
  if(any(family$family%in%.gamlss.bi.list)){
    out$bd <- bd  # denominator for binomial family
  }
  ##*****************************************************************************************
  saveParam <- control$save  # Whether all information should be saved on exit. 
  # If save=FALSE only the global deviance & AIV are saved
  
  ## Output for mean model
  if ("mu"%in%names(family$parameters)){
    out <- c(out, mu=parameterOut(what="mu", save=saveParam))
  } else {
    out$mu.df <- 0
  }
  # define now the overall degrees of freedom for the fit and for the residuals
  out$df.fit <- out$mu.df  # total df (for all parameters)
  out$df.residual <- noObs-out$mu.df  # residual df left after model is fitted
  out$pen <- out$mu.pen  # sum of quadratic penalties (for all parameters)
  
  ## Output for dispersion model
  if ("sigma"%in%names(family$parameters)){
    out <- c(out, sigma = parameterOut(what="sigma", save=saveParam))  
    out$df.fit <- out$mu.df + out$sigma.df  # total df (for all parameters)
    out$df.residual <- noObs-out$mu.df-out$sigma.df  # residual df left after model is fitted
    out$pen <- out$mu.pen + out$sigma.pen  # sum of quadratic penalties (for all parameters)
  }
  
  ## Output for nu model
  if ("nu"%in%names(family$parameters)){
    out <- c(out, nu = parameterOut(what="nu", save=saveParam)) 
    out$df.fit <- out$mu.df+out$sigma.df+out$nu.df  # total df (for all parameters)
    out$df.residual <- noObs-out$mu.df-out$sigma.df-out$nu.df  # residual df left after model is fitted
    out$pen <- out$mu.pen + out$sigma.pen + out$nu.pen  # sum of quadratic penalties (for all parameters)
  }
  
  ## Output for tau model
  if ("tau"%in%names(family$parameters)){
    out <- c(out, tau = parameterOut(what="tau", save=saveParam))       
    out$df.fit <- out$mu.df+out$sigma.df+out$nu.df+out$tau.df  # total df (for all parameters)
    out$df.residual <- noObs-out$mu.df-out$sigma.df- out$nu.df -out$tau.df  # residual df left after model is fitted
    out$pen <- out$mu.pen + out$sigma.pen + out$nu.pen + out$tau.pen  # sum of quadratic penalties (for all parameters)
  }
  
  out$P.deviance <- out$G.deviance+out$pen  # penalized deviance
  out$aic <- G.dev+2*out$df.fit  # aic (k=2)
  out$sbc <- G.dev+log(noObs)*out$df.fit  # sbc (k=log(n))
  class(out) <- c("gamlss","gam","glm","lm")  # change the class of the gamlss object
  out
}
## the END of gamlss

