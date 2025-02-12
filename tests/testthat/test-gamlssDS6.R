# Test-smk-gamlssDS6 ################################################

#
# Set up
#

set.standard.disclosure.settings()
set.random.seed.setting(379)

#
# Tests
#

test_that("output_gamlssDS6_innerit1_autostep0", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss1.rda"))
  outputDS1 <- gamlssDS1(formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         sigma.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         nu.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         tau.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis",
                         family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", mu.fix=FALSE, sigma.fix=FALSE, 
                         nu.fix=FALSE, tau.fix=FALSE, global.mean=mean(gamlss1$e3_bw), global.sd=NULL, 
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","), autostep=TRUE)
  smoother.xl <- outputDS1$smoother.xmin - 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  smoother.xr <- outputDS1$smoother.xmax + 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  outputDS2 <- gamlssDS2(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(c(0,0)), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  beta.vect.next <- as.vector(solve(outputDS2$matrix) %*% outputDS2$vector)
  outputDS3 <- gamlssDS3(parameter="mu", smoother=1, family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$mu.coefSmo[[1]]$coef)[1])), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  lambda <- outputDS1$mod.gamlss.ds$mu.coefSmo[[1]]$lambda
  D.mat <- diff(diag(20+3), diff=2)
  gamma.vect.update <- as.vector(solve(outputDS3$matrix + lambda*t(D.mat)%*%D.mat) %*% outputDS3$vector)
  outputDS4 <- gamlssDS4(parameter="mu", smoother=1, family="BCPEleft_parenthesisright_parenthesis", data="gamlss1",
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  outputDS5 <- gamlssDS5(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  outputDS6 <- gamlssDS6(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","),
                         autostep=TRUE, inner.iteration.count=1, autostep.count=0)
  Z.mat <- bbase(x=gamlss1$e3_gac_None, xl=smoother.xl, xr=smoother.xr, ndx=20, deg=3)
  s <- as.matrix(Z.mat %*% gamma.vect.update)
  eta.old <- BCPE()$mu.linkfun(temp_mu)
  eta <- as.vector(temp_mod.gamlss.ds$mu.x %*% beta.vect.next + base::rowSums(as.matrix(s)))
  fv <- BCPE()$mu.linkinv(eta)
  dr <- 1/BCPE()$mu.dr(eta)
  dv <- sum(BCPE()$G.dev.incr(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau))
  dldp <- BCPE()$dldm(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- BCPE()$d2ldm2(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2/(dr*dr))
  wt <- ifelse(wt>1e+10,1e+10,wt)
  wt <- ifelse(wt<1e-10,1e-10,wt)
  wv <- eta+dldp/(dr*wt)
  pen <-  sum(eta*wt*(wv-eta))
  errorMessage <- NULL
  valid <- BCPE()$mu.valid(temp_mu)
  if (is.na(!valid)){
    errorMessage <- "Fitted values in the inner iteration out of range."
  } 
  
  # output
  expect_equal(outputDS6$dv, dv, tolerance=1e-07)
  expect_equal(outputDS6$pen, pen, tolerance=1e-07)
  expect_equal(outputDS6$errorMessage3, NULL)
  
})


test_that("output_gamlssDS6_innerit2_autostep1", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss1.rda"))
  outputDS1 <- gamlssDS1(formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         sigma.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         nu.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         tau.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis",
                         family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", mu.fix=FALSE, sigma.fix=FALSE, 
                         nu.fix=FALSE, tau.fix=FALSE, global.mean=mean(gamlss1$e3_bw), global.sd=NULL, 
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","), autostep=TRUE)
  smoother.xl <- outputDS1$smoother.xmin - 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  smoother.xr <- outputDS1$smoother.xmax + 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  outputDS2 <- gamlssDS2(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(c(0,0)), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  beta.vect.next <- as.vector(solve(outputDS2$matrix) %*% outputDS2$vector)
  outputDS3 <- gamlssDS3(parameter="mu", smoother=1, family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$mu.coefSmo[[1]]$coef)[1])), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  lambda <- outputDS1$mod.gamlss.ds$mu.coefSmo[[1]]$lambda
  D.mat <- diff(diag(20+3), diff=2)
  gamma.vect.update <- as.vector(solve(outputDS3$matrix + lambda*t(D.mat)%*%D.mat) %*% outputDS3$vector)
  outputDS4 <- gamlssDS4(parameter="mu", smoother=1, family="BCPEleft_parenthesisright_parenthesis", data="gamlss1",
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  outputDS5 <- gamlssDS5(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  eta.old <- BCPE()$mu.linkfun(temp_mu)
  s.old <- temp_mu.s
  outputDS6 <- gamlssDS6(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","),
                         autostep=TRUE, inner.iteration.count=2, autostep.count=1)
  Z.mat <- bbase(x=gamlss1$e3_gac_None, xl=smoother.xl, xr=smoother.xr, ndx=20, deg=3)
  s <- as.matrix(Z.mat %*% gamma.vect.update)
  eta <- as.vector(temp_mod.gamlss.ds$mu.x %*% beta.vect.next + base::rowSums(as.matrix(s)))
  eta <- (eta.old + eta)/2
  s <- (s.old + s)/2
  fv <- BCPE()$mu.linkinv(eta)
  dr <- 1/BCPE()$mu.dr(eta)
  dv <- sum(BCPE()$G.dev.incr(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau))
  dldp <- BCPE()$dldm(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- BCPE()$d2ldm2(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2/(dr*dr))
  wt <- ifelse(wt>1e+10,1e+10,wt)
  wt <- ifelse(wt<1e-10,1e-10,wt)
  wv <- eta+dldp/(dr*wt)
  pen <-  sum(eta*wt*(wv-eta))
  errorMessage <- NULL
  valid <- BCPE()$mu.valid(temp_mu)
  if (is.na(!valid)){
    errorMessage <- "Fitted values in the inner iteration out of range."
  } 
  
  # output
  expect_equal(outputDS6$dv, dv, tolerance=1e-07)
  expect_equal(outputDS6$pen, pen, tolerance=1e-07)
  expect_equal(outputDS6$errorMessage3, NULL)
  
  # assignments
  expect_equal(temp_mu.s, s, tolerance=1e-07)
  expect_equal(temp_mu, fv, tolerance=1e-07)
  
})

test_that("output_gamlssDS6_innerit1_autostep2", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss1.rda"))
  outputDS1 <- gamlssDS1(formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         sigma.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         nu.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis", 
                         tau.formula="gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis",
                         family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", mu.fix=FALSE, sigma.fix=FALSE, 
                         nu.fix=FALSE, tau.fix=FALSE, global.mean=mean(gamlss1$e3_bw), global.sd=NULL, 
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","), autostep=TRUE)
  smoother.xl <- outputDS1$smoother.xmin - 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  smoother.xr <- outputDS1$smoother.xmax + 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  outputDS2 <- gamlssDS2(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(c(0,0)), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  beta.vect.next <- as.vector(solve(outputDS2$matrix) %*% outputDS2$vector)
  outputDS3 <- gamlssDS3(parameter="mu", smoother=1, family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$mu.coefSmo[[1]]$coef)[1])), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  lambda <- outputDS1$mod.gamlss.ds$mu.coefSmo[[1]]$lambda
  D.mat <- diff(diag(20+3), diff=2)
  gamma.vect.update <- as.vector(solve(outputDS3$matrix + lambda*t(D.mat)%*%D.mat) %*% outputDS3$vector)
  outputDS4 <- gamlssDS4(parameter="mu", smoother=1, family="BCPEleft_parenthesisright_parenthesis", data="gamlss1",
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  outputDS5 <- gamlssDS5(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","))
  eta.old <- BCPE()$mu.linkfun(temp_mu)
  s.old <- temp_mu.s
  outputDS6 <- gamlssDS6(parameter="mu", family="BCPEleft_parenthesisright_parenthesis", data="gamlss1", 
                         mu.beta.vect=paste0(as.character(beta.vect.next), collapse=","), sigma.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         nu.beta.vect=paste0(as.character(c(0,0)), collapse=","), tau.beta.vect=paste0(as.character(c(0,0)), collapse=","), 
                         mu.gamma.vect=paste0(as.character(gamma.vect.update), collapse=","),
                         sigma.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse=","),
                         nu.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse=","),
                         tau.gamma.vect=paste0(as.character(rep(0, times=dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse=","),
                         smoother.names=outputDS1$smoother.names, smoother.xl=paste0(as.character(smoother.xl), collapse=","),
                         smoother.xr=paste0(as.character(smoother.xr), collapse=","),
                         control=paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse=","),
                         i.control=paste0(as.character(c(0.001, 50, 30, 0.001)), collapse=","),
                         autostep=TRUE, inner.iteration.count=1, autostep.count=2)
  Z.mat <- bbase(x=gamlss1$e3_gac_None, xl=smoother.xl, xr=smoother.xr, ndx=20, deg=3)
  s <- as.matrix(Z.mat %*% gamma.vect.update)
  eta <- as.vector(temp_mod.gamlss.ds$mu.x %*% beta.vect.next + base::rowSums(as.matrix(s)))
  eta <- (eta.old*(2**2 - 1) - eta)/(2**2-2)
  s <- (s.old*(2**2 - 1) - s)/(2**2-2)
  fv <- BCPE()$mu.linkinv(eta)
  dr <- 1/BCPE()$mu.dr(eta)
  dv <- sum(BCPE()$G.dev.incr(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau))
  dldp <- BCPE()$dldm(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- BCPE()$d2ldm2(gamlss1$e3_bw, fv, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2/(dr*dr))
  wt <- ifelse(wt>1e+10,1e+10,wt)
  wt <- ifelse(wt<1e-10,1e-10,wt)
  wv <- eta+dldp/(dr*wt)
  pen <-  sum(eta*wt*(wv-eta))
  errorMessage <- NULL
  valid <- BCPE()$mu.valid(temp_mu)
  if (is.na(!valid)){
    errorMessage <- "Fitted values in the inner iteration out of range."
  } 
  
  # output
  expect_equal(outputDS6$dv, dv, tolerance=1e-07)
  expect_equal(outputDS6$pen, pen, tolerance=1e-07)
  expect_equal(outputDS6$errorMessage3, NULL)
  
  # assignments
  expect_equal(temp_mu.s, s, tolerance=1e-07)
  expect_equal(temp_mu, fv, tolerance=1e-07)
  
})



#
# Done
#


