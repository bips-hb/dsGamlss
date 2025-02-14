# Test-discctrl-gamlssDS1 ################################################

#
# Set up
#

set.standard.disclosure.settings()
set.random.seed.setting(379)

#
# Tests
#

test_that("discctrl_gamlssDS1", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss1.rda"))
  gamlss1_red <- gamlss1[1:28, ]
  gamlss1_red$binary_disc <- c(rep(0, times = 26), 1, 1)
  outputDS1 <- gamlssDS1(
    formula = "gamlss1_red$e3_bw tilde_symbol pbleft_parenthesisgamlss1_red$e3_gac_Noneright_parenthesis",
    sigma.formula = "gamlss1_red$e3_bw tilde_symbol binary_disc",
    nu.formula = "gamlss1_red$e3_bw tilde_symbol pbleft_parenthesisgamlss1_red$e3_gac_Noneright_parenthesis",
    tau.formula = "gamlss1_red$e3_bw tilde_symbol pbleft_parenthesisgamlss1_red$e3_gac_Noneright_parenthesis",
    family = "BCPEleft_parenthesisright_parenthesis", data = "gamlss1_red", mu.fix = FALSE, sigma.fix = FALSE,
    nu.fix = FALSE, tau.fix = FALSE, global.mean = mean(gamlss1_red$e3_bw), global.sd = NULL,
    control = paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse = ","),
    i.control = paste0(as.character(c(0.001, 50, 30, 0.001)), collapse = ","), autostep = TRUE
  )
  expect_true(is.na(outputDS1$mod.gamlss.ds))
  expect_true(is.na(outputDS1$G.dev))
  expect_true(is.na(outputDS1$dim.mu.x))
  expect_true(is.na(outputDS1$dim.sigma.x))
  expect_true(is.na(outputDS1$dim.nu.x))
  expect_true(is.na(outputDS1$dim.tau.x))
  expect_true(is.na(outputDS1$smoother.names))
  expect_true(is.na(outputDS1$smoother.xmin))
  expect_true(is.na(outputDS1$smoother.xmax))
  expect_equal(outputDS1$y.invalid, 0)
  expect_equal(outputDS1$mu.par.invalid, c(0, 0))
  expect_equal(outputDS1$sigma.par.invalid, c(0, 1))
  expect_equal(outputDS1$nu.par.invalid, c(0, 0))
  expect_equal(outputDS1$tau.par.invalid, c(0, 0))
  expect_equal(outputDS1$gamlss.saturation.invalid, 1)
  expect_equal(outputDS1$errorMessage, "Study data or applied model invalid for this source")
  expect_equal(length(ls(pattern = "^temp_")), 0)
})

#
# Done
#


# Test-smk-gamlssDS1 ################################################

#
# Set up
#

set.standard.disclosure.settings()
set.random.seed.setting(379)

#
# Tests
#

test_that("outputDS1_gamlssDS1", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss1.rda"))
  mod.gamlss <- base::suppressWarnings(gamlss::gamlss(
    formula = e3_bw ~ pb(e3_gac_None), sigma.formula = e3_bw ~ pb(e3_gac_None),
    nu.formula = e3_bw ~ pb(e3_gac_None), tau.formula = e3_bw ~ pb(e3_gac_None),
    family = BCPE(), data = gamlss1, mu.fix = FALSE, sigma.fix = FALSE, nu.fix = FALSE,
    tau.fix = FALSE, control = gamlss::gamlss.control(n.cyc = 1, trace = FALSE, autostep = TRUE),
    i.control = gamlss::glim.control(cyc = 1, bf.cyc = 1)
  ))
  outputDS1 <- gamlssDS1(
    formula = "gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis",
    sigma.formula = "gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis",
    nu.formula = "gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis",
    tau.formula = "gamlss1$e3_bw tilde_symbol pbleft_parenthesisgamlss1$e3_gac_Noneright_parenthesis",
    family = "BCPEleft_parenthesisright_parenthesis", data = "gamlss1", mu.fix = FALSE, sigma.fix = FALSE,
    nu.fix = FALSE, tau.fix = FALSE, global.mean = mean(gamlss1$e3_bw), global.sd = NULL,
    control = paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse = ","),
    i.control = paste0(as.character(c(0.001, 50, 30, 0.001)), collapse = ","), autostep = TRUE
  )
  set.seed(379)

  # mod.gamlss.ds
  mod.gamlss.ds <- outputDS1$mod.gamlss.ds
  expect_length(mod.gamlss.ds, length(mod.gamlss))
  expect_equal(mod.gamlss.ds$family, mod.gamlss$family)
  expect_equal(mod.gamlss.ds$parameters, mod.gamlss$parameters)
  expect_equal(mod.gamlss.ds$y, "The response variable is not disclosed!")
  expect_equal(mod.gamlss.ds$control[c(-2, -9)], mod.gamlss$control[c(-2, -9)])
  expect_equal(mod.gamlss.ds$control$n.cyc, 20)
  expect_equal(mod.gamlss.ds$control$trace, TRUE)
  expect_equal(mod.gamlss.ds$weights, mod.gamlss$weights)
  expect_equal(mod.gamlss.ds$G.deviance, mod.gamlss$G.deviance, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$N, mod.gamlss$N)
  expect_equal(mod.gamlss.ds$rqres, mod.gamlss$rqres)
  expect_equal(mod.gamlss.ds$iter, mod.gamlss$iter)
  expect_equal(mod.gamlss.ds$type, mod.gamlss$type)
  expect_equal(mod.gamlss.ds$method, mod.gamlss$method)
  expect_equal(mod.gamlss.ds$contrasts, mod.gamlss$contrasts)
  expect_equal(mod.gamlss.ds$converged, mod.gamlss$converged)
  expect_equal(mod.gamlss.ds$residuals, "The residuals of the model are not disclosed!")
  expect_equal(mod.gamlss.ds$noObs, mod.gamlss$noObs)
  expect_equal(mod.gamlss.ds$mu.fv, "The fitted values of the mu model are not disclosed!")
  expect_equal(mod.gamlss.ds$mu.lp, "The linear predictors of the mu model are not disclosed!")
  expect_equal(mod.gamlss.ds$mu.wv, "The working variable of the mu model are not disclosed!")
  expect_equal(mod.gamlss.ds$mu.wt, "The working weights of the mu model are not disclosed!")
  expect_equal(mod.gamlss.ds$mu.link, mod.gamlss$mu.link)
  expect_equal(mod.gamlss.ds$mu.x, "The design matrix of the mu model is not disclosed!")
  expect_equal(mod.gamlss.ds$mu.qr, "The QR decomposition of the mu model is not disclosed!")
  expect_equal(mod.gamlss.ds$mu.coefficients, mod.gamlss$mu.coefficients, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$mu.offset, mod.gamlss$mu.offset)
  expect_equal(mod.gamlss.ds$mu.xlevels, mod.gamlss$mu.xlevels)
  expect_equal(mod.gamlss.ds$mu.formula, ~1, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$mu.df, mod.gamlss$mu.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.nl.df, mod.gamlss$mu.nl.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.s, "The smoothing fitted values of the mu model are not disclosed!")
  expect_equal(mod.gamlss.ds$mu.var, "The variances for the smoothing fitted values of the mu model are not disclosed!")
  expect_length(mod.gamlss.ds$mu.coefSmo[[1]], length(mod.gamlss$mu.coefSmo[[1]]))
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$coef, mod.gamlss$mu.coefSmo[[1]]$coef, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$fv, "The smoothing fitted values of the mu model are not disclosed!")
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$lambda, mod.gamlss$mu.coefSmo[[1]]$lambda, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$edf, mod.gamlss$mu.coefSmo[[1]]$edf, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$sigb, mod.gamlss$mu.coefSmo[[1]]$sigb, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$sige, mod.gamlss$mu.coefSmo[[1]]$sige, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$method, mod.gamlss$mu.coefSmo[[1]]$method)
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$knots, mod.gamlss$mu.coefSmo[[1]]$knots, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$mu.coefSmo[[1]]$fun, "The function for the knots of the mu model is not disclosed!")
  expect_equal(mod.gamlss.ds$mu.pen, mod.gamlss$mu.pen, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$df.fit, mod.gamlss$df.fit, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$pen, mod.gamlss$pen, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$df.residual, mod.gamlss$df.residual, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.fv, "The fitted values of the sigma model are not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.lp, "The linear predictors of the sigma model are not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.wv, "The working variable of the sigma model are not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.wt, "The working weights of the sigma model are not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.link, mod.gamlss$sigma.link)
  expect_equal(mod.gamlss.ds$sigma.x, "The design matrix of the sigma model is not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.qr, "The QR decomposition of the sigma model is not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.coefficients, mod.gamlss$sigma.coefficients, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$sigma.offset, mod.gamlss$sigma.offset)
  expect_equal(mod.gamlss.ds$sigma.xlevels, mod.gamlss$sigma.xlevels)
  expect_equal(mod.gamlss.ds$sigma.formula, ~1, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$sigma.df, mod.gamlss$sigma.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.nl.df, mod.gamlss$sigma.nl.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.s, "The smoothing fitted values of the sigma model are not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.var, "The variances for the smoothing fitted values of the sigma model are not disclosed!")
  expect_length(mod.gamlss.ds$sigma.coefSmo[[1]], length(mod.gamlss$sigma.coefSmo[[1]]))
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$coef, mod.gamlss$sigma.coefSmo[[1]]$coef, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$fv, "The smoothing fitted values of the sigma model are not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$lambda, mod.gamlss$sigma.coefSmo[[1]]$lambda, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$edf, mod.gamlss$sigma.coefSmo[[1]]$edf, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$sigb, mod.gamlss$sigma.coefSmo[[1]]$sigb, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$sige, mod.gamlss$sigma.coefSmo[[1]]$sige, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$method, mod.gamlss$sigma.coefSmo[[1]]$method)
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$knots, mod.gamlss$sigma.coefSmo[[1]]$knots, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sigma.coefSmo[[1]]$fun, "The function for the knots of the sigma model is not disclosed!")
  expect_equal(mod.gamlss.ds$sigma.pen, mod.gamlss$sigma.pen, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.fv, "The fitted values of the nu model are not disclosed!")
  expect_equal(mod.gamlss.ds$nu.lp, "The linear predictors of the nu model are not disclosed!")
  expect_equal(mod.gamlss.ds$nu.wv, "The working variable of the nu model are not disclosed!")
  expect_equal(mod.gamlss.ds$nu.wt, "The working weights of the nu model are not disclosed!")
  expect_equal(mod.gamlss.ds$nu.link, mod.gamlss$nu.link)
  expect_equal(mod.gamlss.ds$nu.x, "The design matrix of the nu model is not disclosed!")
  expect_equal(mod.gamlss.ds$nu.qr, "The QR decomposition of the nu model is not disclosed!")
  expect_equal(mod.gamlss.ds$nu.coefficients, mod.gamlss$nu.coefficients, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$nu.offset, mod.gamlss$nu.offset)
  expect_equal(mod.gamlss.ds$nu.xlevels, mod.gamlss$nu.xlevels)
  expect_equal(mod.gamlss.ds$nu.formula, ~1, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$nu.df, mod.gamlss$nu.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.nl.df, mod.gamlss$nu.nl.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.s, "The smoothing fitted values of the nu model are not disclosed!")
  expect_equal(mod.gamlss.ds$nu.var, "The variances for the smoothing fitted values of the nu model are not disclosed!")
  expect_length(mod.gamlss.ds$nu.coefSmo[[1]], length(mod.gamlss$nu.coefSmo[[1]]))
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$coef, mod.gamlss$nu.coefSmo[[1]]$coef, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$fv, "The smoothing fitted values of the nu model are not disclosed!")
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$lambda, mod.gamlss$nu.coefSmo[[1]]$lambda, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$edf, mod.gamlss$nu.coefSmo[[1]]$edf, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$sigb, mod.gamlss$nu.coefSmo[[1]]$sigb, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$sige, mod.gamlss$nu.coefSmo[[1]]$sige, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$method, mod.gamlss$nu.coefSmo[[1]]$method)
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$knots, mod.gamlss$nu.coefSmo[[1]]$knots, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$nu.coefSmo[[1]]$fun, "The function for the knots of the nu model is not disclosed!")
  expect_equal(mod.gamlss.ds$nu.pen, mod.gamlss$nu.pen, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.fv, "The fitted values of the tau model are not disclosed!")
  expect_equal(mod.gamlss.ds$tau.lp, "The linear predictors of the tau model are not disclosed!")
  expect_equal(mod.gamlss.ds$tau.wv, "The working variable of the tau model are not disclosed!")
  expect_equal(mod.gamlss.ds$tau.wt, "The working weights of the tau model are not disclosed!")
  expect_equal(mod.gamlss.ds$tau.link, mod.gamlss$tau.link)
  expect_equal(mod.gamlss.ds$tau.x, "The design matrix of the tau model is not disclosed!")
  expect_equal(mod.gamlss.ds$tau.qr, "The QR decomposition of the tau model is not disclosed!")
  expect_equal(mod.gamlss.ds$tau.coefficients, mod.gamlss$tau.coefficients, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$tau.offset, mod.gamlss$tau.offset)
  expect_equal(mod.gamlss.ds$tau.xlevels, mod.gamlss$tau.xlevels)
  expect_equal(mod.gamlss.ds$tau.formula, ~1, ignore_attr = TRUE)
  expect_equal(mod.gamlss.ds$tau.df, mod.gamlss$tau.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.nl.df, mod.gamlss$tau.nl.df, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.s, "The smoothing fitted values of the tau model are not disclosed!")
  expect_equal(mod.gamlss.ds$tau.var, "The variances for the smoothing fitted values of the tau model are not disclosed!")
  expect_length(mod.gamlss.ds$tau.coefSmo[[1]], length(mod.gamlss$tau.coefSmo[[1]]))
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$coef, mod.gamlss$tau.coefSmo[[1]]$coef, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$fv, "The smoothing fitted values of the tau model are not disclosed!")
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$lambda, mod.gamlss$tau.coefSmo[[1]]$lambda, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$edf, mod.gamlss$tau.coefSmo[[1]]$edf, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$sigb, mod.gamlss$tau.coefSmo[[1]]$sigb, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$sige, mod.gamlss$tau.coefSmo[[1]]$sige, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$method, mod.gamlss$tau.coefSmo[[1]]$method)
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$knots, mod.gamlss$tau.coefSmo[[1]]$knots, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$tau.coefSmo[[1]]$fun, "The function for the knots of the tau model is not disclosed!")
  expect_equal(mod.gamlss.ds$tau.pen, mod.gamlss$tau.pen, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$P.deviance, mod.gamlss$P.deviance, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$aic, mod.gamlss$aic, tolerance = 1e-07)
  expect_equal(mod.gamlss.ds$sbc, mod.gamlss$sbc, tolerance = 1e-07)

  # other outputDS1
  expect_equal(outputDS1$G.dev, 8576.428, tolerance = 1e-07)
  expect_equal(outputDS1$dim.mu.x, dim(mod.gamlss$mu.x))
  expect_equal(outputDS1$dim.sigma.x, dim(mod.gamlss$sigma.x))
  expect_equal(outputDS1$dim.nu.x, dim(mod.gamlss$nu.x))
  expect_equal(outputDS1$dim.tau.x, dim(mod.gamlss$tau.x))
  expect_equal(outputDS1$smoother.names, "gamlss1$e3_gac_None")
  expect_equal(outputDS1$smoother.xmin, min(gamlss1$e3_gac_None) - abs(stats::rnorm(n = 1, mean = 0, sd = sqrt(0.25 * stats::var(gamlss1$e3_gac_None)))))
  expect_equal(outputDS1$smoother.xmax, max(gamlss1$e3_gac_None) + abs(stats::rnorm(n = 1, mean = 0, sd = sqrt(0.25 * stats::var(gamlss1$e3_gac_None)))))
  expect_equal(outputDS1$y.invalid, 0)
  expect_equal(outputDS1$mu.par.invalid, c(0, 0))
  expect_equal(outputDS1$sigma.par.invalid, c(0, 0))
  expect_equal(outputDS1$nu.par.invalid, c(0, 0))
  expect_equal(outputDS1$tau.par.invalid, c(0, 0))
  expect_equal(outputDS1$gamlss.saturation.invalid, 0)
  expect_equal(outputDS1$errorMessage, "No errors")

  # assignments
  expect_equal(temp_mu, (gamlss1$e3_bw + mean(gamlss1$e3_bw)) / 2)
  expect_equal(temp_sigma, rep(0.1, nrow(gamlss1)))
  expect_equal(temp_nu, rep(1, nrow(gamlss1)))
  expect_equal(temp_tau, rep(2, nrow(gamlss1)))
  expect_equal(temp_mu.s, as.matrix(rep(0, times = nrow(gamlss1))))
  expect_equal(temp_sigma.s, as.matrix(rep(0, times = nrow(gamlss1))))
  expect_equal(temp_nu.s, as.matrix(rep(0, times = nrow(gamlss1))))
  expect_equal(temp_tau.s, as.matrix(rep(0, times = nrow(gamlss1))))
})

#
# Done
#
