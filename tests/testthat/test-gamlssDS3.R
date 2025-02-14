# Test-smk-gamlssDS3 ################################################

#
# Set up
#

set.standard.disclosure.settings()
set.random.seed.setting(379)

#
# Tests
#

test_that("output_gamlssDS3", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss1.rda"))
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
  smoother.xl <- outputDS1$smoother.xmin - 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  smoother.xr <- outputDS1$smoother.xmax + 0.01 * (outputDS1$smoother.xmax - outputDS1$smoother.xmin)
  outputDS2 <- gamlssDS2(
    parameter = "mu", family = "BCPEleft_parenthesisright_parenthesis", data = "gamlss1",
    mu.beta.vect = paste0(as.character(c(0, 0)), collapse = ","), sigma.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    nu.beta.vect = paste0(as.character(c(0, 0)), collapse = ","), tau.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    control = paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse = ","),
    i.control = paste0(as.character(c(0.001, 50, 30, 0.001)), collapse = ",")
  )
  beta.vect.next <- as.vector(solve(outputDS2$matrix) %*% outputDS2$vector)
  outputDS3 <- gamlssDS3(
    parameter = "mu", smoother = 1, family = "BCPEleft_parenthesisright_parenthesis", data = "gamlss1",
    mu.beta.vect = paste0(as.character(beta.vect.next), collapse = ","), sigma.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    nu.beta.vect = paste0(as.character(c(0, 0)), collapse = ","), tau.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    mu.gamma.vect = paste0(as.character(rep(0, times = dim(outputDS1$mod.gamlss.ds$mu.coefSmo[[1]]$coef)[1])), collapse = ","),
    sigma.gamma.vect = paste0(as.character(rep(0, times = dim(outputDS1$mod.gamlss.ds$sigma.coefSmo[[1]]$coef)[1])), collapse = ","),
    nu.gamma.vect = paste0(as.character(rep(0, times = dim(outputDS1$mod.gamlss.ds$nu.coefSmo[[1]]$coef)[1])), collapse = ","),
    tau.gamma.vect = paste0(as.character(rep(0, times = dim(outputDS1$mod.gamlss.ds$tau.coefSmo[[1]]$coef)[1])), collapse = ","),
    smoother.names = outputDS1$smoother.names, smoother.xl = paste0(as.character(smoother.xl), collapse = ","),
    smoother.xr = paste0(as.character(smoother.xr), collapse = ","),
    control = paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse = ","),
    i.control = paste0(as.character(c(0.001, 50, 30, 0.001)), collapse = ",")
  )
  Z.mat <- bbase(x = gamlss1$e3_gac_None, xl = smoother.xl, xr = smoother.xr, ndx = 20, deg = 3)
  eta <- BCPE()$mu.linkfun(temp_mu)
  dr <- 1 / BCPE()$mu.dr(eta)
  dldp <- BCPE()$dldm(gamlss1$e3_bw, temp_mu, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- BCPE()$d2ldm2(gamlss1$e3_bw, temp_mu, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2 / (dr * dr))
  wt <- ifelse(wt > 1e+10, 1e+10, wt)
  wt <- ifelse(wt < 1e-10, 1e-10, wt)
  wv <- eta + dldp / (dr * wt)
  partial.residuals <- wv - temp_mod.gamlss.ds$mu.x %*% beta.vect.next - base::rowSums(as.matrix(temp_mu.s[, -1]))

  vector <- t(Z.mat) %*% (wt * partial.residuals)
  matrix <- t(Z.mat) %*% diag(wt) %*% Z.mat

  # output
  expect_equal(outputDS3$matrix, matrix, ignore_attr = TRUE)
  expect_equal(outputDS3$vector, vector, ignore_attr = TRUE)
  expect_equal(outputDS3$sumofsquares, 0)
  expect_equal(outputDS3$sumofweights, 1)

  # assignments
  expect_equal(temp_mu.s, as.matrix(rep(0, times = nrow(gamlss1))))
})

#
# Done
#
