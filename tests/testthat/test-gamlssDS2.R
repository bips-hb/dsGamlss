# Test-discctrl-gamlssDS2 ################################################

#
# Set up
#

set.standard.disclosure.settings()

#
# Tests
#

test_that("discctrl_gamlssDS2", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss1.rda"))
  gamlss1_red <- gamlss1[1:28, ]
  gamlss1_red$binary_disc <- c(rep(0, times = 26), 1, 1)
  temp_mod.gamlss.ds <- base::suppressWarnings(gamlss::gamlss(
    formula = e3_bw ~ pb(e3_gac_None), sigma.formula = e3_bw ~ binary_disc,
    nu.formula = e3_bw ~ pb(e3_gac_None), tau.formula = e3_bw ~ pb(e3_gac_None),
    family = BCPE(), data = gamlss1_red, mu.fix = FALSE, sigma.fix = FALSE, nu.fix = FALSE,
    tau.fix = FALSE, control = gamlss::gamlss.control(n.cyc = 1, trace = FALSE, autostep = TRUE),
    i.control = gamlss::glim.control(cyc = 1, bf.cyc = 1)
  ))
  temp_mu <- (gamlss1_red$e3_bw + mean(gamlss1_red$e3_bw)) / 2
  temp_sigma <- rep(0.1, nrow(gamlss1_red))
  temp_nu <- rep(1, nrow(gamlss1_red))
  temp_tau <- rep(2, nrow(gamlss1_red))
  temp_mu.s <- as.matrix(rep(0, times = nrow(gamlss1_red)))
  temp_sigma.s <- as.matrix(rep(0, times = nrow(gamlss1_red)))
  temp_nu.s <- as.matrix(rep(0, times = nrow(gamlss1_red)))
  temp_tau.s <- as.matrix(rep(0, times = nrow(gamlss1_red)))
  outputDS2 <- gamlssDS2(
    parameter = "mu", family = "BCPEleft_parenthesisright_parenthesis", data = "gamlss1_red",
    mu.beta.vect = paste0(as.character(c(0, 0)), collapse = ","), sigma.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    nu.beta.vect = paste0(as.character(c(0, 0)), collapse = ","), tau.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    control = paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse = ","),
    i.control = paste0(as.character(c(0.001, 50, 30, 0.001)), collapse = ",")
  )

  expect_true(is.na(outputDS2$matrix))
  expect_true(is.na(outputDS2$vector))
  expect_true(is.na(outputDS2$dv))
  expect_equal(outputDS2$disclosure.risk, 1)
  expect_length(outputDS2$errorMessage2, 3)
  expect_equal(outputDS2$errorMessage2[1], "ERROR: Model has too many parameters, there is a possible risk of disclosure - please simplify model")
  expect_equal(outputDS2$errorMessage2[2], "ERROR: at least one column in sigma.x matrix is binary with one category less than filter threshold for table cell size")
  expect_equal(outputDS2$errorMessage2[3], "MODEL FAILED: model or data invalid, matrix and vector destroyed")
  expect_equal(length(ls(pattern = "^temp_")), 0)
})



# Test-smk-gamlssDS2 ################################################

#
# Set up
#

set.standard.disclosure.settings()
set.random.seed.setting(379)

#
# Tests
#

test_that("output_gamlssDS2", {
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
  outputDS2 <- gamlssDS2(
    parameter = "mu", family = "BCPEleft_parenthesisright_parenthesis", data = "gamlss1",
    mu.beta.vect = paste0(as.character(c(0, 0)), collapse = ","), sigma.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    nu.beta.vect = paste0(as.character(c(0, 0)), collapse = ","), tau.beta.vect = paste0(as.character(c(0, 0)), collapse = ","),
    control = paste0(as.character(c(0.001, 20, 1, 1, 1, 1, Inf)), collapse = ","),
    i.control = paste0(as.character(c(0.001, 50, 30, 0.001)), collapse = ",")
  )
  eta <- BCPE()$mu.linkfun(temp_mu)
  dr <- 1 / BCPE()$mu.dr(eta)
  dldp <- BCPE()$dldm(gamlss1$e3_bw, temp_mu, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- BCPE()$d2ldm2(gamlss1$e3_bw, temp_mu, temp_sigma, temp_nu, temp_tau)
  d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
  wt <- -(d2ldp2 / (dr * dr))
  wt <- ifelse(wt > 1e+10, 1e+10, wt)
  wt <- ifelse(wt < 1e-10, 1e-10, wt)
  wv <- eta + dldp / (dr * wt)
  partial.residuals <- wv - base::rowSums(temp_mu.s)
  vector <- t(temp_mod.gamlss.ds$mu.x) %*% (wt * partial.residuals)
  matrix <- t(temp_mod.gamlss.ds$mu.x) %*% diag(wt) %*% temp_mod.gamlss.ds$mu.x

  expect_equal(outputDS2$matrix, matrix, ignore_attr = TRUE)
  expect_equal(outputDS2$vector, vector, ignore_attr = TRUE)
  expect_equal(outputDS2$dv, 8576.428, tolerance = 1e-07)
  expect_equal(outputDS2$disclosure.risk, 0)
  expect_equal(outputDS2$errorMessage2, "No errors")
})

#
# Done
#
