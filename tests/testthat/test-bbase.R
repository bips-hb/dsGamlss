# Test-smk-bbase ################################################

#
# Set up
#

#
# Tests
#

test_that("output_bbase", {
  x <- seq(1, 100, by = 1)
  xl <- min(x) - 0.01 * (max(x) - min(x))
  xr <- max(x) + 0.01 * (max(x) - min(x))
  bbase.ds <- bbase(x = x, xl = xl, xr = xr, ndx = 20, deg = 3)
  bbase <- attr(gamlss::pb(x), "X")
  expect_equal(bbase.ds, bbase, tolerance = 1e-07)
})

#
# Done
#
