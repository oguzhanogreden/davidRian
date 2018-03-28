# Input tests
# These rather redundant but I find it safer to keep an eye on Rcpp <-> R interface.

test_that("Wrong inputs fail ddc().", {
  expect_error(ddc(1))                # missing phi
  expect_error(ddc(1, rep(1.5, 11)))  # too many phis
  expect_error(ddc(1, "a"))           # character phi
  expect_error(ddc("a", 1.5))         # character x
})

test_that("Wrong inputs fail rdc().", {
  expect_error(rdc(1))                # missing phi
  expect_error(rdc(1, rep(1.5, 11)))  # too many phis
  expect_error(rdc(1, "a"))           # character phi
  expect_error(rdc("a", 1.5))         # character x
})

test_that("Wrong inputs fail dcGrad().", {
  expect_error(dcGrad(1))                # missing phi
  expect_error(dcGrad(1, rep(1.5, 11)))  # too many phis
  expect_error(dcGrad(1, "a"))           # character phi
  expect_error(dcGrad("a", 1.5))         # character x
})

test_that("Infinities returns 0 in ddc().", {
  expect_that(ddc(Inf, 1.5), equals(0))   # Inf should be 0.
  expect_that(ddc(-Inf, 1.5), equals(0))  # -Inf should be 0.
})