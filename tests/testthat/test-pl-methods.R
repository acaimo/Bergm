test_that("ergmAPL, bergmC and evidence estimators run and agree roughly", {
  skip_on_cran()
  data(florentine, package = "ergm")
  f <- flomarriage ~ edges + kstar(2)
  apl <- ergmAPL(f, aux.iters = 500, n.aux.draws = 50, ladder = 20, seed = 1)
  expect_s3_class(apl, "ergmAPL")
  expect_true(is.finite(apl$logC))
  expect_equal(dim(apl$W), c(2, 2))

  cj <- evidenceCJ(f, main.iters = 3000, burn.in = 500, V.proposal = 2.5, apl = apl,
                   verbose = FALSE)
  pp <- evidencePP(f, main.iters = 2000, burn.in = 200, V.proposal = 2.5, apl = apl,
                   verbose = FALSE)
  expect_true(is.finite(cj$log.evidence))
  expect_true(is.finite(pp$log.evidence))
  expect_lt(abs(cj$log.evidence - pp$log.evidence), 1.5)
  expect_output(summary(cj), "Log-evidence")
  ev <- evidence("PP", formula = f, main.iters = 500, burn.in = 100, apl = apl, verbose = FALSE)
  expect_true(is.finite(ev$log.evidence))

  cc <- bergmC(f, burn.in = 500, main.iters = 2000, aux.iters = 500, n.aux.draws = 100,
               map.iters = 6, V.proposal = 2.5, init = apl$Theta_MLE, verbose = FALSE)
  expect_s3_class(cc, "bergm")
  expect_output(summary(cc), "Posterior Density")
  expect_equal(unname(colMeans(cc$Theta)), unname(cc$Theta_star), tolerance = 0.5)
})

test_that("ergmAPL recovers the exact likelihood of a dyad-independent model", {
  skip_on_cran()
  data(florentine, package = "ergm")
  f <- flomarriage ~ edges + absdiff("wealth")
  th <- coef(ergm(f))
  m <- ergmMPLE(f, output = "matrix")
  eta <- drop(m$predictor %*% th)
  exact <- sum(m$weights * (m$response * eta - log1p(exp(eta))))
  a <- ergmAPL(f, mle = th, seed = 1)
  expect_equal(a$ll_true, exact, tolerance = 1e-6)
  expect_lt(abs(a$logC), 1e-6)
})
