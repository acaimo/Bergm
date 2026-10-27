test_that("bergm returns a well-formed object and is reproducible", {
  data(florentine, package = "ergm")
  f <- flomarriage ~ edges + kstar(2)
  a <- bergm(f, burn.in = 20, main.iters = 60, aux.iters = 200, seed = 11, verbose = FALSE)
  b <- bergm(f, burn.in = 20, main.iters = 60, aux.iters = 200, seed = 11, verbose = FALSE)
  expect_s3_class(a, "bergm")
  expect_equal(as.matrix(a$Theta), as.matrix(b$Theta))
  expect_equal(dim(a$Theta), c(60 * a$nchains, 2))
  expect_equal(length(a$chains), a$nchains)
  expect_named(coef(a), c("edges", "kstar2"))
  expect_equal(dim(vcov(a)), c(2, 2))
  expect_true(all(is.finite(a$ess)))
  out <- capture.output(tab <- summary(a, returnTable = TRUE))
  expect_true(all(c("Mean", "SD", "ESS", "Rhat") %in% colnames(tab)))
})

test_that("thinning, offsets, single-parameter models and custom priors work", {
  data(florentine, package = "ergm")
  a <- bergm(flomarriage ~ edges + kstar(2), burn.in = 10, main.iters = 40, thin = 4,
             aux.iters = 100, seed = 1, verbose = FALSE)
  expect_equal(nrow(a$Theta), 10 * a$nchains)

  o <- bergm(flomarriage ~ edges + kstar(2) + offset(triangle), offset.coef = -0.5,
             burn.in = 10, main.iters = 30, aux.iters = 100, seed = 1, verbose = FALSE)
  expect_true(all(o$Theta[, 3] == -0.5))
  expect_equal(unname(o$fixed), c(FALSE, FALSE, TRUE))
  expect_output(summary(o), "fixed")

  e <- bergm(flomarriage ~ edges, burn.in = 10, main.iters = 30, aux.iters = 100,
             seed = 1, verbose = FALSE)
  expect_equal(ncol(e$Theta), 1)

  tp <- bergm(flomarriage ~ edges + kstar(2), burn.in = 10, main.iters = 30, aux.iters = 100,
              seed = 1, verbose = FALSE,
              log.prior = function(th) if (th["kstar2"] > 0) -Inf else 0)
  expect_true(all(tp$Theta[, "kstar2"] <= 0))
})

test_that("delayed acceptance runs and simulates fewer networks", {
  data(florentine, package = "ergm")
  d <- bergm(flomarriage ~ edges + kstar(2), burn.in = 40, main.iters = 100,
             aux.iters = 100, seed = 3, delayed.acceptance = TRUE, verbose = FALSE)
  expect_false(is.null(d$delayed.acceptance))
  expect_lt(d$n.sims[["main"]], 100 * d$nchains)
})

test_that("curved terms are supported", {
  data(florentine, package = "ergm")
  cfit <- suppressWarnings(
    bergm(flomarriage ~ edges + gwesp(fixed = FALSE), burn.in = 10,
          main.iters = 20, aux.iters = 200, seed = 1, verbose = FALSE))
  expect_equal(ncol(cfit$Theta), 3)
})

test_that("a run can be continued from a previous fit", {
  data(florentine, package = "ergm")
  a <- bergm(flomarriage ~ edges + kstar(2), burn.in = 10, main.iters = 20,
             aux.iters = 100, seed = 1, verbose = FALSE)
  b <- bergm(flomarriage ~ edges + kstar(2), burn.in = 0, main.iters = 20,
             aux.iters = 100, seed = 2, startVals = a, verbose = FALSE)
  expect_equal(b$nchains, a$nchains)
  expect_equal(b$gamma, a$gamma)
  expect_error(bergm(flomarriage ~ edges, startVals = a, main.iters = 5, verbose = FALSE),
               "different model")
})

test_that("decay parameters are sampled positive under a log-normal prior", {
  data(florentine, package = "ergm")
  f <- flomarriage ~ edges + gwesp(fixed = FALSE)
  a <- suppressWarnings(bergm(f, burn.in = 20, main.iters = 60, aux.iters = 200,
                              seed = 1, verbose = FALSE))
  expect_true(all(a$Theta[, "gwesp.decay"] > 0))
  expect_equal(unname(a$prior$decay$meanlog), log(0.5))
  expect_true(is.na(a$prior$mean[["gwesp.decay"]]))

  b <- suppressWarnings(bergm(f, burn.in = 20, main.iters = 60, aux.iters = 200, seed = 1,
                              prior.mean = c(-2, 0), prior.sigma = c(4, 4),
                              prior.decay = list(meanlog = 0, sdlog = 0.3),
                              verbose = FALSE))
  expect_equal(unname(b$prior$mean[1:2]), c(-2, 0))
  expect_lt(abs(median(log(b$Theta[, "gwesp.decay"]))), 1)

  expect_error(bergm(f, startVals = c(-2, 0, -1), main.iters = 5, verbose = FALSE),
               "positive")
})

test_that("the log-normal decay prior is exact on the log scale", {
  data(florentine, package = "ergm")
  setup <- Bergm:::.ergm_setup(flomarriage ~ edges + gwesp(fixed = FALSE))
  dec <- Bergm:::.decay_index(setup)
  expect_equal(unname(dec), c(FALSE, FALSE, TRUE))
  pr <- Bergm:::.make_prior(NULL, NULL, 3, rep(TRUE, 3), setup$coef.names, dec = dec)
  phi <- c(-1, 0.2, log(0.8))
  expect_equal(pr$logdens(phi),
               sum(dnorm(phi[1:2], 0, 10, log = TRUE)) +
                 dlnorm(0.8, log(0.5), 0.75, log = TRUE) + log(0.8))
  lp <- function(th) dexp(th[["gwesp.decay"]], 1, log = TRUE)
  pc <- Bergm:::.make_prior(NULL, NULL, 3, rep(TRUE, 3), setup$coef.names,
                            log.prior = lp, dec = dec)
  expect_equal(pc$logdens(phi), dexp(0.8, 1, log = TRUE) + log(0.8))
})

test_that("an unidentified decay parameter recovers its prior", {
  skip_on_cran()
  data(florentine, package = "ergm")
  f <- flomarriage ~ edges + gwesp(fixed = FALSE)
  fit <- function(...) suppressWarnings(
    bergm(f, prior.mean = c(-1.6, 0), prior.sigma = c(1, 1e-6), burn.in = 300,
          main.iters = 2000, aux.iters = 100, nchains = 8, seed = 1,
          V.proposal = c(0.0025, 1e-9, 0.01), verbose = FALSE, ...))
  a <- fit(prior.decay = list(meanlog = 0.5, sdlog = 0.3))
  lt <- log(a$Theta[, "gwesp.decay"])
  expect_lt(abs(mean(lt) - 0.5), 0.05)
  expect_lt(abs(sd(lt) / 0.3 - 1), 0.15)
  lp <- function(th) dnorm(th[["edges"]], -1.6, 1, log = TRUE) +
    dnorm(th[["gwesp"]], 0, 1e-3, log = TRUE) +
    dgamma(th[["gwesp.decay"]], shape = 3, rate = 4, log = TRUE)
  g <- fit(log.prior = lp)$Theta[, "gwesp.decay"]
  expect_lt(abs(mean(g) / 0.75 - 1), 0.1)
  expect_lt(abs(sd(g) / sqrt(3 / 16) - 1), 0.15)
})
