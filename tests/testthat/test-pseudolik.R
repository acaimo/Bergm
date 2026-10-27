test_that("compiled pseudolikelihood matches the binomial formula", {
  data(florentine, package = "ergm")
  pl <- Bergm:::.pl_setup(flomarriage ~ edges + kstar(2))
  th <- c(-1.3, 0.05)
  eta <- drop(pl$X %*% th)
  ref <- sum(dbinom(pl$w * pl$y, pl$w, plogis(eta), log = TRUE))
  expect_equal(Bergm:::.pl_loglik(pl, th), ref, tolerance = 1e-10)

  W <- matrix(c(1.1, 0.2, 0, 0.9), 2)
  adj <- list(W = W, Theta_MLE = c(-1.6, 0), Theta_PL = c(-1.65, 0.01))
  Th <- rbind(th, c(-1, -0.1))
  b <- Bergm:::.adjpl_batch(pl, Th, adj, grad = TRUE)
  tr <- function(t) drop(W %*% (t - adj$Theta_MLE) + adj$Theta_PL)
  expect_equal(b$logpl[1], Bergm:::.pl_loglik(pl, tr(th)), tolerance = 1e-10)
  num <- sapply(1:2, function(k) {
    e <- c(0, 0); e[k] <- 1e-6
    (Bergm:::.pl_loglik(pl, tr(th + e)) - Bergm:::.pl_loglik(pl, tr(th - e))) / 2e-6
  })
  expect_equal(b$grad[1, ], num, tolerance = 1e-5)
})

test_that("pseudo-posterior mode is a stationary point", {
  data(florentine, package = "ergm")
  pl <- Bergm:::.pl_setup(flomarriage ~ edges + kstar(2))
  pr <- Bergm:::.gauss_prior(NULL, NULL, 2, pl$names)
  m <- Bergm:::.pl_mode(pl, prior = pr)
  g <- Bergm:::.pl_grad(pl, m) - drop(pr$prec %*% (m - pr$mean))
  expect_lt(max(abs(g)), 1e-6)
})

test_that("compiled RWM targets the right distribution (Gaussian check)", {
  data(florentine, package = "ergm")
  pl <- Bergm:::.pl_setup(flomarriage ~ edges)
  pr <- Bergm:::.gauss_prior(0, 100, 1, pl$names)
  set.seed(1)
  r <- Bergm:::.rwm_pl(pl, pr, matrix(0.1), init = -1.5, burn.in = 1000, main.iters = 40000)
  grid <- seq(-3, 0, length.out = 4001)
  lp <- sapply(grid, function(t) Bergm:::.pl_loglik(pl, t)) + dnorm(grid, 0, 10, log = TRUE)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  expect_equal(mean(r$draws), sum(grid * w), tolerance = 0.02)
  expect_equal(sd(r$draws), sqrt(sum(grid^2 * w) - sum(grid * w)^2), tolerance = 0.05)
})
