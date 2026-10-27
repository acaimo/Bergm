test_that("bergmM imputes missing ties per chain", {
  data(florentine, package = "ergm")
  set.seed(1)
  flomarriage[3, ] <- NA
  flomarriage[, 3] <- NA
  m <- bergmM(flomarriage ~ edges + kstar(2), burn.in = 10, main.iters = 30,
              aux.iters = 100, nImp = 3, seed = 2, verbose = FALSE)
  expect_s3_class(m, "bergm")
  expect_length(m$impNets, 3)
  am <- as.matrix(m$impNets[[1]])
  obs <- as.matrix(flomarriage)
  ok <- !is.na(obs)
  expect_equal(am[ok], obs[ok])
  expect_false(anyNA(am))
})

test_that("bgof returns distributions with the expected shapes", {
  data(florentine, package = "ergm")
  a <- bergm(flomarriage ~ edges + kstar(2), burn.in = 10, main.iters = 30,
             aux.iters = 100, seed = 1, verbose = FALSE)
  g <- bgof(a, sample.size = 8, aux.iters = 200, plot = FALSE, seed = 1)
  n <- network.size(flomarriage)
  expect_s3_class(g, "bgof")
  expect_equal(dim(g$sim.degree), c(n, 8))
  expect_equal(nrow(g$sim.esp), n - 1)
  expect_equal(unname(colSums(g$sim.degree)), rep(1, 8))
  expect_equal(unname(g$obs.model), c(20, 47))
  s <- simulate(a, nsim = 3, aux.iters = 200, output = "stats", seed = 1)
  expect_equal(dim(s), c(3, 2))
  nets <- simulate(a, nsim = 2, aux.iters = 200, seed = 1)
  expect_s3_class(nets[[1]], "network")
})

test_that("bergmM imputes missing nodal attributes", {
  skip_on_cran()
  data(lazega, package = "Bergm")
  y <- lazega
  pr <- y %v% "Practice"
  pr[c(2, 9)] <- NA
  set.vertex.attribute(y, "Practice", pr)
  imputeData <- data.frame(Practice = factor(pr), Office = factor(y %v% "Office"),
                           Age = y %v% "Age")
  m <- bergmM(y ~ edges + nodematch("Practice"), burn.in = 2, main.iters = 6,
              aux.iters = 200, imputeData = imputeData, attributeNames = "Practice",
              nImp = 2, seed = 1, verbose = FALSE)
  expect_length(m$impAttr, 2)
  expect_false(anyNA(m$impAttr[[1]]$Practice))
})

test_that("bergmM imputes the right dyads in a bipartite network", {
  skip_if_not_installed("latentnet")
  data(davis, package = "latentnet")
  y <- davis
  obs <- as.matrix(y, matrix.type = "adjacency")
  y[1, 19] <- NA
  y[5, 25] <- NA
  m <- bergmM(y ~ edges + gwb2degree(0.7, fixed = TRUE), burn.in = 10, main.iters = 30,
              aux.iters = 200, nImp = 2, seed = 1, verbose = FALSE)
  imp <- as.matrix(m$impNets[[1]], matrix.type = "adjacency")
  keep <- !is.na(as.matrix(y, matrix.type = "adjacency"))
  expect_equal(imp[keep], obs[keep])
  expect_equal(dim(imp), dim(obs))
})
