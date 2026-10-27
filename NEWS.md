# Bergm 6.0.0

A performance and flexibility release. The user-facing API of Bergm 5 is
preserved (same functions, arguments and returned elements); new behaviour is
opt-in through new arguments, except where noted under "Changed defaults".

## Speed

* `bergm()` / `bergmM()`: every auxiliary network is drawn by one call to
  **ergm** compiled sampler on a network state built once per fit, instead of
  going through `ergm_MCMC_sample()`; the Gaussian prior uses a cached
  Cholesky factor and proposal noise is drawn in blocks.
* `bergmC()`: the pseudo-posterior is sampled by a compiled random-walk
  Metropolis kernel,
  and the posterior mode is located by a Newton-type stochastic approximation
  with a warm-started auxiliary chain (`map.method = "newton"`, default), which
  converges in about ten iterations (instead of hundreds of Robbins-Monro steps).
* `evidenceCJ()` / `evidencePP()`: compiled sampling of the (tempered)
  adjusted pseudo-posteriors and compiled, vectorised evaluation of the
  adjusted pseudolikelihood and its gradient for the estimators and the
  zero-variance control variates.
* `ergmAPL()`: bridge sampling from the dyad-independent submodel needs far
  fewer simulated networks.
* `bgof()`: goodness-of-fit statistics are monitored inside **ergm** sampler
  with optional
  parallel simulation (`cores`).
* `bergmM()`: imputation steps act directly on compiled network states
  instead of rebuilding networks and models with `simulate()`.

## New features

* `bergm()`: parallel chains (`cores`), burn-in adaptation of the ADS move
  factor (`adapt`, `target.accept`), delayed acceptance with a Gaussian
  surrogate fitted on the burn-in population (`delayed.acceptance`,
  `da.inflate`), user-defined priors (`log.prior`, e.g. truncated or
  heavy-tailed priors), reproducible runs (`seed`), and `...` forwarded to
  `control.simulate.formula()` to choose the network proposal.
* Curved ERGM terms (e.g. `gwesp(fixed = FALSE)`) and one-parameter models are
  now supported by `bergm()` and `bergmM()`.
* Decay parameters of curved terms (`*.decay`, `altkstar.lambda`) are sampled
  on the log scale, so they are always positive, with a log-normal prior set by
  the new argument `prior.decay` (default `meanlog = log(0.5)`, `sdlog = 0.75`).
  `prior.mean` and `prior.sigma` apply to the other parameters and can be
  given for all parameters (decay entries ignored) or for the non-decay ones
  only. A custom `log.prior` receives the decay on its natural scale and the
  Jacobian of the log transformation is added automatically.
* `bergm()`/`bergmM()` set ergm's `gw.cutoff` term option to `n - 2` (at most
  1000) unless it is given through `term.options`, so geometrically weighted
  statistics are exact for every simulated network. With the default cut-off
  of 30, a proposal in a near-degenerate region could make the auxiliary
  simulation stop with an error, and a chain that reached such a region could
  not leave it. Any remaining simulation that exceeds an ergm limit counts as
  a rejected proposal; the number is reported in `$n.failed`, with a warning.
* Fits keep the per-chain samples (`$chains`, an `mcmc.list`), per-chain
  acceptance rates, potential scale reduction factors (`$Rhat`, also shown by
  `summary()`) and the number of simulated networks (`$n.sims`).
* New methods: `print()`, `coef()`, `vcov()`, `as.mcmc()`, `as.mcmc.list()` and
  `simulate()` (posterior predictive networks or statistics) for `bergm`
  objects; `print()` and `plot()` for `bgof` objects; `print()` for `ergmAPL`
  objects.
* `ergmAPL()` returns an `ergmAPL` object that can be reused by
  `evidenceCJ()` / `evidencePP()` through their new `apl` argument, e.g. to
  compare evidence estimators without recomputing the adjustment. It also
  accepts a precomputed MLE (`mle`).
* `bgof()` additionally returns the posterior predictive distribution of the
  model statistics with p-values, and plotting can be turned off (`plot = FALSE`).
* `prior.sigma` and `V.proposal` accept a scalar, a vector (diagonal) or a matrix.
* `plot()` and `bgof()` draw in grayscale: chains in shades of gray, observed
  statistics as a thick black line, 5%/95% posterior predictive bands dashed.

* `bergm(cores > 1)` with delayed acceptance: the surrogate's mean was
  overwritten by a loop index after the first parallel update, so almost every
  proposal was rejected at the first stage and the chains stopped moving.
* `bergm(cores > 1)`: distributing a single proposal to the workers failed with
  "invalid number of intervals".
* `bergmM()`: missing ties in bipartite networks were mapped to the wrong
  dyads (and the attribute-imputation predictors used the wrong matrix).
* `bergm()`/`bergmM()`: `saveEveryX` larger than the run length no longer
  fails; `main.iters < thin` is rejected with a clear message.
* `evidenceCJ()`: the posterior ordinate is computed on the log scale, which
  avoids underflow of the proposal densities in higher dimensions.
* `bgof()`/`simulate()`: posterior predictive simulation uses the same
  `gw.cutoff` as the fit; bipartite networks give a clear error in `bgof()`;
  failed simulations report the actual ergm limit that was exceeded.

## Changed defaults

* `bergm()`: `delayed.acceptance = TRUE`. With a burn-in shorter than 10
  iterations it is switched off silently; set `delayed.acceptance = FALSE`
  for the Bergm 5 kernel.
* `bergm()`/`bergmM()`: `adapt = TRUE` (the move factor `gamma` is tuned during
  burn-in only), chains start overdispersed around the MPLE using its
  covariance, `cut.reject = FALSE` and, for `bergmM()`, `imputeAllItr = TRUE`.
  Set `adapt = FALSE` to reproduce the Bergm 5 sampler kernel.
* `bergmC()`: `map.method = "newton"`; use `map.method = "rm"` for the Bergm 5
  Robbins-Monro recursion.

## Dependencies

* `MCMCpack`, `mvtnorm`, `matrixcalc`, `Matrix`, `Rglpk`, `statnet.common` and
  `stringr` are no longer needed; `Rcpp` is used for the compiled kernels.
