# **Bergm:** Bayesian Exponential Random Graph Models

**Bergm** provides a comprehensive framework for Bayesian parameter estimation
and model selection for exponential random graph models using advanced
computational algorithms. It can also supply graphical Bayesian goodness-of-fit
procedures that address the issue of model adequacy and missing data
imputation.

## What's new

**Bergm 6** is a rewrite of the computational core that keeps the **Bergm 5** interface
(same functions, arguments and returned elements) while making every function
faster and more flexible. See `NEWS.md` for the full list.

* **Faster.** Auxiliary networks are drawn by single calls to ergm's compiled
  sampler on pre-built network states; pseudolikelihood-based methods
  (`bergmC`, `ergmAPL`, `evidenceCJ`, `evidencePP`) run their MCMC and
  pseudolikelihood evaluations in compiled code; `bgof` monitors
  goodness-of-fit statistics inside the sampler.
* **More efficient samplers.** Burn-in adaptation of the ADS move factor,
  optional delayed acceptance that skips network simulation for proposals a
  cheap surrogate already rejects, parallel population updates (`cores`), and a
  Newton-type stochastic approximation for the posterior mode in `bergmC`.
* **More flexible.** Curved terms (`gwesp(fixed = FALSE)`), one-parameter
  models, user-defined priors (`log.prior`), network proposals via `...`,
  reusable pseudolikelihood adjustments (`apl =`), posterior predictive
  simulation (`simulate()`), per-chain output with R-hat, and `coef()`,
  `vcov()`, `as.mcmc()` methods.

### Speed-ups over Bergm 5.0.7

Wall-clock time on a 2-core machine (ergm 4.12.0), average of two runs; "ESS/s
gain" compares effective samples per second for the worst-mixing parameter.
Full table, settings and validation against Bergm 5 in `benchmarks/` (in the
source repository).

| Function, data, (Bergm 6 option) | Bergm 5 | Bergm 6 | ESS/s gain |
|---|---:|---:|---:|
| `bergm()`, Florentine | 9.4 s | 3.6 s | 2.7x |
| `bergm()`, Lazega | 45.4 s | 31.9 s | 1.4x |
| `bergm()`, Lazega, (delayed acceptance) | 45.4 s | 18.2 s | 1.9x |
| `bergm()`, Faux Mesa High, (delayed acceptance) | 44.2 s | 19.1 s | 2.7x |
| `bergmM()`, Florentine with missing nodes | 62.2 s | 5.4 s | 10.0x |
| `bergmC()`, Florentine | 17.9 s | 1.0 s | 14.8x |
| `bergmC()`, Lazega | 33.1 s | 6.0 s | 5.1x |
| `evidencePP()`, Florentine | 18.1 s | 2.2 s | 8.1x |
| `evidencePP()`, Lazega | 77.2 s | 13.3 s | 5.8x |
| `bgof()`, Florentine | 3.9 s | 0.2 s | |
| `bgof()`, Lazega | 7.5 s | 2.0 s | |

```r
library(Bergm)
data(florentine)

# Curved ERGM, delayed acceptance (default) and 4 workers
data(lazega)

fit <- bergm(lazega ~ edges + nodematch("Office") + gwesp(0.5, fixed = FALSE),
              burn.in = 200, main.iters = 2000, aux.iters = 3000,
              seed = 1)

# Goodness of fit and posterior predictive networks
bgof(fit, sample.size = 100, aux.iters = 5000)
nets <- simulate(fit, nsim = 10)

# Model selection
a1 <- ergmAPL(flomarriage ~ edges + kstar(2))
a2 <- ergmAPL(flomarriage ~ edges + triangle)
evidencePP(flomarriage ~ edges + kstar(2), apl = a1)$log.evidence -
  evidencePP(flomarriage ~ edges + triangle, apl = a2)$log.evidence
```

### How to cite **Bergm**

Caimo, A., Bouranis, L., Krause, R., and Friel, N. (2022). Statistical Network Analysis with Bergm. *Journal of Statistical Software*, 104(1), 1–23. doi: https://doi.org/10.18637/jss.v104.i01.
