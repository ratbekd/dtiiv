# dtiiv

Distributional Tendency Intrinsic Instrumental Variables (DT-IIV) for R.

`dtiiv` estimates the slope of an endogenous regressor in a linear model
**without an external instrument**. It replaces the exclusion restriction
with a restriction on the joint distribution of the regressor's exogenous
part and the error, and uses the skewness and kurtosis of the data to find
the slope at which that restriction holds.

## Installation

```r
# from a local copy of this folder
install.packages("path/to/dtiiv", repos = NULL, type = "source")
```

## Example

```r
library(dtiiv)
d <- dtiiv_simulate(n = 3000, beta = 1, clusters = 100)

fit <- dtiiv(y ~ x, data = d,
             controls = ~ w,                 # exogenous controls
             cluster  = ~ cluster,           # cluster-robust inference
             sets     = c("VT", "mean", "error", "both"))
fit           # CUE estimates, J tests and confidence sets
summary(fit)  # plus shape diagnostics and the roots of the variance tendency
plot(fit)     # test statistic / critical value against the slope
confint(fit, set = "mean")
dtiiv_test(fit, beta0 = c(0, 1))
```

## Moment sets

| Set | Restriction | Identifies when |
|---|---|---|
| `"VT"` | variance tendency, E[u² s] = 0 | the error is skewed |
| `"mean"` | E[u₀ \| s₀] = 0 | the exogenous part of x is non-Gaussian |
| `"bounded"` | E[u₀ \| s₀] = 0, bounded functions of a robust z | as `"mean"`; recommended for heavy-tailed x |
| `"error"` | E[s₀ \| u₀] = 0 | the error is non-Gaussian |
| `"both"` | both | either is non-Gaussian |
| `"two-sided"`, `"one-sided"`, `"four"`, `"five"`, `"variance"` | polynomial sets of the paper | see `?dtiiv_moments` |

User-defined moments: `dtiiv_moment()`.

## Inference

* Confidence sets invert a score test with first-step corrections for the
  projections on the controls, the OLS slope and the scales. They are robust
  to weak identification and may be unions of intervals or empty.
* Variances are cluster-robust. The default critical values are
  Hotelling-type F values with a Satterthwaite estimate of the cluster degrees
  of freedom; `crit = "chisq"` gives chi-squared values.
* Point estimates are continuously updated GMM with a sandwich standard error,
  and the J test checks the overidentifying restrictions.

The computations reproduce the paper's replication code (Python core in
`DT_IIV_replication/code/lib/iiv/core.py`) to numerical precision.

## Limits

* One endogenous regressor.
* Both the regressor and the outcome binary: not identified (the function
  stops). A binary regressor alone: the mean tendencies hold mechanically at
  OLS.
* No sampling weights yet.
