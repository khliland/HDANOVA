# Log-Likelihood for HDANOVA Objects

Extract log-likelihood from fitted HDANOVA objects. For fixed-effects
and MoM ('r()' with `REML = NULL`) models, the log-likelihood is
computed from the residual sum of squares. For REML/ML mixed models
('r()' with `REML = TRUE` or `FALSE`), the log-likelihood is extracted
directly from the stored lme4 model objects.

## Usage

``` r
# S3 method for class 'hdanova'
logLik(object, ...)
```

## Arguments

- object:

  A fitted `hdanova` object.

- ...:

  Reserved for generic compatibility.

## Value

An object of class `logLik` with the computed log-likelihood value,
degrees of freedom ('df'), and the number of observations ('nobs').

## Details

The log-likelihood is computed as follows:

**Fixed-effects and MoM models:** The multivariate Gaussian
log-likelihood is \$\$\ell_j = -\frac{n}{2}\log(2\pi) -
\frac{n}{2}\log(\sigma_j^2) - \frac{1}{2\sigma_j^2}\sum_i e\_{ij}^2\$\$
where \\e\_{ij}\\ is the residual for observation \\i\\ and response
\\j\\, and \\\sigma_j^2 = SSE_j / n\\ is the ML variance estimate for
response column \\j\\. The total log-likelihood is \\\ell = \sum_j
\ell_j\\.

**REML/ML mixed models:** The log-likelihood is the sum of the
individual response-specific log-likelihoods from the fitted lme4
models, \$\$\ell = \sum_j \ell_j\$\$ where \\\ell_j\\ is the
log-likelihood of the lme4::lmerMod fit for response \\j\\.

## See also

Model fitting:
[`hdanova`](https://khliland.github.io/HDANOVA/reference/hdanova.md).
Information criteria: [`AIC`](https://rdrr.io/r/stats/AIC.html) and
[`BIC`](https://rdrr.io/r/stats/AIC.html).

## Examples

``` r
data(candies)
mod <- hdanova(assessment ~ candy + assessor, data = candies)
ll <- logLik(mod)
print(ll)
#> 'log Lik.' -3388.011 (df=135)

if (FALSE) { # \dontrun{
# For mixed models:
mod_reml <- hdanova(assessment ~ candy + r(assessor), data = candies, REML = TRUE)
ll_reml <- logLik(mod_reml)
print(ll_reml)
} # }
```
