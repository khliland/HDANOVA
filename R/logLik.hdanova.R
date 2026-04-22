#' Log-Likelihood for HDANOVA Objects
#'
#' @description Extract log-likelihood from fitted HDANOVA objects. For fixed-effects
#' and MoM ('r()' with \code{REML = NULL}) models, the log-likelihood is computed
#' from the residual sum of squares. For REML/ML mixed models ('r()' with
#' \code{REML = TRUE} or \code{FALSE}), the log-likelihood is extracted directly
#' from the stored lme4 model objects.
#'
#' @param object A fitted \code{hdanova} object.
#' @param ... Reserved for generic compatibility.
#'
#' @details The log-likelihood is computed as follows:
#'
#' \strong{Fixed-effects and MoM models:} The multivariate Gaussian log-likelihood is
#' \deqn{\ell_j = -\frac{n}{2}\log(2\pi) - \frac{n}{2}\log(\sigma_j^2) - \frac{1}{2\sigma_j^2}\sum_i e_{ij}^2}
#' where \eqn{e_{ij}} is the residual for observation \eqn{i} and response \eqn{j},
#' and \eqn{\sigma_j^2 = SSE_j / n} is the ML variance estimate for response
#' column \eqn{j}. The total log-likelihood is \eqn{\ell = \sum_j \ell_j}.
#'
#' \strong{REML/ML mixed models:} The log-likelihood is the sum of the individual
#' response-specific log-likelihoods from the fitted lme4 models,
#' \deqn{\ell = \sum_j \ell_j}
#' where \eqn{\ell_j} is the log-likelihood of the lme4::lmerMod fit for response \eqn{j}.
#'
#' @return An object of class \code{logLik} with the computed log-likelihood value,
#' degrees of freedom ('df'), and the number of observations ('nobs').
#'
#' @seealso Model fitting: \code{\link{hdanova}}.
#' Information criteria: \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}}.
#'
#' @examples
#' data(candies)
#' mod <- hdanova(assessment ~ candy + assessor, data = candies)
#' ll <- logLik(mod)
#' print(ll)
#'
#' \dontrun{
#' # For mixed models:
#' mod_reml <- hdanova(assessment ~ candy + r(assessor), data = candies, REML = TRUE)
#' ll_reml <- logLik(mod_reml)
#' print(ll_reml)
#' }
#'
#' @method logLik hdanova
#' @export
logLik.hdanova <- function(object, ...){
  more <- object$more
  is_reml_or_ml <- is.logical(more$REML)

  if(!is_reml_or_ml){
    # Fixed-effects or MoM model: compute from SSE/residuals
    Y <- object$Y
    residuals <- object$residuals
    n <- nrow(residuals)
    p <- ncol(residuals)
    
    # Number of parameters per response
    # For a univariate model: intercept + fixed effects
    # We use the dimension of coefficients
    coefs <- object$coefficients
    if(is.matrix(coefs)){
      n_params <- nrow(coefs)
    } else {
      n_params <- length(coefs)
    }
    
    # Residual degrees of freedom for each response
    df_resid_per_response <- n - n_params
    
    # Sum of squared residuals per response
    SSE <- colSums(residuals^2)
    
    # Maximum likelihood estimator of variance (dividing by n, not n-p)
    sigma2_ml <- SSE / n
    
    # Multivariate Gaussian log-likelihood for each response
    ll_per_response <- -n/2 * log(2*pi) - n/2 * log(sigma2_ml) - n/2
    
    # Total log-likelihood: sum across responses
    ll_total <- sum(ll_per_response)
    
    # Degrees of freedom: total number of fixed-effect parameters
    df_total <- n_params * p
    
    # Total observations across all responses
    nobs <- n * p
    
  } else {
    # REML/ML mixed model: extract from lme4 objects
    models <- object$models
    if(!is.list(models) || length(models) == 0){
      stop("No fitted lme4 models found in object$models.")
    }
    
    # Sum the log-likelihood from each response-specific model
    ll_per_response <- sapply(models, function(m){
      as.numeric(stats::logLik(m))
    })
    ll_total <- sum(ll_per_response)
    
    # Degrees of freedom: number of parameters in first model (all should be same)
    # For lmerMod, df includes fixed effects and variance parameters
    df_total <- 0
    for(m in models){
      # Count fixed effects
      df_total <- df_total + length(lme4::fixef(m))
      # Count variance/covariance parameters
      # This includes residual variance and random effect variances
      vc <- lme4::VarCorr(m)
      # Number of variance components: 1 for residual + number of random effect groups
      df_total <- df_total + 1 + length(vc) - 1
    }
    
    # Total observations
    nobs <- nrow(object$model.frame) * ncol(object$Y)
  }

  # Create logLik object
  attr(ll_total, "df") <- df_total
  attr(ll_total, "nobs") <- nobs
  class(ll_total) <- "logLik"
  
  ll_total
}

