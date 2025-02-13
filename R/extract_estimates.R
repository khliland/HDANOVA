#' @title Extract estimates for a given factor combination
#'
#' @description Extracts and sums the LS estimates for a given factor combination
#' from an object of class \code{hdanova}. If \code{add_residuals} is \code{TRUE},
#' the residuals are added to the LS estimates. If \code{remove_factors} is \code{TRUE},
#' the returned matrix is the data with chosen estimates subtracted.
#'
#' @param object \code{asca} object.
#' @param factors \code{vector} of factor names or numbers.
#' @param subtract \code{logical} subtract the estimates from the data (default = FALSE).
#' @param add_residuals \code{logical} add residuals to the estimates (default = FALSE).
#'
#' @returns A matrix of the extracted estimates.
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic HDANOVA model with two factors and interaction
#' mod <- hdanova(assessment ~ candy * assessor, data=candies)
#'
#' # Extract estimates for the interaction
#' inter <- extract_estimates(mod, c("assessor:candy"))
#'
#' # Visualize the interaction effect
#' image(t(inter), main="Interaction effect", xlab="Attribute", ylab="Sample")
#'
#' @export
extract_estimates <- function(object, factors, subtract=FALSE, add_residuals=FALSE){
  if(object$add_error && length(factors) > 1)
    stop("If errors are added to LS means, e.g., in APCA, only one factor can be used.")
  estimate <- do.call("+", object$LS[factors])
  if(add_residuals){
    estimate <- estimate + object$residuals
  }
  if(subtract){
    estimate <- object$Y - estimate
  }
  return(estimate)
}
