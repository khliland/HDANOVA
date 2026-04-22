#' Predict Methods for ASCA-family Objects
#'
#' Reconstructs the HDANOVA base on \code{newdata} via \code{predict.hdanova()}
#' and then applies the class-specific decomposition step:
#' \code{sca()} for ASCA/APCA/MSCA and \code{pls()} for APLS.
#' By default, decomposition is done by projection onto training component
#' spaces. A refit mode is also available.
#'
#' @param object A fitted \code{asca}, \code{apca}, \code{msca}, \code{apls}, or
#' \code{limmpca} object.
#' @param newdata A data frame containing variables from the original model formula.
#' @param decomposition Decomposition mode: \code{"project"} (default) projects
#'   onto training component spaces; \code{"refit"} recomputes decomposition on
#'   predicted LS matrices.
#' @param ... Reserved for compatibility; forwarded to \code{predict.hdanova()}.
#'
#' @return A predicted object of the same high-level class as \code{object}.
#'
#' @seealso Base prediction engine: \code{\link{predict.hdanova}}.
#' Related model constructors: \code{\link{asca}}, \code{\link{apca}},
#' \code{\link{apls}}, \code{\link{msca}} and \code{\link{limmpca}}.
#'
#' @examples
#' data(candies)
#' test_idx  <- seq(3, nrow(candies), by = 3)
#' train_idx <- setdiff(seq_len(nrow(candies)), test_idx)
#' candies_train <- candies[train_idx, ]
#' candies_test  <- candies[test_idx, ]
#'
#' mod_asca <- asca(assessment ~ candy * assessor, data = candies_train)
#' pred_asca <- predict(mod_asca, newdata = candies_test)
#' scoreplot(mod_asca, factor="candy", legend=TRUE)
#' with(pred_asca$projected, points(candy[,1], candy[,2], pch="x", cex=0.8,
#'                                  col=as.numeric(pred_asca$model.frame$candy)))
#'
#' pred_asca_refit <- predict(mod_asca, newdata = candies_test, decomposition = "refit")
#'
#' mod_apca <- apca(assessment ~ candy + assessor, data = candies_train)
#' pred_apca <- predict(mod_apca, newdata = candies_test)
#'
#' mod_msca <- msca(assessment ~ candy, data = candies_train)
#' pred_msca <- predict(mod_msca, newdata = candies_test)
#'
#' mod_apls <- apls(assessment ~ candy + assessor, data = candies_train)
#' pred_apls <- predict(mod_apls, newdata = candies_test)
#'
#' mod_limmpca <- limmpca(assessment ~ candy + r(assessor),
#'                        data = candies_train, pca.in = 3)
#' pred_limmpca <- predict(mod_limmpca, newdata = candies_test)
#'
#' @name predict_asca_family
NULL

.hda_predict_base <- function(object, newdata, ...){
  pred <- predict.hdanova(object, newdata = newdata, ...)
  pred$permute <- NULL
  pred
}

.hda_loadings_for_ls <- function(loadings_out, train_object, ncol_ls){
  if(nrow(loadings_out) == ncol_ls)
    return(loadings_out)

  pca_in <- train_object$more$pca.in
  if(!is.null(pca_in) && !identical(pca_in, 0) &&
     !is.null(train_object$Ypca) && !is.null(train_object$Ypca$pca) &&
     !is.null(train_object$Ypca$pca$loadings)){
    pca_loads <- train_object$Ypca$pca$loadings[, seq_len(pca_in), drop = FALSE]
    if(nrow(loadings_out) == nrow(pca_loads) && ncol_ls == ncol(pca_loads))
      return(t(pca_loads) %*% loadings_out)
  }

  stop("Training loadings are incompatible with predicted LS dimensions for projection.")
}

.hda_project_sca <- function(train_object, pred_object){
  scores <- loadings <- projected <- singulars <- list()
  for(eff in names(train_object$scores)){
    if(!eff %in% names(train_object$loadings))
      next
    L_out <- train_object$loadings[[eff]]
    if(identical(eff, "Residuals")){
      L_work <- .hda_loadings_for_ls(L_out, train_object, ncol(pred_object$residuals))
      scores[[eff]] <- pred_object$residuals %*% L_work
      projected[[eff]] <- pred_object$residuals %*% L_work
    } else {
      if(!eff %in% names(pred_object$LS))
        next
      L_work <- .hda_loadings_for_ls(L_out, train_object, ncol(pred_object$LS[[eff]]))
      scores[[eff]] <- pred_object$LS[[eff]] %*% L_work
      if(!is.null(pred_object$error[[eff]]))
        projected[[eff]] <- pred_object$error[[eff]] %*% L_work
      else
        projected[[eff]] <- pred_object$residuals %*% L_work
    }
    attr(scores[[eff]], "explvar") <- attr(train_object$scores[[eff]], "explvar")
    loadings[[eff]] <- L_out
    if(!is.null(train_object$singulars[[eff]]))
      singulars[[eff]] <- train_object$singulars[[eff]]
  }
  pred_object$scores <- scores
  pred_object$loadings <- loadings
  pred_object$projected <- projected
  if(length(singulars) > 0)
    pred_object$singulars <- singulars
  class(pred_object) <- c("asca", setdiff(class(pred_object), "asca"))
  pred_object
}

.hda_project_pls <- function(train_object, pred_object){
  scores <- loadings <- projected <- list()
  for(eff in names(train_object$scores)){
    if(!eff %in% names(train_object$loadings))
      next
    L_out <- train_object$loadings[[eff]]
    if(eff == "Residuals"){
      L_work <- .hda_loadings_for_ls(L_out, train_object, ncol(pred_object$residuals))
      scores[[eff]] <- pred_object$residuals %*% L_work
      projected[[eff]] <- pred_object$residuals %*% L_work
    } else {
      if(!eff %in% names(pred_object$LS))
        next
      L_work <- .hda_loadings_for_ls(L_out, train_object, ncol(pred_object$LS[[eff]]))
      scores[[eff]] <- pred_object$LS[[eff]] %*% L_work
      projected[[eff]] <- pred_object$residuals %*% L_work + scores[[eff]]
    }
    attr(scores[[eff]], "explvar") <- attr(train_object$scores[[eff]], "explvar")
    loadings[[eff]] <- L_out
  }
  pred_object$scores <- scores
  pred_object$loadings <- loadings
  pred_object$projected <- projected
  class(pred_object) <- c("apls", "asca", setdiff(class(pred_object), c("apls", "asca")))
  pred_object
}

.hda_rebuild_msca_within <- function(object){
  between <- object$model.frame[[2]]
  object$scores.within <- list()
  object$explvar.within <- matrix(0, nlevels(between), ncol(object$scores[[2]]))
  for(i in seq_len(nlevels(between))){
    keep_i <- between == levels(between)[i]
    object$scores.within[[i]] <- object$scores[[2]][keep_i, ]
    for(j in seq_len(ncol(object$scores[[2]]))){
      res_i <- object$residuals[keep_i, , drop = FALSE]
      fit_i <- tcrossprod(object$scores[[2]][keep_i, j], object$loadings[[2]][, j])
      object$explvar.within[i, j] <- 100 * (1 - sum((res_i - fit_i)^2) / sum(res_i^2))
    }
  }
  names(object$scores.within) <- levels(between)
  dimnames(object$explvar.within) <- list(levels(between), colnames(object$scores[[2]]))
  for(i in seq_len(nlevels(between)))
    attr(object$scores.within[[i]], "explvar") <- object$explvar.within[i, ]
  object
}

#' @rdname predict_asca_family
#' @method predict asca
#' @export
predict.asca <- function(object, newdata, decomposition = c("project", "refit"), ...){
  decomposition <- match.arg(decomposition)
  pred <- .hda_predict_base(object, newdata = newdata, ...)
  if(identical(decomposition, "project"))
    pred <- .hda_project_sca(object, pred)
  else
    pred <- sca(pred)
  pred$call <- match.call()
  pred
}

#' @rdname predict_asca_family
#' @method predict apca
#' @export
predict.apca <- function(object, newdata, decomposition = c("project", "refit"), ...){
  decomposition <- match.arg(decomposition)
  pred <- .hda_predict_base(object, newdata = newdata, ...)
  if(identical(decomposition, "project"))
    pred <- .hda_project_sca(object, pred)
  else
    pred <- sca(pred)
  pred$call <- match.call()
  class(pred) <- c("apca", setdiff(class(pred), "apca"))
  pred
}

#' @rdname predict_asca_family
#' @method predict msca
#' @export
predict.msca <- function(object, newdata, decomposition = c("project", "refit"), ...){
  decomposition <- match.arg(decomposition)
  pred <- .hda_predict_base(object, newdata = newdata, ...)
  if(identical(decomposition, "project"))
    pred <- .hda_project_sca(object, pred)
  else
    pred <- sca(pred)
  pred <- .hda_rebuild_msca_within(pred)
  pred$call <- match.call()
  class(pred) <- c("msca", setdiff(class(pred), "msca"))
  pred
}

#' @rdname predict_asca_family
#' @method predict apls
#' @export
predict.apls <- function(object, newdata, decomposition = c("project", "refit"), ...){
  decomposition <- match.arg(decomposition)
  pred <- .hda_predict_base(object, newdata = newdata, ...)
  if(identical(decomposition, "project"))
    pred <- .hda_project_pls(object, pred)
  else
    pred <- pls(pred)
  pred$call <- match.call()
  class(pred) <- c("apls", setdiff(class(pred), "apls"))
  pred
}

#' @rdname predict_asca_family
#' @method predict limmpca
#' @export
predict.limmpca <- function(object, newdata, decomposition = c("project", "refit"), ...){
  decomposition <- match.arg(decomposition)
  pred <- .hda_predict_base(object, newdata = newdata, ...)
  if(identical(decomposition, "project"))
    pred <- .hda_project_sca(object, pred)
  else
    pred <- sca(pred)
  pred$call <- match.call()
  class(pred) <- c("limmpca", setdiff(class(pred), "limmpca"))
  pred
}
