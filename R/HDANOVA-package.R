#' @aliases HDANOVA
#' @title High-Dimensional Analysis of Variance
#'
#' @description Included methods:
#' * ASCA+ - Analysis of Variance Simultaneous Component Analysis
#' * APCA+ - ANOVA Principal Component Analysis
#' * LiMM-PCA - Linear Mixed Model PCA
#' * MSCA - Multilevel Simultaneous Component Analysis
#' * PC-ANOVA - Principal Component Analysis of Variance
#' * PRC - Principal Response Curves
#' * PERMANOVA - Permutation Based MANOVA
#'
#' @importFrom mixlm lm
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{hdanova}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#' @docType package
#' @name HDANOVA
#' @keywords internal
"_PACKAGE"
