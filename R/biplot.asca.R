#' Biplot for ASCA models
#'
#' @param x \code{asca} object.
#' @param factor Factor number or name.
#' @param comps \code{integer} vector of selected components.
#' @param xlim \code{numeric} vector of length 2 for x-axis limits of the loadings.
#' @param ylim \code{numeric} vector of length 2 for y-axis limits of the loadings.
#' @param col \code{vector} of colours for score axes and loading axes and points/texts.
#' @param expand \code{numeric} expansion for the scores, defaulting to 1.
#' @param labels optional. If \code{"names"}, row names are used as labels.
#' If \code{"numbers"}, row numbers are used as labels. (Can also be a vector of labels.)
#' @param legendpos \code{character} position of legend.
#' @param ... Additional arguments to \code{plot} and \code{scoreplot}.
#'
#' @returns No return, only a plot.
#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics par text
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic ASCA model with two factors and interaction
#' mod <- asca(assessment ~ candy * assessor, data=candies)
#'
#' # Biplot
#' biplot(mod)
#'
#' # Biplot with named loadings
#' biplot(mod, labels="names")
#'
#' @export
biplot.asca <- function(x, factor=1, comps=1:2, xlim = NULL, ylim = NULL,
                        col="darkgray", expand=1, labels, legendpos, ...){
  # If length comps is 1, stop and warn
  if(length(comps) != 2){
    stop("Exactly two components are needed for a biplot.")
  }
  loads <- x$loadings[[factor]][,comps]
  scors <- x$projected[[factor]][,comps]
  unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)),
                                  abs(max(x, na.rm = TRUE)))
  rangL1 <- unsigned.range(loads[, 1L])
  rangL2 <- unsigned.range(loads[, 2L])
  rangS1 <- unsigned.range(scors[, 1L])
  rangS2 <- unsigned.range(scors[, 2L])
  if (missing(xlim) && missing(ylim))
    xlim <- ylim <- rangL1 <- rangL2 <- range(rangL1, rangL2)
  else if (missing(xlim))
    xlim <- rangL1
  else if (missing(ylim))
    ylim <- rangL2
  ratio <- max(rangS1/rangL1, rangS2/rangL2)/expand
  if(missing(labels)){
    plot(loads, pch=20, axes = FALSE, col=col, cex=0.8, xlab="", ylab="",
         xlim=xlim, ylim=ylim, ...)
  } else {
    plot(loads, pch=20, axes = FALSE, col=col, cex=0.8, xlab="", ylab="",
         xlim=xlim, ylim=ylim, type="n", ...)
    ## Set up point/tick mark labels
    if (length(labels) == 1) {
      labels <- switch(match.arg(labels, c("names", "numbers")),
                       names = {
                         if (is.null(rnames <- rownames(loads))) {
                           stop("The loadings have no row names.")
                         } else {
                           rnames
                         }},
                       numbers = 1:nrow(loads)
      )
    }
    labels <- as.character(labels)
    type <- "n"
    text(loads[,1], loads[,2], labels = labels, col = col[1L], ...)
  }
  axis(3, col = col[1L], col.axis = col[1L])
  axis(4, col = col[1L], col.axis = col[1L])
  par(new = TRUE)
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  class(scors) <- "scores"
  if(missing(legendpos))
    scoreplot(x, factor=factor, comps=comps, xlim=xlim*ratio, ylim=ylim*ratio, ...)
  else
    scoreplot(x, factor=factor, comps=comps, xlim=xlim*ratio, ylim=ylim*ratio, legendpos=legendpos, ...)
}
