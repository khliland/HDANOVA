library(car)
library(mixlm)
library(pracma)
library(lme4)
library(pls)

source("R/ML_variance_partition.R")
source("R/hdanova_denominator.R")
source("R/hdanova.R")
source("R/utilities.R")
source("R/asca_results.R")
source("R/sca.R")
source("R/pls.R")
source("R/signflip.R")
source("R/asca_plots.R")
source("R/biplot.asca.R")
source("R/permutation.R")
source("R/asca.R")
source("R/apca.R")
source("R/pcanova.R")
source("R/msca.R")
source("R/apls.R")
source("R/limmpca.R")

data(candies, package = "HDANOVA")

run_case <- function(name, expr){
  ok <- TRUE
  msg <- ""
  res <- tryCatch(expr, error = function(e){ ok <<- FALSE; msg <<- conditionMessage(e); NULL })
  if(ok){
    cat(name, ": OK\n", sep = "")
    if(is.list(res) && !is.null(res$fit.type)) cat("  fit.type=", res$fit.type, "\n", sep = "")
    if(is.list(res) && !is.null(res$more$effect_source)) cat("  effect_source=", res$more$effect_source, "\n", sep = "")
    if(is.list(res) && !is.null(res$permute$effect_source)) cat("  permute_effect_source=", res$permute$effect_source, "\n", sep = "")
  } else {
    cat(name, ": ERROR\n", sep = "")
    cat("  message=", msg, "\n", sep = "")
  }
  invisible(res)
}

cat("Rehearsal: canonical hdanova/permutation engine route\n")
obj_asca <- run_case("asca()", asca(assessment ~ candy + assessor, data = candies))
obj_asca_p <- run_case("asca(permute=20)", asca(assessment ~ candy + assessor, data = candies, permute = 20))
obj_apca <- run_case("apca()", apca(assessment ~ candy + assessor, data = candies))
obj_pcan <- run_case("pcanova()", pcanova(assessment ~ candy + assessor, data = candies))
obj_msca <- run_case("msca()", msca(assessment ~ candy + assessor, data = candies))
obj_apls <- run_case("apls()", apls(assessment ~ candy + assessor, data = candies))
obj_limm <- run_case("limmpca()", limmpca(assessment ~ candy + r(assessor), data = candies, pca.in = 2, aug_error = 0.05, use_ED = TRUE, REML = TRUE))

if(!is.null(obj_asca)){
  run_case("scoreplot(asca)", { pdf("/tmp/rehearsal_asca_scoreplot.pdf"); scoreplot(obj_asca, factor = "candy"); dev.off(); TRUE })
  run_case("biplot(asca)", { pdf("/tmp/rehearsal_asca_biplot.pdf"); biplot(obj_asca, factor = "candy"); dev.off(); TRUE })
}
if(!is.null(obj_msca))
  run_case("scoreplot(msca)", { pdf("/tmp/rehearsal_msca_scoreplot.pdf"); scoreplot(obj_msca, factor = "Residuals"); dev.off(); TRUE })
