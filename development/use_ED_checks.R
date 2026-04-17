library(car)
library(mixlm)
library(pracma)
library(lme4)

source("R/ML_variance_partition.R")
source("R/hdanova_denominator.R")
legacy_env <- new.env(parent = .GlobalEnv)
source("development/hdanova_legacy.R", local = legacy_env)
legacy_hdanova <- legacy_env$hdanova
source("R/hdanova.R")

data(candies, package = "HDANOVA")

cat("== use_ED REML mixed parity ==\n")
legacy <- legacy_hdanova(
  assessment ~ candy + r(assessor),
  data = candies,
  aug_error = 0.05,
  use_ED = TRUE,
  REML = TRUE
)
ported <- hdanova(
  assessment ~ candy + r(assessor),
  data = candies,
  aug_error = 0.05,
  use_ED = TRUE,
  REML = TRUE
)

stopifnot(!is.null(legacy$ED), !is.null(ported$ED))
stopifnot(identical(dim(legacy$ED), dim(ported$ED)))
max_ed_diff <- max(abs(legacy$ED - ported$ED))
cat("max ED absolute difference:", max_ed_diff, "\n")
stopifnot(max_ed_diff < 1e-12)

cat("ported ED rows:", paste(rownames(ported$ED), collapse = ", "), "\n")
cat("ported fit type:", ported$fit.type, "\n")
cat("ported effect source:", ported$more$effect_source, "\n")

cat("\n== LS_aug numeric use_ED parity (REML mixed) ==\n")
common_num_eff <- intersect(names(legacy$more$LS_aug), names(ported$more$LS_aug))
for(eff in common_num_eff){
  max_diff <- max(abs(legacy$more$LS_aug[[eff]] - ported$more$LS_aug[[eff]]))
  cat(eff, " max LS_aug absolute difference:", max_diff, "\n")
  stopifnot(max_diff < 1e-12)
}

cat("\n== LS_aug denominator parity (REML mixed) ==\n")
legacy_den <- legacy_hdanova(
  assessment ~ candy + r(assessor),
  data = candies,
  aug_error = "denominator",
  REML = TRUE
)
ported_den <- hdanova(
  assessment ~ candy + r(assessor),
  data = candies,
  aug_error = "denominator",
  REML = TRUE
)
common_eff <- intersect(names(legacy_den$more$LS_aug), names(ported_den$more$LS_aug))
for(eff in common_eff){
  max_diff <- max(abs(legacy_den$more$LS_aug[[eff]] - ported_den$more$LS_aug[[eff]]))
  cat(eff, " max LS_aug absolute difference:", max_diff, "\n")
  stopifnot(max_diff < 1e-12)
}

cat("\n== use_ED scope guardrails ==\n")

# use_ED only affects numeric aug_error in REML/ML mixed fits for hdanova.
warn_den <- FALSE
obj_den <- withCallingHandlers(
  hdanova(
    assessment ~ candy + r(assessor),
    data = candies,
    aug_error = "denominator",
    use_ED = TRUE,
    REML = TRUE
  ),
  warning = function(w){
    if(grepl("currently only affects numeric 'aug_error'", conditionMessage(w), fixed = TRUE))
      warn_den <<- TRUE
    invokeRestart("muffleWarning")
  }
)
stopifnot(isTRUE(warn_den))
stopifnot(is.null(obj_den$ED))
cat("non-numeric aug_error guardrail warning observed and ED omitted: TRUE\n")

warn_fixed <- FALSE
obj_fixed <- withCallingHandlers(
  hdanova(
    assessment ~ candy + assessor,
    data = candies,
    aug_error = 0.05,
    use_ED = TRUE
  ),
  warning = function(w){
    if(grepl("currently only affects numeric 'aug_error'", conditionMessage(w), fixed = TRUE))
      warn_fixed <<- TRUE
    invokeRestart("muffleWarning")
  }
)
stopifnot(isTRUE(warn_fixed))
stopifnot(is.null(obj_fixed$ED))
cat("non-REML mixed guardrail warning observed and ED omitted: TRUE\n")

cat("\nAll use_ED checks passed.\n")
