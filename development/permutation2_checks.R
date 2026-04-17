library(car)
library(mixlm)
library(pracma)
library(lme4)

source("R/ML_variance_partition.R")
source("R/hdanova_denominator.R")
source("R/hdanova.R")
source("R/utilities.R")
source("R/permutation.R")
source("R/asca.R")
source("R/sca.R")

data(candies, package = "HDANOVA")

cat("== Fixed aligned path ==\n")
cand_unbal <- candies[-seq(1, nrow(candies), by = 4), , drop = FALSE]
m_legacy <- hdanova(assessment ~ candy * assessor, data = cand_unbal, SStype = "II")
m_align2 <- hdanova(assessment ~ candy * assessor, data = cand_unbal, SStype = "II",
                     respect_SStype = TRUE)
m_align3 <- hdanova(assessment ~ candy * assessor, data = cand_unbal, SStype = "III",
                     respect_SStype = TRUE)
cat("legacy keeps regression LS:",
    max(abs(m_legacy$LS[["candy"]] - m_legacy$LS_regression[["candy"]])), "\n")
cat("aligned uses SStype LS:",
    max(abs(m_align2$LS[["candy"]] - m_align2$LS_SStype[["candy"]])), "\n")
cat("aligned SStype sensitivity candy II-III:",
    max(abs(m_align2$LS[["candy"]] - m_align3$LS[["candy"]])), "\n")

p2_fixed <- permutation(m_align2, permute = 25, perm.type = "approximate")
cat("fixed permutation finite pvalues:", all(is.finite(p2_fixed$permute$pvalues)), "\n")

cat("\n== Mixed MoM aligned path ==\n")
m_mix <- hdanova(assessment ~ candy * r(assessor), data = candies,
                  SStype = "II", unrestricted = FALSE, respect_SStype = TRUE)
p2_mix_a <- permutation(m_mix, permute = 25, perm.type = "approximate")
p2_mix_e <- permutation(m_mix, permute = 25, perm.type = "exact")
cat("mixed approximate finite pvalues:", all(is.finite(p2_mix_a$permute$pvalues)), "\n")
cat("mixed exact effect names:", paste(names(p2_mix_e$permute$ssqa), collapse = ", "), "\n")
cat("mixed exact candy exchangeable blocks:",
    p2_mix_e$permute$exchangeable$candy$exchangeable_blocks, "\n")

cat("\n== Mixed REML permutation metadata ==\n")
m_reml <- asca(assessment ~ candy + r(assessor), data = candies,
                             SStype = "III", REML = TRUE, REML_ssq_method = "wald")
warned_reml_perm <- FALSE
p2_reml <- withCallingHandlers(
    permutation(m_reml, permute = 25, perm.type = "approximate"),
    warning = function(w){
        if(grepl("Permutation in REML/ML currently uses regression-projection SSQ statistics",
                         conditionMessage(w), fixed = TRUE))
            warned_reml_perm <<- TRUE
        invokeRestart("muffleWarning")
    }
)
cat("reml permutation finite pvalues:", all(is.finite(p2_reml$permute$pvalues)), "\n")
cat("reml permutation effect source:", p2_reml$permute$effect_source, "\n")
cat("reml permutation ssq_method:", p2_reml$permute$ssq_method, "\n")
cat("reml permutation warning observed:", warned_reml_perm, "\n")

cat("\n== Mixed REML aligned permutation ==\n")
warned_reml_align <- FALSE
p2_reml_align_a <- withCallingHandlers(
    permutation(m_reml, permute = 25, perm.type = "approximate", respect_SStype = TRUE),
    warning = function(w){
        if(grepl("uses SS-type-aligned QR permutation statistics", conditionMessage(w), fixed = TRUE))
            warned_reml_align <<- TRUE
        invokeRestart("muffleWarning")
    }
)
p2_reml_align_e <- withCallingHandlers(
    permutation(m_reml, permute = 25, perm.type = "exact", respect_SStype = TRUE),
    warning = function(w){
        if(grepl("uses SS-type-aligned QR permutation statistics", conditionMessage(w), fixed = TRUE))
            warned_reml_align <<- TRUE
        if(grepl("produced only", conditionMessage(w), fixed = TRUE))
            invokeRestart("muffleWarning")
        else
            invokeRestart("muffleWarning")
    }
)
cat("reml aligned approximate finite pvalues:", all(is.finite(p2_reml_align_a$permute$pvalues) | is.na(p2_reml_align_a$permute$pvalues)), "\n")
cat("reml aligned exact finite/NA pvalues:", all(is.finite(p2_reml_align_e$permute$pvalues) | is.na(p2_reml_align_e$permute$pvalues)), "\n")
cat("reml aligned effect source (approx):", p2_reml_align_a$permute$effect_source, "\n")
cat("reml aligned effect source (exact):", p2_reml_align_e$permute$effect_source, "\n")
cat("reml aligned warning observed:", warned_reml_align, "\n")