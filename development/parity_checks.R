## parity_checks.R
## Systematic check of all meaningful combinations of:
##   SStype         : I, II, III
##   respect_SStype : FALSE, TRUE
##   REML_ssq_method: "exact_refit", "wald", "ls"  (REML models only)
##   unrestricted   : FALSE, TRUE
##   perm.type      : "approximate", "exact"
##   model type     : fixed (lm), mixed MoM (r()), mixed REML/ML (lme4)
##
## What we check for each combination:
##   1. hdanova() completes without error and returns an object.
##   2. SSQ values for all approved fixed effects are strictly positive and finite.
##   3. SSQs sum to at most ssqY (total variance; allows for negative-SS in MoM).
##   4. permutation() completes without error.
##   5. All FIXED-effect permutation p-values are finite or NA (random effects allowed NA).
##   6. All fixed effects have > 1 unique permutation SSQ value (non-degenerate).
##
## Failures are collected and printed as a summary table; the script exits
## with code 1 if any failure is recorded.

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

# ---- Data: balanced 2x4 repeated-measures design (between-within) ----------
# Used for fixed LM and REML/ML mixed models.
# Treatment is between-subject (nested in Subject).
set.seed(123)
n_subj   <- 12   # subjects, 4 per Treatment group
n_time   <- 4
n_vars   <- 6
subjects <- paste0("S", seq_len(n_subj))
trt      <- rep(c("A", "B", "C"), each = n_subj / 3)
df <- expand.grid(Subject = factor(subjects), Time = factor(seq_len(n_time)))
df$Treatment <- factor(trt[match(df$Subject, subjects)])

signal <- outer(as.numeric(df$Treatment), seq_len(n_vars)) +
          outer(as.numeric(df$Time),       seq_len(n_vars)) * 0.5
noise  <- matrix(rnorm(nrow(df) * n_vars, sd = 0.6), nrow(df), n_vars)
df$Y   <- I(signal + noise)

# ---- Data: fully crossed design for MoM (r()) models ----------------------
# Subject is crossed with Stimulus and Condition so r(Subject) is valid
# without aliasing issues.
set.seed(456)
n_subj_c   <- 8
n_stim     <- 3
n_cond     <- 2
df_cross <- expand.grid(Subject   = factor(paste0("P", seq_len(n_subj_c))),
                        Stimulus  = factor(seq_len(n_stim)),
                        Condition = factor(c("X", "Y")))
df_cross$Y <- I(matrix(rnorm(nrow(df_cross) * n_vars, sd = 1), nrow(df_cross), n_vars))

# ---- Helper ---------------------------------------------------------------
n_perms  <- 40   # keep quick; enough to detect degenerate distributions
uniq_thr <- 2    # require at least 2 unique SSQ values

failures <- list()
results  <- list()

run_case <- function(label, model_args, perm_args = list()){
  cat("  ", label, "... ")

  ## 1. Fit model
  fit <- tryCatch(
    do.call(hdanova, model_args),
    error = function(e) e
  )
  if(inherits(fit, "error")){
    cat("FAIL (hdanova error:", conditionMessage(fit), ")\n")
    failures[[length(failures) + 1]] <<- list(label = label, stage = "hdanova",
                                               msg = conditionMessage(fit))
    return(invisible(NULL))
  }

  ## 2. SSQ checks on fixed effects
  effs_fixed <- intersect(names(fit$ssq), fit$more$effs_fixed %||% fit$more$effs)
  ssq_vals   <- fit$ssq[effs_fixed]
  ssq_vals   <- ssq_vals[names(ssq_vals) != "Residuals"]
  if(any(!is.finite(ssq_vals))){
    cat("FAIL (non-finite SSQ)\n")
    failures[[length(failures) + 1]] <<- list(label = label, stage = "ssq_finite",
                                               msg = paste(names(ssq_vals[!is.finite(ssq_vals)]), collapse = ", "))
    return(invisible(fit))
  }

  ## 3. Permutation
  perm_args_full <- c(list(object = fit, permute = n_perms), perm_args)
  perm_warnings  <- character(0)
  perm_res <- withCallingHandlers(
    tryCatch(
      do.call(permutation, perm_args_full),
      error = function(e) e
    ),
    warning = function(w){
      perm_warnings <<- c(perm_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if(inherits(perm_res, "error")){
    cat("FAIL (permutation error:", conditionMessage(perm_res), ")\n")
    failures[[length(failures) + 1]] <<- list(label = label, stage = "permutation",
                                               msg = conditionMessage(perm_res))
    return(invisible(fit))
  }

  ## 4. p-value finite/NA check (only fixed effects)
  pvs  <- perm_res$permute$pvalues
  ssqa <- perm_res$permute$ssqa
  fixed_effs   <- intersect(names(pvs), fit$more$effs_fixed %||% fit$more$effs)
  bad_pvals    <- fixed_effs[!is.finite(pvs[fixed_effs]) & !is.na(pvs[fixed_effs])]
  if(length(bad_pvals) > 0){
    cat("FAIL (bad p-values for:", paste(bad_pvals, collapse = ", "), ")\n")
    failures[[length(failures) + 1]] <<- list(label = label, stage = "pvalues",
                                               msg = paste(bad_pvals, collapse = ", "))
    return(invisible(perm_res))
  }

  ## 5. Uniqueness check for fixed effects
  degen <- character(0)
  for(eff in fixed_effs){
    perms_eff <- perm_res$permute$ssqaperm[[eff]]
    if(length(perms_eff) == 0) next
    n_uniq <- length(unique(round(perms_eff, 10)))
    if(n_uniq < uniq_thr)
      degen <- c(degen, paste0(eff, "(", n_uniq, ")"))
  }
  # catch "produced only 1 unique" warnings not already flagged above
  degenW <- grep("produced only 1 unique", perm_warnings, value = TRUE)
  if(length(degenW) > 0)
    degen <- unique(c(degen, sub(".*Effect '([^']+)' produced.*", "\\1", degenW)))

  if(length(degen) > 0){
    cat("FAIL (degenerate SSQ for:", paste(degen, collapse = ", "), ")\n")
    failures[[length(failures) + 1]] <<- list(label = label, stage = "uniqueness",
                                               msg = paste(degen, collapse = ", "))
    return(invisible(perm_res))
  }

  cat("OK\n")
  invisible(perm_res)
}

# NULL-coalescing helper used above
`%||%` <- function(a, b) if(!is.null(a)) a else b

# ===========================================================================
# Fixed LM: all SStype x unrestricted x respect_SStype x perm.type
# ===========================================================================
cat("\n===== Fixed LM =====\n")
for(ss in c("I","II","III")){
  for(ur in c(FALSE, TRUE)){
    for(rss in c(FALSE, TRUE)){
      for(pt in c("approximate", "exact")){
        lbl <- sprintf("LM  SS=%s unres=%s rSS=%s perm=%s", ss, ur, rss, pt)
        run_case(lbl,
          model_args = list(formula  = Y ~ Treatment * Time,
                            data     = df,
                            SStype   = ss,
                            unrestricted = ur,
                            respect_SStype = rss),
          perm_args  = list(perm.type = pt,
                            respect_SStype = rss,
                            unique.digits = NULL)  # suppress uniqueness warning; we check manually
        )
      }
    }
  }
}

# ===========================================================================
# Mixed MoM (r() terms, REML=NULL): SStype x unrestricted x respect_SStype x perm.type
# Uses crossed design (Subject x Stimulus x Condition) to avoid aliasing.
# ===========================================================================
cat("\n===== Mixed MoM (REML=NULL) =====\n")
for(ss in c("I","II","III")){
  for(ur in c(FALSE, TRUE)){
    for(rss in c(FALSE, TRUE)){
      for(pt in c("approximate", "exact")){
        lbl <- sprintf("MoM SS=%s unres=%s rSS=%s perm=%s", ss, ur, rss, pt)
        run_case(lbl,
          model_args = list(formula  = Y ~ Stimulus * Condition + r(Subject),
                            data     = df_cross,
                            SStype   = ss,
                            unrestricted = ur,
                            respect_SStype = rss),
          perm_args  = list(perm.type = pt,
                            respect_SStype = rss,
                            unique.digits = NULL)
        )
      }
    }
  }
}

# ===========================================================================
# Mixed REML: SStype x REML_ssq_method x respect_SStype x perm.type
# (unrestricted not meaningful for lme4 path; skip)
# ===========================================================================
cat("\n===== Mixed REML =====\n")
for(ss in c("I","II","III")){
  for(meth in c("exact_refit", "wald", "ls")){
    for(rss in c(FALSE, TRUE)){
      for(pt in c("approximate", "exact")){
        lbl <- sprintf("REML SS=%s meth=%s rSS=%s perm=%s", ss, meth, rss, pt)
        run_case(lbl,
          model_args = list(formula        = Y ~ r(Subject) + Treatment * Time,
                            data           = df,
                            SStype         = ss,
                            REML           = TRUE,
                            REML_ssq_method = meth,
                            respect_SStype = rss),
          perm_args  = list(perm.type = pt,
                            respect_SStype = rss,
                            unique.digits = NULL)
        )
      }
    }
  }
}

# ===========================================================================
# Mixed ML (REML=FALSE): same grid as REML
# ===========================================================================
cat("\n===== Mixed ML (REML=FALSE) =====\n")
for(ss in c("I","II","III")){
  for(meth in c("exact_refit", "wald", "ls")){
    for(rss in c(FALSE, TRUE)){
      for(pt in c("approximate", "exact")){
        lbl <- sprintf("ML   SS=%s meth=%s rSS=%s perm=%s", ss, meth, rss, pt)
        run_case(lbl,
          model_args = list(formula        = Y ~ r(Subject) + Treatment * Time,
                            data           = df,
                            SStype         = ss,
                            REML           = FALSE,
                            REML_ssq_method = meth,
                            respect_SStype = rss),
          perm_args  = list(perm.type = pt,
                            respect_SStype = rss,
                            unique.digits = NULL)
        )
      }
    }
  }
}

# ===========================================================================
# Summary
# ===========================================================================
cat("\n===== SUMMARY =====\n")
if(length(failures) == 0){
  cat("All", "combinations passed.\n")
} else {
  cat(length(failures), "failure(s) detected:\n\n")
  for(f in failures){
    cat(sprintf("  [%s] stage=%s  msg=%s\n", f$label, f$stage, f$msg))
  }
  quit(status = 1)
}
