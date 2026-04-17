library(car)
library(mixlm)
library(pracma)
library(lme4)

source("R/ML_variance_partition.R")
source("R/hdanova_denominator.R")
source("R/hdanova.R")
source("R/utilities.R")
source("R/asca_results.R")
source("R/sca.R")
source("R/asca.R")

set.seed(20260416)

simulate_small_random_case <- function(n_subject = 18, n_rep = 2, n_var = 6,
                                       sigma_subject = 0.18, sigma_error = 0.35){
  design <- expand.grid(
    Subject = factor(seq_len(n_subject)),
    Trt1 = factor(c("A", "B")),
    Trt2 = factor(c("C", "D")),
    Rep = seq_len(n_rep)
  )
  design <- design[order(design$Subject, design$Trt1, design$Trt2, design$Rep), ]

  trt1_num <- ifelse(design$Trt1 == "B", 1, -1)
  trt2_num <- ifelse(design$Trt2 == "D", 1, -1)
  int_num <- trt1_num * trt2_num

  beta_trt1 <- seq(1.0, 1.5, length.out = n_var)
  beta_trt2 <- seq(0.8, 1.3, length.out = n_var)
  beta_int <- seq(0.5, 0.9, length.out = n_var)

  subject_effect <- matrix(rnorm(n_subject * n_var, sd = sigma_subject),
                           nrow = n_subject, ncol = n_var)
  error_term <- matrix(rnorm(nrow(design) * n_var, sd = sigma_error),
                       nrow = nrow(design), ncol = n_var)

  Y <- matrix(0, nrow = nrow(design), ncol = n_var)
  for(j in seq_len(n_var)){
    Y[, j] <- beta_trt1[j] * trt1_num +
      beta_trt2[j] * trt2_num +
      beta_int[j] * int_num +
      subject_effect[as.integer(design$Subject), j] +
      error_term[, j]
  }

  colnames(Y) <- paste0("Var", seq_len(n_var))
  design$Y <- I(Y)
  design
}

term_ssq <- function(object, terms){
  ssq <- object$ssq[terms]
  unname(ssq)
}

relative_diff <- function(x, y){
  abs(x - y) / pmax(abs(y), 1e-12)
}

compare_models <- function(data, method){
  mod_mixed <- asca(
    Y ~ Trt1 * Trt2 + r(Subject),
    data = data,
    REML = TRUE,
    SStype = "III",
    REML_ssq_method = method
  )
  mod_fixed <- asca(
    Y ~ Trt1 * Trt2,
    data = data,
    SStype = "III",
    respect_SStype = TRUE
  )

  total_y_ssq <- sum(scale(data$Y, scale = FALSE)^2)
  common_terms <- c("Trt1", "Trt2", "Trt1:Trt2")
  mixed_fixed_ssq <- term_ssq(mod_mixed, common_terms)
  fixed_ssq <- term_ssq(mod_fixed, common_terms)

  list(
    method = method,
    mixed = mod_mixed,
    fixed = mod_fixed,
    total_y_ssq = total_y_ssq,
    mixed_total_ratio = sum(mod_mixed$ssq) / total_y_ssq,
    fixed_total_ratio = sum(mod_fixed$ssq) / total_y_ssq,
    mixed_fixed_ssq = mixed_fixed_ssq,
    fixed_ssq = fixed_ssq,
    random_fraction = unname(mod_mixed$ssq["Subject"]) / total_y_ssq,
    fixed_relative_diff = setNames(relative_diff(mixed_fixed_ssq, fixed_ssq), common_terms)
  )
}

report_case <- function(result){
  cat("\n== Method:", result$method, "==\n")
  cat("Total Y SSQ:", round(result$total_y_ssq, 4), "\n")
  cat("Mixed-model SSQ / total Y SSQ:", round(result$mixed_total_ratio, 4), "\n")
  cat("Fixed-model SSQ / total Y SSQ:", round(result$fixed_total_ratio, 4), "\n")
  cat("Random-effect fraction of total Y SSQ:", round(result$random_fraction, 4), "\n")
  print(round(cbind(mixed = result$mixed_fixed_ssq,
                    fixed = result$fixed_ssq,
                    rel_diff = result$fixed_relative_diff), 4))
}

dat <- simulate_small_random_case()
methods <- c("exact_refit", "wald", "ls")
results <- lapply(methods, function(method) compare_models(dat, method))
names(results) <- methods

invisible(lapply(results, report_case))

for(method in methods){
  res <- results[[method]]
  stopifnot(is.finite(res$mixed_total_ratio), res$mixed_total_ratio > 0.90, res$mixed_total_ratio < 1.05)
  stopifnot(is.finite(res$fixed_total_ratio), res$fixed_total_ratio > 0.90, res$fixed_total_ratio < 1.05)
  stopifnot(is.finite(res$random_fraction), res$random_fraction >= 0, res$random_fraction < 0.05)
  stopifnot(all(is.finite(res$fixed_relative_diff)))
  stopifnot(all(res$fixed_relative_diff < 0.15))
}

cat("\nAll SSQ sanity checks passed.\n")