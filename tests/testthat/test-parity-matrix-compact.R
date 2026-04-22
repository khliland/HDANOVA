test_that("compact parity matrix: fixed and REML paths are numerically sane", {
  data(candies)

  candies_unbal <- candies[-seq(1, nrow(candies), by = 4), , drop = FALSE]

  cases <- list(
    list(
      label = "fixed-II-approx",
      formula = assessment ~ candy * assessor,
      data = candies_unbal,
      fit_args = list(SStype = "II", respect_SStype = TRUE),
      perm_args = list(permute = 8, perm.type = "approximate", respect_SStype = TRUE)
    ),
    list(
      label = "fixed-III-exact",
      formula = assessment ~ candy * assessor,
      data = candies_unbal,
      fit_args = list(SStype = "III", respect_SStype = TRUE),
      perm_args = list(permute = 8, perm.type = "exact", respect_SStype = TRUE)
    ),
    list(
      label = "reml-II-approx",
      formula = assessment ~ candy + r(assessor),
      data = candies,
      fit_args = list(SStype = "II", REML = TRUE, REML_ssq_method = "wald", respect_SStype = TRUE),
      perm_args = list(permute = 8, perm.type = "approximate", respect_SStype = TRUE)
    ),
    list(
      label = "reml-III-exact",
      formula = assessment ~ candy + r(assessor),
      data = candies,
      fit_args = list(SStype = "III", REML = TRUE, REML_ssq_method = "wald", respect_SStype = TRUE),
      perm_args = list(permute = 8, perm.type = "exact", respect_SStype = TRUE)
    )
  )

  for (case in cases) {
    fit_call <- c(list(formula = case$formula, data = case$data), case$fit_args)
    mod <- do.call(hdanova, fit_call)

    fixed_names <- if (!is.null(mod$more$effs_fixed)) mod$more$effs_fixed else mod$more$effs
    fixed_names <- setdiff(fixed_names, "Residuals")

    expect_true(all(is.finite(unname(mod$ssq[intersect(names(mod$ssq), fixed_names)]))),
                info = case$label)

    perm_call <- c(list(object = mod), case$perm_args,
                   list(unique.digits = NULL, unique.frac = NULL, exhaustive.warn = FALSE))
    p <- suppressWarnings(do.call(permutation, perm_call))

    p_fixed <- p$permute$pvalues[intersect(names(p$permute$pvalues), fixed_names)]
    expect_true(all(is.finite(p_fixed) | is.na(p_fixed)), info = case$label)

    ssqa_fixed <- p$permute$ssqa[intersect(names(p$permute$ssqa), fixed_names)]
    expect_true(all(is.finite(ssqa_fixed) | is.na(ssqa_fixed)), info = case$label)
  }
})
