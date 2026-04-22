test_that("use_ED warns and is inactive for non-numeric aug_error", {
  data(candies)

  expect_warning(
    mod <- hdanova(
      assessment ~ candy + r(assessor),
      data = candies,
      aug_error = "denominator",
      use_ED = TRUE,
      REML = TRUE
    ),
    "currently only affects numeric 'aug_error'"
  )

  expect_null(mod$ED)
})

test_that("use_ED warns and is inactive for fixed-effects fits", {
  data(candies)

  expect_warning(
    mod <- hdanova(
      assessment ~ candy + assessor,
      data = candies,
      aug_error = 0.05,
      use_ED = TRUE
    ),
    "currently only affects numeric 'aug_error'"
  )

  expect_null(mod$ED)
})

test_that("use_ED is active for numeric aug_error in REML mixed fits", {
  data(candies)

  mod <- hdanova(
    assessment ~ candy + r(assessor),
    data = candies,
    aug_error = 0.05,
    use_ED = TRUE,
    REML = TRUE
  )

  expect_true(is.matrix(mod$ED))
  expect_gt(nrow(mod$ED), 0)
  expect_gt(ncol(mod$ED), 0)
  expect_true(all(is.finite(mod$ED)))
})

test_that("permutation warns for REML + respect_SStype and records SStype metadata", {
  data(candies)

  mod <- asca(
    assessment ~ candy + r(assessor),
    data = candies,
    SStype = "III",
    REML = TRUE,
    REML_ssq_method = "wald"
  )

  expect_warning(
    p <- permutation(mod, permute = 10, perm.type = "approximate", respect_SStype = TRUE),
    "uses SS-type-aligned QR permutation statistics"
  )

  expect_identical(p$permute$effect_source, "SStype")
  expect_identical(p$permute$ssq_method, "wald")
  expect_true(all(is.finite(p$permute$pvalues) | is.na(p$permute$pvalues)))
})

test_that("permutation default path keeps regression effect source", {
  data(candies)

  mod <- asca(
    assessment ~ candy + r(assessor),
    data = candies,
    SStype = "III",
    REML = TRUE,
    REML_ssq_method = "wald"
  )

  p <- permutation(mod, permute = 10, perm.type = "approximate", respect_SStype = FALSE)

  expect_identical(p$permute$effect_source, "regression")
  expect_identical(p$permute$ssq_method, "wald")
  expect_true(all(is.finite(p$permute$pvalues) | is.na(p$permute$pvalues)))
})
