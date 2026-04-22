test_that("REML without r() warns and falls back to fixed-effects", {
  data(candies)

  expect_warning(
    mod <- hdanova(assessment ~ candy + assessor, data = candies, REML = TRUE),
    "no random-effects r\\(\\) terms were found"
  )

  expect_true(is.null(mod$more$REML))
  expect_true(inherits(mod, "hdanova"))
})
