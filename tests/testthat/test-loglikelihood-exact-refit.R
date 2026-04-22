test_that("exact_refit mixed-model SSQ stores expected public metadata", {
  data(candies)

  mod <- hdanova(
    assessment ~ candy + r(assessor),
    data = candies,
    REML = TRUE,
    REML_ssq_method = "exact_refit"
  )

  expect_identical(mod$more$ssq_method, "exact_refit")
  expect_identical(names(mod$ssq), c("candy", "assessor", "Residuals"))
  expect_identical(names(mod$dfNum), c("candy", "assessor", "Residuals"))
  expect_identical(names(mod$dfDenom), c("candy", "assessor", "Residuals"))
  expect_identical(unname(mod$denoms), c(3, 3, NA_real_))
  expect_true(all(is.finite(unname(mod$ssq))))
  expect_true(all(unname(mod$ssq) > 0))
  expect_true(all(is.finite(unname(mod$dfNum[1:2]))))
  expect_true(all(is.finite(unname(mod$dfDenom[1:2]))))
})

test_that("exact_refit respects SS type and REML/ML mode in public SSQ output", {
  data(candies)

  mod_reml_i <- hdanova(
    assessment ~ candy + r(assessor),
    data = candies,
    REML = TRUE,
    SStype = 1,
    REML_ssq_method = "exact_refit"
  )
  mod_reml_ii <- hdanova(
    assessment ~ candy + r(assessor),
    data = candies,
    REML = TRUE,
    SStype = 2,
    REML_ssq_method = "exact_refit"
  )
  mod_reml_iii <- hdanova(
    assessment ~ candy + r(assessor),
    data = candies,
    REML = TRUE,
    SStype = 3,
    REML_ssq_method = "exact_refit"
  )
  mod_ml_ii <- hdanova(
    assessment ~ candy + r(assessor),
    data = candies,
    REML = FALSE,
    SStype = 2,
    REML_ssq_method = "exact_refit"
  )

  expect_equal(unname(mod_reml_i$ssq[c("candy", "assessor")]),
               unname(mod_reml_ii$ssq[c("candy", "assessor")]),
               tolerance = 1e-8)
  expect_equal(unname(mod_reml_i$ssq["Residuals"]),
               unname(mod_reml_ii$ssq["Residuals"]),
               tolerance = 1e-8)

  expect_lt(unname(mod_reml_iii$ssq["candy"]), unname(mod_reml_ii$ssq["candy"]))
  expect_equal(unname(mod_reml_iii$ssq[c("assessor", "Residuals")]),
               unname(mod_reml_ii$ssq[c("assessor", "Residuals")]),
               tolerance = 1e-8)

  expect_equal(unname(mod_ml_ii$ssq["candy"]),
               unname(mod_reml_ii$ssq["candy"]),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(
    unname(mod_ml_ii$ssq[c("assessor", "Residuals")]),
    unname(mod_reml_ii$ssq[c("assessor", "Residuals")]),
    tolerance = 1e-8
  )))
})