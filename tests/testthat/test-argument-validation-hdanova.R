test_that("hdanova validates core scalar arguments", {
  data(candies)

  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, unrestricted = c(TRUE, FALSE)),
    "'unrestricted' must be TRUE or FALSE"
  )
  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, respect_SStype = NA),
    "'respect_SStype' must be TRUE or FALSE"
  )
  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, use_ED = "yes"),
    "'use_ED' must be TRUE or FALSE"
  )
  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, REML = 1),
    "'REML' must be NULL, TRUE or FALSE"
  )
})

test_that("hdanova validates aug_error domain and type", {
  data(candies)

  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, aug_error = "foo"),
    "'aug_error' must be"
  )
  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, aug_error = -0.1),
    "Numeric 'aug_error' must be in \\[0,1\\]"
  )
  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, aug_error = 1.1),
    "Numeric 'aug_error' must be in \\[0,1\\]"
  )
})

test_that("hdanova rejects unsupported formula modes and incompatible compression", {
  data(candies)

  expect_error(
    hdanova(assessment ~ candy + (1 | assessor), data = candies),
    "supports mixlm-style random effects r\\(\\), not lme4 '\\\\|' notation"
  )
  expect_error(
    hdanova(assessment ~ candy + assessor, data = candies, pca.in = 2, pls.in = 2),
    "Only one of 'pca.in' and 'pls.in' can be active at a time"
  )
})
