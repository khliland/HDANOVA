test_that("permutation validates control arguments", {
  data(candies)
  mod <- asca(assessment ~ candy + assessor, data = candies)

  expect_error(
    permutation(mod, unique.digits = "12"),
    "'unique.digits' must be a single numeric value or NULL"
  )
  expect_error(
    permutation(mod, unique.frac = 1.5),
    "'unique.frac' must be a single numeric value in \\[0, 1\\] or NULL"
  )
  expect_error(
    permutation(mod, exhaustive.warn = NA),
    "'exhaustive.warn' must be TRUE or FALSE"
  )
})

test_that("permutation respects object default when respect_SStype is NULL", {
  data(candies)
  mod <- asca(assessment ~ candy + assessor, data = candies, respect_SStype = TRUE)

  p <- suppressWarnings(permutation(mod, permute = 10, perm.type = "approximate", respect_SStype = NULL))
  expect_identical(p$permute$effect_source, "SStype")
})

test_that("exact permutation emits exhaustive warning on tiny design", {
  tiny <- data.frame(
    g = factor(c("A", "A", "B", "B")),
    y = I(matrix(c(1, 2, 3, 4), ncol = 1))
  )

  mod <- asca(y ~ g, data = tiny)

  expect_warning(
    permutation(mod, permute = 1000, perm.type = "exact", respect_SStype = TRUE,
                unique.digits = NULL, unique.frac = NULL),
    "uses exhaustive permutation"
  )
})
