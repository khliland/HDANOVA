test_that("hdanova accepts typical single contrast aliases", {
  data(candies)

  mod_ref <- hdanova(
    assessment ~ candy + assessor,
    data = candies,
    contrasts = "reference"
  )
  expect_equal(unname(mod_ref$contrasts), "contr.reference")

  mod_treat <- hdanova(
    assessment ~ candy + assessor,
    data = candies,
    contrasts = "treatment"
  )
  expect_equal(unname(mod_treat$contrasts), "contr.treatment")

  mod_sum <- hdanova(
    assessment ~ candy + assessor,
    data = candies,
    contrasts = "sum"
  )
  expect_equal(unname(mod_sum$contrasts), "contr.sum")

  mod_ref_full <- hdanova(
    assessment ~ candy + assessor,
    data = candies,
    contrasts = "contr.reference"
  )
  expect_equal(unname(mod_ref_full$contrasts), "contr.reference")
})

test_that("hdanova accepts named per-factor contrast vectors", {
  data(candies)

  mod_named <- hdanova(
    assessment ~ candy + assessor,
    data = candies,
    contrasts = c(candy = "reference", assessor = "sum")
  )

  expect_equal(unname(mod_named$contrasts), c("contr.reference", "contr.sum"))
  expect_equal(names(mod_named$contrasts), c("candy", "assessor"))
})

test_that("hdanova validates contrast vector values", {
  data(candies)

  expect_error(
    hdanova(
      assessment ~ candy + assessor,
      data = candies,
      contrasts = c(candy = "sum", assessor = "badcoding")
    ),
    "Each value in 'contrasts' must be one of"
  )

  expect_error(
    hdanova(
      assessment ~ candy + assessor,
      data = candies,
      contrasts = c("sum", NA_character_)
    ),
    "'contrasts' must be a character value or character vector without NA"
  )
})
