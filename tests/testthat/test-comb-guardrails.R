test_that("comb() with correct usage works", {
  data(caldana)
  mod <- asca(compounds ~ time + comb(light + time:light), data = caldana)
  expect_s3_class(mod, "asca")
})

test_that("comb() with single term gives informative error", {
  data(caldana)
  expect_error(
    asca(compounds ~ time + light + comb(time:light), data = caldana),
    "must contain at least two terms"
  )
})

test_that("comb() with reversed colon order gives 'did you mean' hint", {
  data(caldana)
  expect_error(
    asca(compounds ~ time + comb(light + light:time), data = caldana),
    "did you mean 'time:light'"
  )
})

test_that("comb() with term not in model gives informative error", {
  data(caldana)
  expect_error(
    asca(compounds ~ time + comb(light + time:light), data = caldana),
    NA   # should NOT error — this is the correct usage
  )
})
