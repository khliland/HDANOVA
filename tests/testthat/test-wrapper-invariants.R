test_that("ASCA-family wrappers return expected classes and core dimensions", {
  data(candies)

  asca_mod <- asca(assessment ~ candy + assessor, data = candies)
  expect_s3_class(asca_mod, "asca")
  expect_equal(nrow(asca_mod$Y), nrow(candies))

  apca_mod <- apca(assessment ~ candy + assessor, data = candies)
  expect_s3_class(apca_mod, "apca")
  expect_true(inherits(apca_mod, "asca"))
  expect_equal(nrow(apca_mod$Y), nrow(candies))

  apls_mod <- apls(assessment ~ candy + assessor, data = candies)
  expect_s3_class(apls_mod, "apls")
  expect_equal(nrow(apls_mod$Y), nrow(candies))

  limmpca_mod <- limmpca(assessment ~ candy + r(assessor), data = candies, pca.in = 3)
  expect_s3_class(limmpca_mod, "limmpca")
  expect_true(inherits(limmpca_mod, "asca"))
  expect_equal(nrow(limmpca_mod$Y), nrow(candies))
})

test_that("pcanova and msca wrappers expose expected structures", {
  data(candies)

  pcanova_mod <- pcanova(assessment ~ candy + assessor, data = candies, ncomp = 2)
  expect_s3_class(pcanova_mod, "pcanova")
  expect_true(length(pcanova_mod$anovas) >= 1)
  expect_equal(nrow(pcanova_mod$Y), nrow(candies))

  msca_mod <- msca(assessment ~ candy, data = candies)
  expect_s3_class(msca_mod, "msca")
  expect_equal(length(msca_mod$scores.within), nlevels(candies$candy))
  expect_equal(nrow(msca_mod$explvar.within), nlevels(candies$candy))
})
