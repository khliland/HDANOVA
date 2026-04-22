test_that("predict.limmpca returns limmpca object", {
  data(candies)

  test_idx <- seq(3, nrow(candies), by = 3)
  train_idx <- setdiff(seq_len(nrow(candies)), test_idx)
  candies_train <- candies[train_idx, ]
  candies_test <- candies[test_idx, ]

  mod <- limmpca(assessment ~ candy + r(assessor),
                 data = candies_train,
                 pca.in = 3,
                 REML = TRUE)

  pred <- predict(mod, newdata = candies_test)

  expect_true(inherits(pred, "limmpca"))
  expect_true(inherits(pred, "asca"))
  expect_equal(nrow(pred$Y), nrow(candies_test))
})
