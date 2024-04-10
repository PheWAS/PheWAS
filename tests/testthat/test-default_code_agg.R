test_that("numeric", {
  expect_equal(default_code_agg(c(3, 5, 6, 9)), 23)
})

test_that("non numeric", {
  expect_equal(default_code_agg(c("hi", "bye")), 2)
})
