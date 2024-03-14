test_that("non chara phenotypes", {
  expect_error(phewasManhattan(test_phewas_no_chara), "Non-character phenotypes passed in, so an accurate
             phecode mapping is not possible.")
})
test_that("not padded phenotypes", {
  expect_error(phewasManhattan(test_phewas_non_padded), "Phenotypes with length <3 observed,ensure they are are 0-padded")
})

#test_that("Working Test Phewas", {
# expect_equal(phewasManhattan(test_phewas_1k), test_phewas_manhattan_1k)

#})
