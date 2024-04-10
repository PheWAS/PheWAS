test_that("multiplication works", {
  expect_equal(restrictPhecodesBySex(sample_data$id.vocab.code.count, sample_data$id.sex), restrict_phe_by_sex_test_1)
})
