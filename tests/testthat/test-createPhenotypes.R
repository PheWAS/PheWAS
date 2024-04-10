library(PheWASmaps)
test_that("No id.sex throws Warning (to be discussed)", {
  expect_warning(createPhenotypes(sample_data$id.vocab.code.count), 'It is recommended to provide id.sex information to help address spurious sex-specific associations.')
})
test_that("Translate False causes - phecode)", {
  expect_error(createPhenotypes(sample_data$id.vocab.code.count, translate = F), "Codes are not translated and vocab is not 'phecode' for every row, but exclusions are to be applied. Ensure that the code column has only phecodes or disable add.phecode.exclusions for accurate results.")
})
test_that('Base Case Works', {
  expect_equal(createPhenotypes(sample_data$id.vocab.code.count, id.sex = sample_data$id.sex), create_pheno_test_1)})
test_that('Exclusions = False works', { 
  expect_equal(createPhenotypes(sample_data$id.vocab.code.count, id.sex = sample_data$id.sex, add.phecode.exclusions = F), create_pheno_test_2)})
test_that('Min Case Count Works', {
  expect_equal(createPhenotypes(sample_data$id.vocab.code.count, id.sex = sample_data$id.sex, min.code.count = 3), create_pheno_test_3)})
test_that('No Sex Works', {
  expect_equal(createPhenotypes(sample_data$id.vocab.code.count), create_pheno_test_4)})
