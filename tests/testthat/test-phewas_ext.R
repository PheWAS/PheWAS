test_that("No Data Causes Error", {
  expect_error(phewas_ext(), "A data frame must be supplied in 'data'")
})
test_that("Translate False causes - phecode)", {
  expect_error(phewas_ext(data = test_data_phenotype_2), "Either phenotypes or outcomes must be passed in.")
})
test_that("Translate False causes - phecode)", {
  expect_error(phewas_ext(names(test_data_phenotype_1)[-1], covariates= 'sex', data = test_data_phenotype_2), "Either genotypes or predictors must be passed in.")
})
test_that('Base Case Works', {
  expect_equal(phewas_ext(names(test_data_phenotype_1)[-1], genotypes = c('rsEXAMPLE'), covariates= 'sex', data = test_data_phenotype_2), test_phewas_1)})
test_that('Base Case Works', {
  expect_equal(phewas_ext(names(test_data_phenotype_1)[-1], genotypes = c('rsEXAMPLE'), covariates= 'sex', data = test_data_phenotype_2, method = 'lrt'), test_phewas_2)})
test_that('Base Case Works', {
  expect_equal(phewas_ext(names(test_data_phenotype_1)[-1], genotypes = c('rsEXAMPLE'), covariates= 'sex', data = test_data_phenotype_2, method = 'logistf'), test_phewas_3)})