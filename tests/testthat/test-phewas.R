 test_that("No Data Causes Error", {
   expect_error(phewas(), "A data frame must be supplied in 'data'")
 })
 test_that("Translate False causes - phecode)", {
   expect_error(phewas(data = test_data_phenotype_2), "Either phenotypes or outcomes must be passed in.")
})
test_that("Translate False causes - phecode)", {
  expect_error(phewas(names(test_data_phenotype_1)[-1], covariates= 'sex', data = test_data_phenotype_2), "Either genotypes or predictors must be passed in.")
})
test_that('Base Case Works', {
expect_equal(phewas(names(test_data_phenotype_1)[-1], genotypes = c('rsEXAMPLE'), covariates= 'sex', data = test_data_phenotype_2)$p[14], 0.90935256)})
#test_that('Base Case Works', {
#expect_equal(phewas_ext(names(test_data_phenotype_1)[-1], genotypes = c('rsEXAMPLE'), covariates= 'sex', data = test_data_phenotype_2, method = 'lrt'), test_phewas_2)})
test_that('Base Case Works', {
expect_equal(phewas(names(test_data_phenotype_1)[-1], genotypes = c('rsEXAMPLE'), covariates= 'sex', data = test_data_phenotype_2, method = 'logistf')$p[14], 0.91162193)})