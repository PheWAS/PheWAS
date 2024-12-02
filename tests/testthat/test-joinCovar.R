phenotype_data <- createPhenotypes(sample_data$id.vocab.code.count, id.sex =
 sample_data$id.sex)
joinCovar(phenotype_data, sample_data$id.sex, sample_data$genotypes)
test_that('Base Works', {
  expect_equal(joinCovar(phenotype_data, sample_data$id.sex, sample_data$genotypes)$rsEXAMPLE[4], 0)})
test_that('Base Works', {
  expect_equal(joinCovar(phenotype_data, sample_data$id.sex, sample_data$genotypes)$`008.52`[4], FALSE)})

test_that('Base Works', {
  expect_equal(joinCovar(phenotype_data, sample_data$id.sex, sample_data$genotypes)$`010`[1], NA)})