library(PheWASmaps)
test_that("Code is not class character)", {
  expect_error(mapCodesToPhecodes(sample_data_no_chara$id.vocab.code.count, vocabulary.map = phecode_map, rollup.map = phecode_rollup_map), "Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")
})

test_that("Null Vocab Map)", {
  expect_error(mapCodesToPhecodes(sample_data$id.vocab.code.count, vocabulary.map = NULL, rollup.map = phecode_rollup_map), "Phecode mapping was not requested, but the vocabulary_id of all codes is not 'phecode'")
})
test_that('Base Case Works', {
  expect_equal((mapCodesToPhecodes(sample_data$id.vocab.code.count, vocabulary.map = phecode_map, rollup.map = phecode_rollup_map)), mCtP_test_1)})