test_that("test pre_check", {
  expect_message(pre_check(counts_input, group_list, tempdir()),"have done")
})
