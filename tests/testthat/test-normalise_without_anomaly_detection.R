
input_file <- testthat::test_path("test_data/example1.xlsx") #CHANGE to example file in your directory
normalised_df <- invisible(final_func(input_file))
row.names(normalised_df) <- NULL # row.names not written using openxlsx in normalisation.R

expected_file <- testthat::test_path("test_data/example1_normalised_without_anomaly_detection.xlsx")
expected_df <- data.frame(readxl::read_excel(expected_file),check.names = FALSE)

test_that("Neutralisation using defaults works", {
  expect_equal(normalised_df,expected_df)
})
