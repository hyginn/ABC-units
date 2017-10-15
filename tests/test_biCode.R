# test_biCode.R
#

context("biCode() utility function tests")  # A set of tests for some
                                            # functionality

test_that("expected input is processed correctly", {  # Related expectations
  expect_equal(biCode("homo sapiens"), "HOMSA")
  expect_equal(biCode("[homo sapiens neanderthaliensis]"), "HOMSA")
  expect_equal(biCode(c("Phascolarctos cinereus", "Macropus rufus")),
               c("PHACI", "MACRU"))
})

test_that("unexpected input is managed", {
  expect_equal(biCode(""), ".....")
  expect_equal(biCode(" "), ".....")
  expect_equal(biCode("123 12"), ".....")
  expect_equal(biCode("h sapiens"), "H..SA")
})

test_that("NA values are preserved", {
  expect_true(is.na((biCode(NA))))
  expect_equal(biCode(c("first", NA, "last")),
               c("FIRST", NA, "LAST."))
})

test_that("Missing argument throws an error", {
  expect_error(biCode(), "argument \"s\" is missing, with no default")
})


# [END]
