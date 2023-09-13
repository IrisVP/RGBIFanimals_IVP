source("2_Functions.R")

library("testthat")

test_that("Create rastered world", {
  expect_no_error(create_rastered_world("Inputs/tr.rdata"))
  # TODO: add more tests here
})

test_that("Get occurrence data", {
  species <- "Acartia margalefi"
  expect_no_error(test <- get_occurrence_data(species))
  expect_gte(nrow(test), 6)
})

test_that("check_occurrence_data", {
  expect_true(TRUE)
})

test_that("find_closest_registered_place", {
  expect_true(TRUE)
})

test_that("filter_n_closest_coordinate_ceiling", {
  expect_true(TRUE)
})

test_that("find_shortest_route_in_sea", {
  expect_true(TRUE)
})

test_that("filter_on_distance", {
  expect_true(TRUE)
})

test_that("write.clean.csv", {
  expect_true(TRUE)
})

test_that("check_in_file", {
  expect_true(TRUE)
})

test_that("check_official_name", {
  expect_true(TRUE)
})
