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
  # test file that is available
  species <- "Acartia margalefi"
  filename <- paste0("OccurrenceData/", species, ".csv")
  expect_true(file.exists(filename))
  expect_no_error(test <- check_occurrence_data(species))
  expect_gte(nrow(test), 6)
  # test animal not existing in the data set
  species <- "Desulfurococcus amylolyticus"
  filename <- paste0("OccurrenceData/", species, ".csv")
  expect_false(file.exists(filename))
  expect_no_error(test <- check_occurrence_data(species))
  expect_true(file.exists(filename))
  expect_gte(nrow(test), 4)
})

test_that("find_closest_registered_place", {
  species <- "Acartia margalefi"
  coordinates <- data.frame("Observatory.ID" = "BelgianCoast",
                            "Longitude" = 2.79178508,
                            "Latitude" = 51.6533567)
  tr <- create_rastered_world("Inputs/tr.rdata")
  # expect_no_condition(find_closest_registered_place(species, coordinates,
  #                                                   tr, "test_outputs/tests.csv"))
  # TODO: add more tests here
  expect_true(TRUE)
})

test_that("filter_n_closest_coordinate_ceiling", {
  expect_true(TRUE)
  # TODO: add more tests here
})

test_that("find_shortest_route_in_sea", {
  row <- data.frame("species" = "Acartia margalefi",
                    "location" = "Roscoff",
                    "value" = 1570)
  samplelocation <- data.frame("location" = "Roscoff",
                               "Longitude" = -3.9665,
                               "Latitude" = 48.7175)
  occurrence_data <- check_occurrence_data(row[1])
  tr <- create_rastered_world("Inputs/tr.rdata")
  filename <- "test_outputs/DistanceOverSeaFORTESTS.csv"
  expect_false(file.exists(filename))
  expect_no_error(find_shortest_route_in_sea(samplelocation, occurrence_data,
                                             tr, row, filename))
  expect_true(file.exists(filename))
  expect_warning(find_shortest_route_in_sea(samplelocation, occurrence_data,
                                            tr, row, filename))
  # TODO: add more tests here
})

test_that("filter_on_distance", {
  expect_true(TRUE)
  # TODO: add more tests here
})

test_that("write.clean.csv", {
  expect_true(TRUE)
  # TODO: add more tests here
})

test_that("check_in_file", {
  expect_true(TRUE)
  # TODO: add more tests here
})

test_that("check_official_name", {
  expect_true(TRUE)
  # TODO: add more tests here
})

rm(list = ls())
