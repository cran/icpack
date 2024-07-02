context("Checks for get_input_icfit function\n")

test_that("entry should be before left end of intervals", {
  # Generate data
  formula <- Surv(x, y, type="interval2") ~ 1
  data <- data.frame(y = c(0, 1, 2, 3, 4) + 1, x = c(0, 1, 2, 3, 4))
  entry <- c(1, 2, 3, 4, 5)
  
  # Assert that function stops when entry is greater than Y[, 1]
  expect_error(get_input_icfit(formula, data, entry), "entry should be before left end of intervals")
})

context("Input checks for icfit\n")

test_that("response needs a formula", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  expect_error(get_input_icfit(Surv(left, right, type="interval2"), data=fake),
               "object 'left' not found")
})



test_that("dimensions of y and X are correct, no covariates", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  tmp <- get_input_icfit(Surv(left, right, type="interval2") ~ 1, data=fake)
  expect_equal(dim(tmp$Ymat), c(3,3))
  expect_equal(dim(tmp$X), c(3,0))
})

test_that("dimensions of y and X are correct, one covariate", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  tmp <- get_input_icfit(Surv(left, right, type="interval2") ~ x1, data=fake)
  expect_equal(dim(tmp$Ymat), c(3,3))
  expect_equal(dim(tmp$X), c(3,1))
})

test_that("dimensions of y and X are correct, two covariates", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  tmp <- get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake)
  expect_equal(dim(tmp$Ymat), c(3,3))
  expect_equal(dim(tmp$X), c(3,2))
})

test_that("warning (from survival package) when left is larger than right", {
  fake <- data.frame(left=c(0,8,5), right=c(6, 2, 9), x1=c(1,2,3), x2=c(4,3,6))
  expect_warning(get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake),
                 "Invalid interval: start > stop, NA created")
})

test_that("observations deleted for which left is larger than right", {
  fake <- data.frame(left=c(0,8,5), right=c(6, 2, 9), x1=c(1,2,3), x2=c(4,3,6))
  tmp <- get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake)
  expect_equal(dim(tmp$Ymat), c(3,3))
  expect_equal(dim(tmp$X), c(2,2))
})

test_that("error when entry is not correct length", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  expect_error(get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake, entry=c(1,2)),
               "Assertion on 'entry' failed: Must have length 3, but has length 2.")
})

test_that("error when entry is character vector", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  expect_error(get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake, entry=c("a","b","c")),
               "undefined columns selected")
})

test_that("error when entry is not numeric", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  expect_error(get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake, entry=c(TRUE,FALSE,TRUE)),
               "Assertion on 'entry' failed: Must be of type 'numeric', not 'logical'.")
})

test_that("error when entry is larger than left", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  expect_error(get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake, entry=c(1,3,5)),
               "entry should be before left end of intervals")
})

test_that("entry is inserted before left and right", {
  fake <- data.frame(left=c(0,2,5), right=c(6, 8, 9), x1=c(1,2,3), x2=c(4,3,6))
  tmp <- get_input_icfit(Surv(left, right, type="interval2") ~ x1+x2, data=fake, entry=c(0,1,2))
  expect_equal(c(tmp$Ymat), c(0, 1, 2, 0, 2, 5, 6, 8, 9))
})