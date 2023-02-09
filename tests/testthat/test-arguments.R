context("test vector arguments")
# here we should see what happens if 0s or negatives are in inputs

# these are basically tests for our checks?

# age checks
test_that("check age for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = age_w_neg, 
         lx = LT$lx, 
         dx = LT$dx, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age_w_na, 
         lx = LT$lx, 
         dx = LT$dx, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age_w_z, 
         lx = LT$lx, 
         dx = LT$dx, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))}
)

# lx checks
test_that("check lx for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = 0:110, 
         lx = lx_w_neg, 
         dx = LT$dx, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         lx = lx_w_na, 
         dx = LT$dx, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         lx = lx_w_z, 
         dx = LT$dx, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))}
)


# ex checks
test_that("check ex for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = 0:110, 
         lx = LT$lx, 
         dx = LT$dx, 
         ex = ex_w_neg, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         lx = LT$lx, 
         dx = LT$dx, 
         ex = ex_w_na, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         lx = LT$lx, 
         dx = LT$dx, 
         ex = ex_w_z ,
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  }
)

# dx checks

test_that("check dx for negatives, NAs, and 0s", {
  expect_error(
    ineq(age = 0:110, 
         lx = LT$lx, 
         dx = dx_w_neg, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         lx = lx, 
         dx = dx_w_na, 
         ex = LT$ex, 
         ax = LT$ax, 
         method = "sd", 
         check = TRUE))
  # 0s are allowable for dx
 
  expect_no_error(
      ineq(age = 0:110, 
           lx = LT$lx, 
           dx = dx_w_z, 
           ex = LT$ex,
           ax = LT$ax, 
           method = "sd", 
           check = TRUE)
      )
  
  
    expect_no_error(
      ineq(age = 0:110, 
           lx = LT$lx, 
           dx = dx_w_z, 
           ex = LT$ex,
           ax = LT$ax, 
           method = "cov", 
           distribution_type = "achieved_age")
    )
  }
)

