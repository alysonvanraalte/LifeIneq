source("setup.R")
# here we should see what happens if 0s or negatives are in inputs

# these are basically tests for our checks?

# age checks
test_that("check age for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = age_w_neg, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age_w_na, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age_w_z, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))}
)

# lx checks
test_that("check lx for negatives, NAs, and irregular order", {
  expect_error(
    ineq_iqr(age = 0:110, 
         lx = lx_w_neg,
         check = TRUE))
  expect_error(
    ineq_iqr(age = 0:110, 
             lx = lx_w_na,
             check = TRUE))
  expect_error(
    ineq_iqr(age = 0:110, 
             lx = lx_w_z,
             check = TRUE))}
)


# ex checks
test_that("check ex for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = 0:110, 
         dx = dx, 
         ex = ex_w_neg, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         dx = dx, 
         ex = ex_w_na, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         dx = dx, 
         ex = ex_w_z ,
         ax = ax, 
         method = "sd", 
         check = TRUE))
  }
)

# dx checks

test_that("check dx for negatives, NAs, and 0s", {
  expect_error(
    ineq(age = 0:110, 
         dx = dx_w_neg, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = 0:110, 
         dx = dx_w_na, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  # 0s are allowable for dx
 
  expect_no_error(
      ineq(age = 0:110, 
           dx = dx_w_z, 
           ex = ex,
           ax = ax, 
           method = "sd", 
           check = TRUE)
      )
  
  
   
  }
)

test_that("distribution_type works", {
  
  # cov
  expect_no_error(
    ineq(age = 0:110, 
         dx = dx, 
         ex = ex,
         ax = ax, 
         method = "cov", 
         distribution_type = "aad")
  )
  expect_no_error(
    ineq(age = 0:110, 
         dx = dx, 
         ex = ex,
         ax = ax, 
         method = "cov", 
         distribution_type = "rl")
  )
  # theil
  expect_no_error(
    ineq(age = 0:110, 
         dx = dx, 
         ex = ex,
         ax = ax, 
         method = "theil", 
         distribution_type = "aad")
  )
  expect_no_error(
    ineq(age = 0:110, 
         dx = dx, 
         ex = ex,
         ax = ax, 
         method = "theil", 
         distribution_type = "rl")
  )
    
})