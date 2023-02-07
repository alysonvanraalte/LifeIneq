
# here we should see what happens if 0s or negatives are in inputs

data(LT)
age_w_z <- age_w_neg <- age_w_na <-age <- 0:110
dx_w_z  <- dx_w_neg  <- dx_w_na  <- dx  <- LT$dx
lx_w_z  <- lx_w_neg  <- lx_w_na  <- lx  <- LT$lx
ex_w_z  <- ex_w_neg  <- ex_w_na  <- ex  <- LT$ex
ax_w_z  <- ax_w_neg  <- ax_w_na  <- ax  <- LT$ax

# with negatives/ all categorical errors
age_w_neg[1]  <- -1
dx_w_neg[10]  <- -1
lx_w_neg[10]  <- -1
ex_w_neg[111] <- -1
ax_w_neg[10]  <- -.01

# with 0s in bad places/ monotonic and accidental errors
age_w_z[2]    <- 0  # non-monotonic error
dx_w_z[10]    <- 0 # no error
lx_w_z[10]    <- 0 # non-monotonic error
ex_w_z[10]    <- 0 # non-monotonic error and denom inf for some
ax_w_z[111]   <- 0 # not sure

# with NAs
age_w_na[1]   <- NA 
dx_w_na[10]   <- NA 
lx_w_na[10]   <- NA 
ex_w_na[10]   <- NA 
ax_w_na[10]   <- NA 

# these are basically tests for our checks?

# age checks
test_that("check age for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = age_w_neg, 
         lx = lx, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age_w_na, 
         lx = lx, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age_w_z, 
         lx = lx, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))}
)

# lx checks
test_that("check lx for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = age, 
         lx = lx_w_neg, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age, 
         lx = lx_w_na, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age, 
         lx = lx_w_z, 
         dx = dx, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))}
)


# ex checks
test_that("check ex for negatives, NAs, and irregular order", {
  expect_error(
    ineq(age = age, 
         lx = lx, 
         dx = dx, 
         ex = ex_w_neg, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age, 
         lx = lx, 
         dx = dx, 
         ex = ex_w_na, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age, 
         lx = lx, 
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
    ineq(age = age, 
         lx = lx, 
         dx = dx_w_neg, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  expect_error(
    ineq(age = age, 
         lx = lx, 
         dx = dx_w_na, 
         ex = ex, 
         ax = ax, 
         method = "sd", 
         check = TRUE))
  # 0s are allowable for dx
  expect_success(
    expect_type(
      ineq(age = age, 
           lx = lx, 
           dx = dx_w_z, 
           ex = ex,
           ax = ax, 
           method = "sd", 
           check = TRUE),'double'))
  
  expect_success(
    expect_type(
      ineq(age = age, 
           lx = lx, 
           dx = dx_w_z, 
           ex = ex,
           ax = ax, 
           method = "cov", 
           distribution_type = "achieved_age",
           check = TRUE),'double'))
  }
)
