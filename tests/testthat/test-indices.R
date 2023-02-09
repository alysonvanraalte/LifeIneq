source("setup.R")

# test rectangular population against expectations.
test_that("rectangle boundary condition", {
  # dxr <- c(rep(0,110),1)
  # lxr <- c(rep(1,110),0)
  # exr <- c(rep(110,110),0)
  # axr <- rep(0,111) # either 0 or 1 here
  
expect_equal(ineq(age = age,
     dx = dxr,
     lx = lxr,
     ex = exr,
     ax = axr,
     method = "edag",
     check = FALSE)[-N],
     rep(0,N-1))
})

# expect length errors
# TODO: add length checking to ineq()
# expect_error(ineq(age = age,
#                   dx = dx,
#                   lx = lx,
#                   ex = ex,
#                   ax = ax[-1],
#                   method = "edag"))


# ------------------------------------------------------------------------#
# missing argument errors. We only check for these 
# preemtively in ineq(). ineq_xyz do not check args.
# can follow this pattern to do all the others. Note diff
# methods need different things, and non-lifetable columns
# tend to have default args.
test_that("missing args caught", {
  expect_error(ineq(age = age,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    method = "edag"))
  expect_error(ineq(age = age,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    method = "sd"))

})

# test wrapper equivalence  (pattern can be followed for other functions) #
# ------------------------------------------------------------------------#
test_that("wrapper equivalence", {
  # sd
  expect_equal(ineq_sd(age = age,
                       dx = dx,
                       lx = lx,
                       ex = ex,
                       ax = ax), 
               ineq(age =  age,
                    dx =  dx,
                    lx =  lx,
                    ex =  ex,
                    ax =  ax,
                    method = "sd"))
  # variance
  expect_equal(ineq_var(age = age,
                       dx = dx,
                       lx = lx,
                       ex = ex,
                       ax = ax), 
               ineq(age = age,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    ax = ax,
                    method = "var"))
  # just comma-separate more statements like this...
})

