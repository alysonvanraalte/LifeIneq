

# 1) test against real conditions, ie data likely to be encountered
# incl abridged, or strange age groups, what about NAs? Should NA and Inf
#  should they nix a whole operation or just the cell in which found?
# length equality, classes.
# 2) boundary condition testing ensures that implementations abide w theory
# 3) other checks for missing arguments




# NAs are OK, but shouldn't propagate?
# compare:
# ineq(age = age,
#      dx = dx,
#      lx = c(NA, lx[-1]),
#      ex = ex,
#      ax = ax,
#      method = "sd")
# # with:
# ineq(age = age,
#      lx = c(NA, lx[-1]),
#      ex = ex,
#      ax = ax,
#      method = "Gini")



# max inequality
# dxmax <- c(1,rep(0,109),1)/2
# 
# 

# max uniformity over age in ineq measures
# should be if rates are constant. Can test for that?

# test rectangular population against expectations.
test_that("rectangle boundary condition", {
expect_equal(ineq(age = 0:110,
     dx = dxr,
     lx = lxr,
     ex = exr,
     ax = axr,
     method = "edag")[-N], rep(0,N-1))
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
  expect_error(ineq(age = 0:110,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    method = "edag"))
  expect_error(ineq(age = 0:110,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    method = "sd"))

})

# test wrapper equivalence  (pattern can be followed for other functions) #
# ------------------------------------------------------------------------#
test_that("wrapper equivalence", {
  # sd
  expect_equal(ineq_sd(age = 0:110,
                       dx = dx,
                       lx = lx,
                       ex = ex,
                       ax = ax), 
               ineq(age = 0:110,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    ax = ax,
                    method = "sd"))
  # variance
  expect_equal(ineq_var(age = 0:110,
                       dx = dx,
                       lx = lx,
                       ex = ex,
                       ax = ax), 
               ineq(age = 0:110,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    ax = ax,
                    method = "var"))
  # just comma-separate more statements like this...
})

