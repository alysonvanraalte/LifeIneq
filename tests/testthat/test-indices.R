data(LT)
age <- 0:110
dx  <- LT$dx
lx  <- LT$lx
ex  <- LT$ex
ax  <- LT$ax

N   <- 111

# rectangular pop
dxr <- c(rep(0,110),1)
lxr <- c(rep(1,110),0)
exr <- c(rep(110,110),0)
axr <- rep(0,111) # either 0 or 1 here

# triangular pop
dxt <- rep(1,111)/111
lxt <- rev(cumsum(rev(dxt)))
axt <- rep(.5,111)
Lxt <- (lxt + c(lxt[-1],0))/2
Txt <- rev(cumsum(rev(Lxt)))
ext <- Txt / lxt

# ----------------------------------------
# Tests are thematically grouped, can follow the pattern

dxNA <- dx * NA

# NAs are OK, but shouldn't propagate?
# compare:
ineq(age = age,
     dx = dx,
     lx = c(NA, lx[-1]),
     ex = ex,
     ax = ax,
     method = "sd")
# with:
ineq(age = age,
     lx = c(NA, lx[-1]),
     ex = ex,
     ax = ax,
     method = "Gini")



# max inequality
# dxmax <- c(1,rep(0,109),1)/2
# 
# 

# max uniformity over age in ineq measures
# should be if rates are constant. Can test for that?

# test rectangular population against expectations.
test_that("rectangle boundary condition", {
expect_equal(ineq(age = age,
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
               ineq(age = age,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    ax = ax,
                    method = "sd"))
  # variance
  expect_equal(ineq_variance(age = age,
                       dx = dx,
                       lx = lx,
                       ex = ex,
                       ax = ax), 
               ineq(age = age,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    ax = ax,
                    method = "variance"))
  # just comma-separate more statements like this...
})

# ------------------------------------------------------------------------#
# test wrapper equivalence  (pattern can be followed for other functions) #
# ------------------------------------------------------------------------#
test_that("wrapper equivalence for sd", {
  expect_equal(ineq_sd(age = age,
                       dx = dx,
                       lx = lx,
                       ex = ex,
                       ax = ax), 
               ineq(age = age,
                    dx = dx,
                    lx = lx,
                    ex = ex,
                    ax = ax,
                    method = "sd"))
})

