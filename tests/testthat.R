library(testthat)
library(LifeIneq)

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


#test_check("LifeIneq")
