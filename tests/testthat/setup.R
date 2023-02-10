data(LT)
data(LTabr)
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
# dxt <- rep(1,111)/111
# lxt <- rev(cumsum(rev(dxt)))
# axt <- rep(.5,111)
# Lxt <- (lxt + c(lxt[-1],0))/2
# Txt <- rev(cumsum(rev(Lxt)))
# ext <- Txt / lxt

# set up assorted defective inputs to catch errors:

age_w_z <- age_w_neg <- age_w_na <- age 
dx_w_z  <- dx_w_neg  <- dx_w_na  <- dx  
lx_w_z  <- lx_w_neg  <- lx_w_na  <- lx  
ex_w_z  <- ex_w_neg  <- ex_w_na  <- ex  
ax_w_z  <- ax_w_neg  <- ax_w_na  <- ax  

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

age5 <- LTabr$Age
dx5  <- LTabr$dx
lx5  <- LTabr$lx
ex5  <- LTabr$ex
ax5  <- LTabr$ax
