library(tidyverse)
library(HMDHFDplus)
#library(ungroup)
#library(DemoDecomp)
options(scipen=3)

# HMD stuff
myHMDusername <- "vanraalte@demogr.mpg.de"
myHMDpassword <- "1538657766"

# an HMD life table as example data
LT <- readHMD("C:\\hmd_statistics\\lt_female\\fltper_5x1\\CAN.fltper_5x1.txt") %>%
        filter(Year==max(Year))

# let's assume our input is mx, and output an HMD-like life table.

# First, the columns that we will need (probably not all)
Age <- LT$Age; mx <- LT$mx; qx <- LT$qx;
ax <- LT$ax; lx <- LT$lx; dx  <- LT$dx; 
Lx <- LT$Lx; Tx <- LT$Tx; ex <- LT$ex

# number of age classes
nages <- dim(LT)[1]

# length of age interval (probably not needed if functions are built on single-year age)
n <- c(diff(Age), ax[length(Age)])
       
# average age at death of those who died in each age interval
axAge <- Age+ax

# average age at death conditional upon surviving to the age interval 
exAge <- Age+ex

#------------------------------------------------------------------------------#
#----- indices of variability given a life table with single age intervals ----#
#------------------------------------------------------------------------------#

#----------------- variance, standard deviation, coefficient of variation

V <- rev(cumsum(rev(dx*(axAge-exAge)^2)))/lx
S <- sqrt(V)
CoV <- S/exAge

# checking
plot(Age,V/V[1],t="l")
lines(Age,S/S[1],col="red")
lines(Age,CoV/CoV[1],col="blue")

#------------------------ edag, H

ex.av <- ex+ax*(c(ex[-1],ex[nages])-ex)
edag <- rev(cumsum(rev(dx*ex.av)))/lx
H <- edag/exAge

lines(Age,edag/edag[1],col="green")
lines(Age,H/H[1],col="darkorange")

#------------------------- Gini (G), AAD



#------------------------- Theil, MLD
# these are also currently functions that give only starting age values

Theil <- sum(dx*(axAge/exAge[1]*(log(axAge/exAge[1]))))/lx[1]
MLD <- sum(dx*(log(exAge[1]/axAge)))/lx[1]

# For all ages

T1 <- rep(NA,length(Age))
for(i in 1:length(Age)){
  T1[i] <- sum(dx[i:nages]*(axAge[i:nages]/exAge[i]*(log(axAge[i:nages]/exAge[i]))))/lx[i]
  T1 <- ifelse(T1<0,0,T1)
  }

MLD1 <- rep(NA,length(Age))
for(i in 1:length(Age)){
  MLD1[i] <- sum(dx[i:nages]*log(exAge[i]/axAge[i:nages])) / lx[i]
  MLD1 <- ifelse(MLD1<0,0,MLD1)
  }

lines(Age,T1/T1[1],col="brown")
lines(Age,MLD1/MLD1[1],col="grey40")


#------------------------- IQR, Median

# a function which fits a spline through the survivorship curve 
# to get a finer grid of age and lx

# for now these are only programmed to give me one result (i.e. iqr at birth). 
# would need to change if it is desired to have conditional iqrs with each age.

intFUN <- function(x, y, n, n.grid=10000){
  fit <- spline(x=x, y=y, n=n.grid)
  xi <- fit$x
  yi <- fit$y
  xn <- xi[which.min(abs(yi-n))]
  return(xn)
}

iqr <- intFUN(x=Age, y=lx, n=25000) - intFUN(x=Age, y=lx, n=75000)  
median <- intFUN(x=Age, y=lx, n=50000)


#-------------- C50 ----------------------------------------#


# First fitting a spline through the lx values to get a finer grid of age and lx
# A bit slow. I have default n=1000, but for the one example I used I got the 
# same result for n=500. For n=250, I got the same C50, 
# but the lower and upper ages were both 0.2 years higher.

fit <- spline(x=Age, y=lx, n=1000)
xi <- fit$x
yi <- fit$y

# finding the difference between all ages and lx values on our larger grid of lxs
gg <- combn(xi,2,diff)
hh <- combn(yi,2,diff)

agecomp <- data.frame(t(combn(xi,2))) 
diffdf <- data.frame(Age1=agecomp[,1],Age2=agecomp[,2],diffage=gg,difflx=-hh) %>%
            filter(difflx>=50000)

mindiff <- which.min(diffdf$diffage)
Age1 <- diffdf[mindiff,]$Age1
Age2 <- diffdf[mindiff,]$Age2
C50 <- Age2-Age1







