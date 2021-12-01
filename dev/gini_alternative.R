
library(tidyverse)
library(LifeIneq)


ineq_gini2 <- function(age, dx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(ex),length(ax))
  
  stopifnot(age_length_equal)
  N <- length(age)
  g_out <- rep(0,N)
  dx       <- dx / sum(dx)
  ad       <- outer(age+ax,age+ax,"-") * lower.tri(diag(N),TRUE)
  pd       <- outer(dx, dx, "*")
  g_out[1] <- sum(abs(ad) * pd) / ex[1]
  for (i in 2:N){
    dx       <- dx[-1]
    age      <- age[-length(age)]
    ax       <- ax[-1]
    dx       <- dx / sum(dx)
    ad       <- outer(age+ax,age+ax,"-") * lower.tri(diag(N+1-i),TRUE)
    pd       <- outer(dx, dx, "*")
    g_out[i] <- sum(abs(ad) * pd) / ex[i]
  }
  
  return(g_out)
}

data(LT)
age <- LT$Age
dx <- LT$dx 
lx <- LT$lx
ex <- LT$ex
ax <- LT$ax

plot(age,ineq_gini2(age,dx,ex,ax))
lines(age,LifeIneq::ineq_gini(age, lx, ex, ax))
ineq_gini2(age,dx,ex,ax) - LifeIneq::ineq_gini(age, lx, ex, ax)
plot(age,ineq_gini2(age,dx,ex,ax) - LifeIneq::ineq_gini(age, lx, ex, ax))
