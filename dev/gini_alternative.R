
library(tidyverse)
library(LifeIneq)


ineq_gini2 <- function(age, dx, ex, ax, distribution_type = c("aad","rle"), check = TRUE){
  
  distribution_type <- match.arg(distribution_type)
  # dx <- dx / sum(dx)
  
  age_constant <- if (distribution_type == "aad"){
    age_constant <- age
  } else {
    age_constant <- age * 0
  }
  if (check){
    my_args <- as.list(environment())
    check_args(my_args)
  }
  
  N        <- length(age)
  g_out    <- rep(0,N)
  dx       <- dx / sum(dx)
  axage    <- age + ax
  denom    <- ex + age_constant
  ad       <- outer(axage,axage,"-") * lower.tri(diag(N),TRUE)
  pd       <- outer(dx, dx, "*")

  for (i in 1:N){
    g_out[i] <- sum(abs(ad) * pd) / (denom[i])
    ad       <- ad[-1, -1, drop = FALSE]
    pd       <- pd[-1, -1, drop = FALSE]
    pd       <- pd / sum(pd)
  }
  
  return(g_out)
}

data(LT)
age <- LT$Age
dx  <- LT$dx 
lx  <- LT$lx
ex  <- LT$ex
ax  <- LT$ax

plot(age,ineq_gini2(age,dx,ex,ax))
lines(age,LifeIneq::ineq_gini(age, lx, ex, ax))
ineq_gini2(age,dx,ex,ax) - LifeIneq::ineq_gini(age, lx, ex, ax)
plot(age,ineq_gini2(age,dx,ex,ax) - ineq_gini(age, lx, ex, ax))
plot(ineq_gini2(age=LT$Age,dx=LT$dx,ex=LT$ex,ax=LT$ax,distribution_type = "rle"))
lines(ineq_gini2(age=LT$Age,dx=LT$dx,ex=LT$ex,ax=LT$ax,distribution_type = "aad"))

lines(ineq_gini(age=LT$Age,lx=LT$lx,ex=LT$ex,ax=LT$ax,distribution_type = "rle"),col="blue")
lines(ineq_gini(age=LT$Age,lx=LT$lx,ex=LT$ex,ax=LT$ax,distribution_type = "aad"),col="red")

ineq_gini2(age5,dx5,ex5,ax5)
lines(age5,ineq_gini2(age5,dx5,ex5,ax5),type='s')
plot(age,ineq_gini2(age,dx,ex,ax))
lines(age5,ineq_gini2(age5,dx5,ex5,ax5))


plot(age5,ineq_gini2(age5,dx5,ex5,ax5))
lines(age5,ineq_gini(age5,lx5,ex5,ax5))

plot(ineq_sd(LT$Age,LT$dx,LT$lx,LT$ex,LT$ax))
dx

HUN
library(tidyverse)

HUN %>% 
  group_by(Year) %>% 
  mutate(sd = ineq_sd(age = Age,
                      dx = dx,
                      lx= lx,ex=ex,ax=ax)) %>% 
  filter(is.na(sd))

HUN %>% 
  filter(Year == 1950, lx < 5)
