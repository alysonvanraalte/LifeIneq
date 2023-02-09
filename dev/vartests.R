
ineq_var <- function(age, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  
  stopifnot(age_length_equal)
  
  axAge <- age + ax
  exAge <- ex + age
  
  # TR: if magrittr is an acceptable dependency, this
  # might be easier to read for some people.
  # residsq <- dx*(axAge-exAge)^2
  # residsq %>%
  #   rev() %>%
  #   cumsum() %>%
  #   rev() %>%
  #   '/'(lx) %>%
  #   round(1) # Alyson added this
  # 
  rev(cumsum(rev(dx * (axAge - exAge)^2))) / lx
}
sum(dx * (axAge - exAge)^2) / lx[1]
library(devtools)
load_all()
data(LT)

ineq_var3 <- function(age, dx, lx, ex, ax){
  age <- age - age[1]
  n   <- length(age)
  out <- rep(NA, n)
  
  for (i in 1:n){
    axi    <- age[1:(n+1-i)] + ax[i:n]
    out[i] <- sum(dx[i:n] * (axi - ex[i])^2) / lx[i]
  }
  out
}


