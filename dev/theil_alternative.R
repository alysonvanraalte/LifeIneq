ineq_theil2 <- function(age, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  
  stopifnot(age_length_equal)
  
  N     <- length(age)
  
  T1 <- rep(NA,N)
  for(i in 1:N){
    axi <- age[1:(N+1-i)] + ax[i:N]
    T1[i] <- sum(dx[i:N] * (axi / ex[i] * log( axi / ex[i]))) / lx[i]
  }
  T1[T1 < 0] <- 0
  return(T1)
}

devtools::load_all()
data(LT)
plot(0:110,ineq_theil(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax),ylim=c(0,.35))
lines(0:110,ineq_theil2(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax))

# coincides with the conditional loop (being born in age x as if it were 0)
# ineq_theil2(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax) -
#   results_conditional[,"theil"]


ineq_mld2 <-  function(age, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  stopifnot(age_length_equal)
  
  N     <- length(age)
  
  MLD <- rep(NA, N )
  for(i in 1: N ){
    axi <- age[1:(N+1-i)] + ax[i:N]
    MLD[i] <- sum(
      dx[i:N]* (log (ex[i]/axi))
    ) / lx[i]
    
  }
  MLD[MLD < 0] <- 0
  return(MLD)
}