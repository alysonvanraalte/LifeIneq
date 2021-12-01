

devtools$load_all()
data(LT)
measures <- c("mad", "var", "sd", "cov", "edag", "H", "theil", "mld", "gini", "aid")
results_current <- results_conditional <- matrix(0, ncol = nrow(meta), nrow = 111, dimnames = list(age=0:110, measures))

for (measure in measures){
  results_current[,measure] <- ineq(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax,method=measure)
}

for (measure in measures){
  for (a in 1:100){
    results_conditional[a, measure] <- ineq(age = LT$Age[1:(112-a)],
                                            dx = LT$dx[a:111],
                                            lx = LT$lx[a:111],
                                            ex = LT$ex[a:111],
                                            ax = LT$ax[a:111],
                                            method = measure)[1]
  }
}
results_current      <- results_current[1:100,]
results_conditional  <- results_conditional[1:100,]

colSums(abs(results_current - results_conditional)) %>% zapsmall()
# cov
# H
# theil
# mld