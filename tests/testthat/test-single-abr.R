# test scales
source("setup.R")
?ineq_var
?ineq_edag
?ineq_H
?ineq_theil
?ineq_mld # strictly relative
?ineq_gini # (have shortfall and attainment)
?ineq_aid
?ineq_quantile
?ineq_iqr
?ineq_cp

# single
ineq_var(age,dx,lx,ex,ax)[1]
ineq_var(age5,dx5,lx5,ex5,ax5)[1]

ineq_sd(age,dx,lx,ex,ax)[1]
ineq_sd(age5,dx5,lx5,ex5,ax5)[1]

ineq_edag(age,dx,lx,ex,ax)[1]
ineq_edag(age5,dx5,lx5,ex5,ax5)[1]

ineq_theil(age,dx,ex,ax)
ineq_theil(age5,dx5,lx5,ex5,ax5)[1]
load_all()

# This is way off; we now issue message if ages not single
plot(age,ineq_mld(age,dx,lx,ex,ax))
lines(age5,ineq_mld(age5,dx5,lx5,ex5,ax5))

ineq_gini(age,lx,ex,ax)[1]
ineq_gini(age5,lx5,ex5,ax5)[1]

ineq_aid(age,lx,ex,ax)[1]
ineq_aid(age5,lx5,ex5,ax5)[1]

plot(age,ineq_quantile(age,lx, .5))
lines(age5,ineq_quantile(age5,lx5, .5))

ineq_iqr(age,lx)
ineq_iqr(age5,lx5)


ineq_cp(age,lx)
ineq_cp(age5,lx5)
