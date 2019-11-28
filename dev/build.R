
library(devtools)

document()


# do a once-off library() action on our package structure.
# without actually having to install it or build it.
load_all()



?ineq_variance
?ineq_edag
?ineq_H
?ineq_Theil
?ineq_MLD
?ineq_Gini
?ineq_AID
?ineq_LTquantile
?ineq_iqr

T1 <- rep(NA,nages)
for(i in 1:nages){
  T1[i] <- dx[i]*(axAge[i]/exAge[i]*(log (axAge[i]/exAge[i])))
}
