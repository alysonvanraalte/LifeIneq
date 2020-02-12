
library(devtools)

document()
check()

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
?ineq_Cp


install_github("alysonvanraalte/LifeIneq")
