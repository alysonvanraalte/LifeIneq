
#' @title ineq_variance
#' @description Calculate a lifetable column for the conditional variance in lifetable ages at death 
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#'   
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector of the lifetable death distribution.
#' @param lx numeric. vector of the lifetable survivorship.
#' @param ex numeric. vector of remaining life expectancy.
#' @param ax numeric. vector of the average time spent in the age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional variance in age at death 
#' V = ineq_variance(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The variance in age at death from birth
#' V[1]
#' # The variance in age at death conditional upon survival to age 10
#' V[11]

ineq_variance <- function(age, dx, lx, ex, ax){
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



#' @title ineq_sd
#' @description Calculate a lifetable column for the conditional standard deviation in lifetable ages at death 
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector containing the lifetable deaths distribution.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional variance in age at death 
#' S = ineq_sd(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The standard deviation in age at death from birth
#' S[1]
#' # The standard deviation in age at death conditional upon survival to age 10
#' S[11]


ineq_sd <- function(age, dx, lx, ex, ax){
   V <- ineq_variance(age, dx, lx, ex, ax)
   sqrt(V) # %>%
   #   round(2) 
}




#' @title ineq_cov
#' @description Calculate a lifetable column for the conditional coefficient of variation in lifetable ages at death 
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector containing the lifetable deaths distribution.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional coefficient of variation in age at death 
#' CoV = ineq_sd(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The coefficient of variation in age at death from birth
#' CoV[1]
#' # The coefficient of variation in age at death conditional upon survival to age 10
#' CoV[11]


ineq_cov <- function(age, dx, lx, ex, ax){
  V <- ineq_variance(age, dx, lx, ex, ax)
  sqrt(V) / (ex + age) # %>%
  #   round(2) 
}



#' @title ineq_edag
#' @description Calculate a lifetable column for the conditional life disparity (\eqn{e^\dagger}) of a population.  
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector containing the lifetable deaths distribution.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional life disparity of a population
#' edag = ineq_edag(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The life disparity from birth
#' edag[1]
#' # The life disparity conditional upon survival to age 10
#' edag[11]
#' \dontrun{
#' plot(0:110, edag, type='l')
#' }

ineq_edag <- function(age, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  
  stopifnot(age_length_equal)
  
  # length of the age interval
  n <- c(diff(age),1)
  explusone <- c(ex[-1],ex[length(age)])
  # the average remaining life expectancy in each age interval 
  # (as opposed to the beginning of the interval)
  # ends up being roughly half of the ex between ages
    ex_average <- ex + ax/n * (explusone-ex)
  
  rev(cumsum(rev(dx*ex_average))) / lx # %>%
  #   round(2) 
}




#' @title ineq_H
#' @description Calculate a lifetable column for the quantity \emph{H}, generally referred to as either the lifetable entropy citep{Keyfitz1977} or the elasticity of life expectancy citep{Leser1955}
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector containing the lifetable deaths distribution.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional \emph{H} values
#' H = ineq_H(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The \emph{H} from birth
#' H[1]
#' # The \emph{H} conditional upon survival to age 10
#' H[11]


ineq_H <- function(age, dx, lx, ex, ax){
  ineq_edag(age, dx, lx, ex, ax) / (ex + age) # %>%
  #   round(4) 
}




#' @title ineq_Theil
#' @description Calculate a lifetable column for the conditional Theil index of inequality in survivorship
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector containing the lifetable deaths distribution.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional Theil indices
#' Theil = ineq_Theil(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The Theil indices from birth
#' Theil[1]
#' # The Theil indices conditional upon survival to age 10
#' Theil[11]


ineq_Theil <- function(age, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  
  stopifnot(age_length_equal)
  
  axAge <- ax + age
  exAge <- ex + age
  nages <- length(age)
  
  T1 <- rep(NA,nages)
    for(i in 1:nages){
      T1[i] <- sum(
                dx[i:nages]*(axAge[i:nages]/exAge[i]*
                               (log (axAge[i:nages]/exAge[i])))
                   ) / lx[i]
      T1 <- ifelse(T1<0,0,T1)
      }
    return(round(T1,4))
}




#' @title ineq_MLD
#' @description Calculate a lifetable column for the conditional mean log deviation index of inequality in survivorship
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector containing the lifetable deaths distribution.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional MLD indices
#' MLD = ineq_MLD(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The MLD from birth
#' MLD[1]
#' # The MLD conditional upon survival to age 10
#' MLD[11]


ineq_MLD <- function(age, dx, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(dx),
                                length(lx),length(ex),
                                length(ax))
  
  stopifnot(age_length_equal)
  
  axAge <- ax + age
  exAge <- ex + age
  nages <- length(age)
  
  MLD1 <- rep(NA,nages)
  for(i in 1:nages){
    MLD1[i] <- sum(
                dx[i:nages]* (log (exAge[i]/axAge[i:nages]))
                   ) / lx[i]
    MLD1 <- ifelse(MLD1<0,0,MLD1) 
  }
  return(round(MLD1,4))
}



#' @title ineq_Gini
#' @description Calculate a lifetable column for the conditional Gini coefficient of inequality in survivorship
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). \ The formula for calculating the Gini was taken from the citet{Shkolnikov2010} spreadsheet, and is a simplification of the formulas described in citet{Shkolnikov2003; Hanada1983}.
#' 
#'
#' @param age numeric. vector of lower age bounds.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional Gini coefficients
#' G = ineq_Gini(age=LT$Age,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The Gini coefficient from birth
#' G[1]
#' # The Gini coefficient conditional upon survival to age 10
#' G[11]


ineq_Gini <- function(age, lx, ex, ax){
  age_length_equal <- all.equal(length(age),length(lx),
                                length(ex),length(ax))
  
  stopifnot(age_length_equal)
  
  
  nages <- length(age)
  # vector of the length of the age interval
  n <- c(diff(age),ax[nages])
  
  # squared survivorship
  lx2 <- lx^2 / lx[1]^2
  lx2plusn <- c(lx2[-1],0)
  
  # the expression that will be integrated and the Gini 
  intlx2 <- lx2plusn*n + ax*(lx2-lx2plusn)*n
  G <- 1 - rev(cumsum(rev(intlx2)))/(ex*lx2)
  G[G<0] <- 0
  return(round(G,4))
}



#' @title ineq_AID
#' @description Calculate a lifetable column for the conditional absolute inter-individual difference in lifespan (AID)
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#' 
#' The formula for calculating the AID was taken from the citet{Shkolnikov2010} spreadsheet, and is a simplification of the formula described in citet{Shkolnikov2003}.
#' 
#'
#' @param age numeric. vector of lower age bounds.
#' @param lx numeric. vector containing the lifetable survivorship.
#' @param ex numeric. vector containing remaining life expectancy.
#' @param ax numeric. vector containing average time spent in age interval of those dying within the interval.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional absolute inter-individual difference in lifespan
#' AID = ineq_AID(age=LT$Age,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The absolute inter-individual difference in lifespan from birth
#' AID[1]
#' # The absolute inter-individual difference in lifespan from age 10
#' AID[11]


ineq_AID <- function(age, lx, ex, ax){
  AID <- ineq_Gini(age=age,lx=lx,ex=ex,ax=ax) * ex
  return(round(AID,2))
}




#' @title ineq_LTquantile
#' @description Calculate quantiles of survivorship from a lifetable
#'
#' @details The age and lx vectors must be the same length. This function estimates any quantile of survivorship for the life table cohort by fitting a cubic spline through the survival curve. If your data have a lower upper age bound then that which corresponds to the quantile that you are interested in, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). 
#' 
#'
#' @param age numeric. vector of lower age bounds.
#' @param lx numeric. vector containing the lifetable survivorship. The starting radix (\eqn{\ell_{0}}) can be 100000, 1, or any value you like.
#' @param quantile numeric. value between zero and one indicating the quantile of interest.
#' @param n.grid numeric. How many points do we want to interpolate? Defaults to 1000, which would work out to about a tenth of a year for a vector of life table survivors from age zero to 100+.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The median age at survival
#' LTmedian <- ineq_LTquantile(age=LT$Age,lx=LT$lx,quantile=0.5,n.grid=1000)
#' LTmedian
#' # The age reached by 90% of the life table cohort (i.e. the top 10%)
#' ineq_LTquantile(age=LT$Age,lx=LT$lx,quantile=0.1,n.grid=1000) 
#' # The difference between the bottom 10 and top 10 percent of survival age
#' ineq_LTquantile(age=LT$Age,lx=LT$lx,quantile=0.1,n.grid=1000) -
#' ineq_LTquantile(age=LT$Age,lx=LT$lx,quantile=0.9,n.grid=1000)


ineq_LTquantile <- function(age,lx,quantile,n.grid=1000){
  lx <- lx/lx[1]
  fit <- spline(x=age, y=lx, n=n.grid)
  agei <- fit$x
  lxi <- fit$y
  age_quantile <- agei[which.min(abs(lxi-quantile))]
  return(round(age_quantile,2))
}



#' @title ineq_iqr
#' @description Calculate the interquartile range survivorship age from a lifetable
#'
#' @details The age and lx vectors must be the same length. This function estimates the interquartile range of the life table cohort by fitting a cubic spline through the survival curve. If your data have a lower upper age bound then that which corresponds to the top quartile of survivors, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). 
#' 
#'
#' @param age numeric. vector of lower age bounds.
#' @param lx numeric. vector containing the lifetable survivorship. The starting radix (\eqn{\ell_{0}}) can be 100000, 1, or any value you like.
#' @param n.grid numeric. How many points do we want to interpolate? Defaults to 1000, which would work out to about a tenth of a year for a vector of life table survivors from age zero to 100+.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The iqr range of survival age
#' LTiqr <- ineq_iqr(age=LT$Age,lx=LT$lx,n.grid=1000)
#' LTiqr


ineq_iqr <- function(age,lx,n.grid=1000){
  q1 <- ineq_LTquantile(age=age,lx=lx,quantile=0.25,n.grid=1000) 
  q3 <- ineq_LTquantile(age=age,lx=lx,quantile=0.75,n.grid=1000)
  round(q1-q3,2)
}




