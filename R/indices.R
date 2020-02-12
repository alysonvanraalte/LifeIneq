
#' @title ineq_variance
#' @description Calculate a lifetable column for the conditional variance in lifetable ages at death 
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in \code{MortalityLaws::MortalityLaw}). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package \code{ungroup} or with a penalized B-spline approach in package \code{MortalitySmooth}). 
#'   
#'
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector of the lifetable death distribution.
#' @param lx numeric. vector of the lifetable survivorship.
#' @param ex numeric. vector of remaining life expectancy.
#' @param ax numeric. vector of the average time spent in the age interval of those dying within the interval.
#'
#' @seealso 
#' \code{MortalityLaws::\link[MortalityLaws]{MortalityLaw}}
#' 
#' \code{ungroup::\link[ungroup]{pclm}}
#' 
#' \code{MortalitySmooth::\link[MortalitySmooth]{Mort1Dsmooth}}
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
#' @inheritParams ineq_variance
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
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
   sqrt(V) 
}




#' @title ineq_cov
#' @description Calculate a lifetable column for the conditional coefficient of variation in lifetable ages at death 
#'
#' @inheritParams ineq_variance
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
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
  sqrt(V) / (ex + age) 
}



#' @title ineq_edag
#' @description Calculate a lifetable column for the conditional life disparity (\eqn{e^\dagger}) of a population.  
#'
#' @inheritParams ineq_variance
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
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
  ex_average <- ex + ax / n * (explusone - ex)
  
  rev(cumsum(rev(dx * ex_average))) / lx 
}




#' @title ineq_H
#' @description Calculate a lifetable column for the quantity \emph{H}, generally referred to as either the lifetable entropy (Keyfitz, 1977) or the elasticity of life expectancy (Leser, 1955).
#'
#' @inheritParams ineq_variance
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
#' @references
#' \insertRef{Keyfitz1977}{LifeIneq}
#' \insertRef{Leser1955}{LifeIneq}
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
  ineq_edag(age, dx, lx, ex, ax) / (ex + age)
}




#' @title ineq_Theil
#' @description Calculate a lifetable column for the conditional Theil index of inequality in survivorship
#'
#' @inheritParams ineq_variance
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
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
    return(T1)
}




#' @title ineq_MLD
#' @description Calculate a lifetable column for the conditional mean log deviation index of inequality in survivorship
#'
#' @inheritParams ineq_variance
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
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
  
  MLD <- rep(NA,nages)
  for(i in 1:nages){
    MLD[i] <- sum(
                dx[i:nages]* (log (exAge[i]/axAge[i:nages]))
                   ) / lx[i]
    MLD <- ifelse(MLD<0,0,MLD) 
  }
  return(MLD)
}



#' @title ineq_Gini
#' @description Calculate a lifetable column for the conditional Gini coefficient of inequality in survivorship
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in \code{MortalityLaws::MortalityLaw}). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package \code{ungroup} or with a penalized B-spline approach in package \code{MortalitySmooth}). 
#' 
#'   The formula for calculating the Gini was taken from the Shkolnikov (2010) spreadsheet, and is a simplification of the formulas described in Shkolnikov (2003) and Hanada (1983).
#' 
#' @inheritParams ineq_variance
#' @inherit ineq_variance seealso
#'
#' @references 
#' \insertRef{hanada1983}{LifeIneq}
#' \insertRef{shkolnikov2003}{LifeIneq}
#' \insertRef{shkolnikov2010}{LifeIneq}
#' 
#' @export
#' 
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
  return(G)
}



#' @title ineq_AID
#' @description Calculate a lifetable column for the conditional absolute inter-individual difference in lifespan (AID)
#' 
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#' 
#' The formula for calculating the AID was taken from the citet{Shkolnikov2010} spreadsheet, and is a simplification of the formula described in citet{Shkolnikov2003}.
#' 
#'
#' @inheritParams ineq_variance
#' @inherit ineq_variance seealso
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
  return(AID)
}




#' @title ineq_LTquantile
#' @description Calculate quantiles of survivorship from a lifetable
#'
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
#'
#' @inheritParams ineq_variance
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
  lx   <- lx / lx[1]
  fit  <- spline(x = age, y = lx, n = n.grid)
  agei <- fit$x
  lxi  <- fit$y
  age_quantile <- agei[which.min(abs(lxi - quantile))]
  return(age_quantile)
}



#' @title ineq_iqr
#' @description Calculate the interquartile range survivorship age from a lifetable
#'
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
#' 
#' @inheritParams ineq_variance
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
  q1-q3
}



#' @title ineq_Cp
#' @description Calculate Kannisto's C-measures from a lifetable
#'
#' @details The age and lx vectors must be the same length. This function estimates the shortest distance between two ages containing p percent of the life table cohort's death. The mechanics behind the function are to fit a cubic spline through the survival curve to estimate surivorship between age intervals. If your data have an upper age bound lower than 110, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws').
#' 
#'  This function is currently very slow because the author is a better demographer than she is programmer. For ideas of how to improve this function, please let email vanraalte@@demogr.mpg.de.
#'  
#'   The concept behind Kannisto's C-measures is found in Kannisto (2000) 
#' 
#'
#' @param age numeric. vector of lower age bounds.
#' @param lx numeric. vector containing the lifetable survivorship. The starting radix (\eqn{\ell_{0}}) can be 100000, 1, or any value you like.
#' @param n.grid numeric. How many points do we want to interpolate? Defaults to 1000, which would work out to about a tenth of a year for a vector of life table survivors from age zero to 100+.
#' @param p numeric. What percent of the life table cohort do you want captured in the C measure? The default is 50.
#'
#'
#'@references 
#' \insertRef{kannisto2000}{LifeIneq}
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The shortest age range containing 50 percent of the deaths
#' C50 <- ineq_Cp(age=LT$Age,lx=LT$lx,p=50,n.grid=1000)
#' C50


ineq_Cp <- function(age,lx,p=50,n.grid=1000){
  # First fitting a spline through the lx values to get a finer grid of age and lx
  # A bit slow. I have default n=1000. If it's too slow, lower n.grid.
  require(dplyr)
  
  fit <- spline(x=age, y=lx, n=n.grid)
  xi <- fit$x
  yi <- fit$y
  
  # two columns with all combinations of possible ages
  age_comparison <- data.frame(t(combn(xi,2)))
  
  # finding the difference between all lx values on our larger grid of lxs 
  diff_lx <- combn(yi,2,diff)
  
  diffdf <- data.frame(Age1=age_comparison[,1],Age2=age_comparison[,2]) %>%
                mutate(diff_age=Age2-Age1,diff_lx=-diff_lx) %>%
                filter(diff_lx >= p / 100 * lx[1])
  
  mindiff <- which.min(diffdf$diff_age)
  Age1 <- diffdf[mindiff,]$Age1
  Age2 <- diffdf[mindiff,]$Age2
  Cp <- Age2-Age1
  return(Cp)
  
}





