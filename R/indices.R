
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




#' @title calculate a survivorship quantile
#' @description Calculate quantiles of survivorship from a lifetable
#'
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
#'
#' @inheritParams ineq_variance
#' @param quantile numeric. value between zero and one indicating the quantile of interest.
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The median age at survival
#' LTmedian <- ineq_quantile(age=LT$Age,lx=LT$lx,quantile=0.5)
#' LTmedian
#' # The age reached by 90% of the life table cohort (i.e. the top 10%)
#' ineq_quantile(age=LT$Age,lx=LT$lx,quantile=0.1) 
#' # The difference between the bottom 10 and top 10 percent of survival age
#' ineq_quantile(age=LT$Age,lx=LT$lx,quantile=0.1) -
#' ineq_quantile(age=LT$Age,lx=LT$lx,quantile=0.9)

ineq_quantile <- function(age,lx,quantile=.5){
  lx   <- lx / lx[1]
  splinefun(age~lx)(quantile)
}


#' @title ineq_iqr
#' @description Calculate the interquartile range survivorship age from a lifetable. Other quantile ranges can also be calculated.
#'
#' @inherit ineq_variance details
#' @inherit ineq_variance seealso
#' 
#' @inheritParams ineq_variance
#' @param upper numeric. upper survival quantile, default .75
#' @param lower numeric. lower survival quantile, defauly .25
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The iqr range of survival age
#' LTiqr <- ineq_iqr(age=LT$Age,lx=LT$lx)
#' LTiqr


ineq_iqr <- function(age, lx, upper = .75, lower = .25){
  q1 <- ineq_quantile(age = age, lx = lx, quantile = lower) 
  q3 <- ineq_quantile(age = age, lx = lx, quantile = upper)
  q1 - q3
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
#' @inheritParams ineq_variance
#' @param p numeric. What proportion of the life table cohort do you want captured in the C measure? The default is .5
#'
#'
#'@references 
#' \insertRef{kannisto2000}{LifeIneq}
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The shortest age range containing half of the deaths
#' (C50 <- ineq_Cp(age=LT$Age,lx=LT$lx,p=.5))

ineq_Cp <- function(age, lx, p = .5){

  stopifnot(length(age) == length(lx))
  stopifnot(p <= 1)
  
  # span minimizer function
  ineq_Cp_min <- function(par = 60, fun, funinv, p = .5){
    q1 <- fun(par)
    q2 <- q1 - p
    # age span
    funinv(q2) - par
  }

  lx  <- lx / lx[1]
  
  a2q <- splinefun(x = age, y = lx, method = "monoH.FC")
  q2a <- splinefun(x = lx, y = age, method = "monoH.FC")
  
  # lower age can't be higher than this:
  tops <- q2a(p*1.0001)
  # first optimize for lower age by minimizing
  # the interval capturing 
  age1 <- optimize(ineq_Cp_min, 
           interval = c(min(age), tops),
           fun = a2q, 
           funinv = q2a,
           p = p)
  # this is the age span
  age1$objective
}

# -------------------------------------
# wrapper function

#' @title calculate a lifespan inequality measure
#' @description Choose from variance, standard deviation \code{sd},IQR, AID, Gini, edagger, or Kannisto's Cp.
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector of the lifetable death distribution.
#' @param lx numeric. vector of the lifetable survivorship.
#' @param ex numeric. vector of remaining life expectancy.
#' @param ax numeric. vector of the average time spent in the age
#' @param method one of \code{c("variance","sd","iqr","AID","Gini","MLD","edag","Cp")}
#' @param ... other optional arguments used by particular methods.
#'
#' @export

ineq <- function(age, dx, lx, ex, ax, method = c("variance","sd","iqr","AID","Gini","MLD","edag","Cp"),...){
  
  # make sure just one
  method         <- match.arg(method)
  # fun is now the function we need
  fun            <- match.fun(paste0("ineq_", method))
  
  # what do we need and what do we have?
  need_args      <- names(formals(fun))
  have_args      <- c(as.list(environment()), list(...))
  names_have_arg <- names(have_args)
 
  # remove unneeded args
  use_args       <- have_args[need_args]
  # remove NULL entries
  use_args       <- use_args[!is.na(names(use_args))]
  
  # warn about unused arguments
  superfluous_args <- 
    names_have_arg[!names_have_arg %in% 
                     c(names(use_args), "need_args","fun","method")]
  
  if (length(superfluous_args) > 0){
    superfluous_args <- paste(superfluous_args,collapse = ", ")
    message("following arguments not used: ",superfluous_args)
  }
  
  # pass in fitlered-down args as list
  do.call(fun, use_args)
}






