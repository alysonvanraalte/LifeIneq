#' @title ineq_mad
#' @description Calculate a lifetable column for the conditional mean absolute deviation in lifetable ages at death. This may be with respect to either conditional life expectancy or conditional median remaining lifespan. 
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in \code{MortalityLaws::MortalityLaw}). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package \code{ungroup} or with a penalized B-spline approach in package \code{MortalitySmooth}). 
#'   
#' @inheritParams ineq_var
#' @param center_type character either, \code{"ex"}, \code{"mean"} (same thing), or \code{"median"}
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
#' # A vector containing the conditional MAD in age at death 
#' M = ineq_mad(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The MAD in age at death from birth
#' M[1]
#' # The MAD in age at death conditional upon survival to age 10
#' M[11]
#' # MAD wrt condition median
#' ineq_mad(age = LT$Age, dx = LT$dx,
#'   lx = LT$lx, ex = LT$ex, ax = LT$ax, 
#'   center_type = "median")

ineq_mad <- function(age, 
                     dx, 
                     lx, 
                     ex, 
                     ax, 
                     center_type = c("ex","mean","median"),
                     check = TRUE){
  # dx <- dx / sum(dx)
  if (check){
    my_args <- as.list(environment())
    check_args(my_args)
  }
  
  
  # center on conditional mean or median
  center_type <- match.arg(center_type)
  if (center_type == "median"){
    ex <- ineq_quantile(age = age, 
                        lx = lx, 
                        quantile = .5)
  }
  

  
  axAge <- age + ax
  exAge <- ex + age

  rev(cumsum(rev(dx * abs(axAge - exAge)))) / lx
}

#' @title ineq_var
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
#' @param check logical. Shall we perform basic checks on input vectors? Default TRUE
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
#' V = ineq_var(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The variance in age at death from birth
#' V[1]
#' # The variance in age at death conditional upon survival to age 10
#' V[11]

ineq_var <- function(age, dx, lx, ex, ax, check = TRUE){

  # dx <- dx / sum(dx)
  if (check){
    my_args <- as.list(environment())
    check_args(my_args)
  }
  
  age0 <- age - age[1]
  n   <- length(age)
  out <- rep(NA, n)
  
  for (i in 1:n){
    axAge    <- age0[1:(n+1-i)] + ax[i:n]
    # dxi <- dx[i:n] / sum(dx[i:n])
    out[i] <- sum(dx[i:n] * (axAge - ex[i])^2) / lx[i]
  }
  out
  
}



#' @title ineq_sd
#' @description Calculate a lifetable column for the conditional standard deviation in lifetable ages at death 
#'
#' @inheritParams ineq_var
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional standard deviation in age at death 
#' S = ineq_sd(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The standard deviation in age at death from birth
#' S[1]
#' # The standard deviation in age at death conditional upon survival to age 10
#' S[11]


ineq_sd <- function(age, dx, lx, ex, ax, check = TRUE){
   # dx <- dx / sum(dx)
   V <- ineq_var(age = age, 
                 dx = dx, 
                 lx = lx, 
                 ex = ex, 
                 ax = ax, 
                 check = check)
   sqrt(V) 
}




#' @title ineq_cov
#' @description Calculate a lifetable column for the conditional coefficient of variation in lifetable ages at death 
#'
#' @inheritParams ineq_var
#' @param distribution_type character. Either \code{"aad"} (age at death) or \code{"rl"} (remaining life)
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional coefficient of variation in age at death 
#' CoV = ineq_cov(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The coefficient of variation in age at death from birth
#' CoV[1]
#' # The coefficient of variation in age at death conditional upon survival to age 10
#' CoV[11]


ineq_cov <- function(age, dx, lx, ex, ax, distribution_type = c("aad","rl"), check = TRUE){
  
  distribution_type <- match.arg(distribution_type)
  # dx <- dx / sum(dx)
  
  age_constant <- ifelse(distribution_type == "aad", age, age * 0)
  V <- ineq_var(age = age, 
                dx = dx, 
                lx = lx, 
                ex = ex, 
                ax = ax, 
                check = check)
  sqrt(V) / (ex + age_constant)
}



#' @title ineq_edag
#' @description Calculate a lifetable column for the conditional life disparity (\eqn{e^\dagger}) of a population.  
#'
#' @inheritParams ineq_var
#' @inherit ineq_var details
#' @inherit ineq_var seealso
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

ineq_edag <- function(age, dx, lx, ex, ax, check = TRUE){
  # dx <- dx / sum(dx)

  if (check){
    my_args <- as.list(environment())
    check_args(my_args)
  }
  
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
#' @description Calculate a lifetable column for the quantity *H*, generally referred to as either the lifetable entropy Keyfitz (1977) or the elasticity of life expectancy Leser (1955).
#'
#' @inheritParams ineq_var
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#' @references
#' \insertRef{keyfitz1977mortality}{LifeIneq}
#' \insertRef{leser1955variations}{LifeIneq}
#' 
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional H values
#' H = ineq_H(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The H from birth
#' H[1]
#' # The H conditional upon survival to age 10
#' H[11]


ineq_H <- function(age, dx, lx, ex, ax){
  # dx <- dx / sum(dx)
  ineq_edag(age, dx, lx, ex, ax) / ex
}




#' @title ineq_theil
#' @description Calculate a lifetable column for the conditional Theil index of inequality in survivorship
#'
#' @inheritParams ineq_var
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional Theil indices
#' Theil = ineq_theil(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The Theil index from birth
#' Theil[1]
#' # The Theil index conditional upon survival to age 10
#' Theil[11]


ineq_theil <- function(age, dx, lx, ex, ax, check = TRUE){
  # dx <- dx / sum(dx)
  if (check){
    my_args <- as.list(environment())
    check_args(my_args)
  }
  
  N     <- length(age)
  
  T1 <- rep(NA,N)
  for(i in 1:N){
    axAge <- age[1:(N+1-i)] + ax[i:N]
    T1[i] <- sum(dx[i:N] * (axAge / ex[i] * log( axAge / ex[i]))) / lx[i]
  }
  T1[T1 < 0] <- 0
  return(T1)
}





#' @title ineq_mld
#' @description Calculate a lifetable column for the conditional mean log deviation index of inequality in survivorship
#'
#' @inheritParams ineq_var
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional MLD indices
#' MLD = ineq_mld(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The MLD from birth
#' MLD[1]
#' # The MLD conditional upon survival to age 10
#' MLD[11]


ineq_mld <-  function(age, dx, lx, ex, ax, check = TRUE){
  # dx <- dx / sum(dx)
  if (check){
    my_args <- as.list(environment())
    check_args(my_args)
  }
  
  N     <- length(age)
  
  MLD <- rep(NA, N )
  for(i in 1: N ){
    # axAge <- age[1:(N+1-i)] + ax[i:N]
    axAge <- age[i:N] - age[i] + ax[i:N]
    dxi <- dx[i:N] / sum(dx[i:N])
    MLD[i] <- sum(
      dxi* (log (ex[i]/axAge))
    ) #/ lx[i]
    
  }
  MLD[MLD < 0] <- 0
  return(MLD)
}


#' @title ineq_gini
#' @description Calculate a lifetable column for the conditional Gini coefficient of inequality in survivorship
#'
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in \code{MortalityLaws::MortalityLaw}). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package \code{ungroup} or with a penalized B-spline approach in package \code{MortalitySmooth}). 
#' 
#'   The formula for calculating the Gini was taken from the Shkolnikov (2010) spreadsheet, and is a simplification of the formulas described in Shkolnikov (2003) and Hanada (1983).
#' 
#' @inheritParams ineq_var
#' @param distribution_type character. Either \code{"aad"} (age at death) or \code{"rl"} (remaining life)
#' @inherit ineq_var seealso
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
#' G = ineq_gini(age=LT$Age,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The Gini coefficient from birth
#' G[1]
#' # The Gini coefficient conditional upon survival to age 10
#' G[11]

# this is the Hanada version, which we didn't 
# retool to have both add and rle versions
# ineq_gini <- function(age, lx, ex, ax, check = TRUE){
#   if (check){
#     my_args <- as.list(environment())
#     check_args(my_args)
#   }
#   
#   nages <- length(age)
#   # vector of the length of the age interval
#   n <- c(diff(age),ax[nages])
#   
#   # squared survivorship
#   lx2 <- lx^2 / lx[1]^2
#   lx2plusn <- c(lx2[-1], 0)
#   
#   # the expression that will be integrated and the Gini 
#   intlx2 <- lx2plusn * n + ax * (lx2 - lx2plusn)
# 
#   G <- 1 - rev(cumsum(rev(intlx2))) / (ex * lx2)
#   G[G<0] <- 0
#   return(G)
# }
ineq_gini <- function(age, dx, ex, ax, distribution_type = c("aad","rl"), check = TRUE){
  
  distribution_type <- match.arg(distribution_type)

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
  g_out    <- rep(0, N)
  dx       <- dx / sum(dx)
  axage    <- age + ax
  denom    <- ex + age_constant
  ad       <- outer(axage, axage,"-") * lower.tri(diag(N),TRUE)
  pd       <- outer(dx, dx, "*")
  
  for (i in 1:N){
    g_out[i] <- sum(abs(ad) * pd) / (denom[i])
    ad       <- ad[-1, -1, drop = FALSE]
    pd       <- pd[-1, -1, drop = FALSE]
    pd       <- pd / sum(pd)
  }
  
  return(g_out)
}



#' @title ineq_aid
#' @description Calculate a lifetable column for the conditional absolute inter-individual difference in lifespan (AID)
#' 
#' @details All input vectors must be the same length. Also, we recommend using input data from a life table by single year of age with a highest age group of at least age 110. If your data have a lower upper age bound, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws'). If your data are abridged, consider first smoothing over age, and calculating a life table by single year of age (for instance by smoothing with a pclm model in package 'ungroup' or with a penalized B-spline approach in package 'MortalitySmooth'). 
#' 
#' The formula for calculating the AID was taken from the Shkolnikov 2010 spreadsheet, and is a simplification of the formula described in Shkolnikov 2003.
#' @references 
#' \insertRef{shkolnikov2010}{LifeIneq} 
#' \insertRef{shkolnikov2003}{LifeIneq} 
#' 
#' @inheritParams ineq_var
#' @inherit ineq_var seealso
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # A vector containing the conditional absolute inter-individual difference in lifespan
#' aid = ineq_aid(age=LT$Age,lx=LT$lx,ex=LT$ex,ax=LT$ax)
#' # The absolute inter-individual difference in lifespan from birth
#' aid[1]
#' # The absolute inter-individual difference in lifespan from age 10
#' aid[11]


ineq_aid <- function(age, lx, ex, ax, check = TRUE){
 
  aid <- ineq_gini(age = age, 
                   lx = lx, 
                   ex = ex, 
                   ax = ax, 
                   distribution_type = "rl",
                   check = check) * ex
  return(aid)
}




#' @title calculate a survivorship quantile 
#' @description Calculate quantiles of survivorship from a lifetable. Not vectorized: this function just calculates for the lowest age in the vector given. 
#' @details Any of the quantile-based functions in 
#'
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#'
#' @inheritParams ineq_var
#' @param quantile numeric. value between zero and one indicating the quantile desired.
#'
#' @export
#' @importFrom stats splinefun
#' @examples 
#'
#' data(LT)
#' # The median age at survival
#' LTmedian <- ineq_quantile_lower(age=LT$Age,lx=LT$lx,quantile=0.5)
#' LTmedian
#' # The age reached by 90% of the life table cohort (i.e. the top 10%)
#' ineq_quantile_lower(age=LT$Age,lx=LT$lx,quantile=0.1) 
#' # The difference between the bottom 10 and top 10 percent of survival age
#' ineq_quantile_lower(age=LT$Age,lx=LT$lx,quantile=0.1) -
#' ineq_quantile_lower(age=LT$Age,lx=LT$lx,quantile=0.9)

ineq_quantile_lower <- function(age, lx, quantile = .5){
  stopifnot(length(age) == length(lx))
  lx   <- lx / lx[1]
  # make this monotonic
  splinefun(age~lx)(quantile)
}

#' @title calculate a conditional survivorship quantile 
#' @description Calculate quantiles of survivorship from a lifetable, returns full lifetable column.
#'
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#'
#' @inheritParams ineq_var
#' @param quantile numeric. value between zero and one indicating the quantile desired.
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


ineq_quantile <- function(age, lx, quantile = .5){
  stopifnot(length(age) == length(lx))
  n    <- length(age)
  # make sure it closes out
  lx   <- c(lx, 0)
  a    <- c(age, max(age) + 1)
  
  
  qs   <- rep(0,n)
  for (i in 1:n){
    qs[i] <- ineq_quantile_lower(age = a[i:(n+1)],
                                 lx = lx[i:(n+1)],
                                 quantile = quantile)
  }
  qs - age
}


#' @title ineq_iqr
#' @description Calculate the interquartile range survivorship age from a lifetable. Other quantile ranges can also be calculated.
#'
#' @inherit ineq_var details
#' @inherit ineq_var seealso
#' 
#' @inheritParams ineq_var
#' @param upper numeric. upper survival quantile, default .75
#' @param lower numeric. lower survival quantile, defauly .25
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The iqr range of survival age
#' iqr <- ineq_iqr(age=LT$Age,lx=LT$lx)
#' iqr


ineq_iqr <- function(age, lx, upper = .75, lower = .25){
  stopifnot(length(age) == length(lx))
  q1 <- ineq_quantile_lower(age = age, lx = lx, quantile = lower) 
  q3 <- ineq_quantile_lower(age = age, lx = lx, quantile = upper)
  q1 - q3
}



#' @title ineq_cp
#' @description Calculate Kannisto's C-measures from a lifetable
#'
#' @details The age and lx vectors must be the same length. This function estimates the shortest distance between two ages containing p percent of the life table cohort's death. The mechanics behind the function are to fit a cubic spline through the survival curve to estimate surivorship between age intervals. If your data have an upper age bound lower than 110, consider extrapolation methods, for instance a parametric Kannisto model (implemented in package 'MortalityLaws').
#' 
#'  This function is currently very slow because the author is a better demographer than she is programmer. For ideas of how to improve this function, please let email vanraalte@@demogr.mpg.de.
#'  
#'   The concept behind Kannisto's C-measures is found in Kannisto (2000) 
#' 
#'
#' @inheritParams ineq_var
#' @param p numeric. What proportion of the life table cohort do you want captured in the C measure? The default is .5
#' @importFrom stats splinefun
#' @importFrom stats optimize
#'
#'@references 
#' \insertRef{kannisto2000}{LifeIneq}
#'
#' @export
#' @examples 
#'
#' data(LT)
#' # The shortest age range containing half of the deaths
#' (C50 <- ineq_cp(age=LT$Age,lx=LT$lx,p=.5))

ineq_cp <- function(age, lx, p = .5){

  stopifnot(length(age) == length(lx))
  stopifnot(p <= 1)
  
  # span minimizer function
  ineq_cp_min <- function(par = 60, fun, funinv, p = .5){
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
  age1 <- optimize(ineq_cp_min, 
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
#' @description Choose from variance \code{var}, standard deviation \code{sd}, coefficient of variation \code{cov}, interquartile range \code{iqr}, gini, absolute interindividual difference (absolute gini) \code{aid}, edagger, Kannisto's compression measure \code{cp}, Leser-Keyfitz entropy \code{H}, theil, mean log deviation \code{mld}, mean absolute deviation \code{mad} (wrt mean or median).
#' @param age numeric. vector of lower age bounds.
#' @param dx numeric. vector of the lifetable death distribution.
#' @param lx numeric. vector of the lifetable survivorship.
#' @param ex numeric. vector of remaining life expectancy.
#' @param ax numeric. vector of the average time spent in the age
#' @param method one of \code{c("var","sd","cov","iqr","aid","gini","mld","edag","cp","theil","H","mad")}
#' @param check logical. Shall we perform basic checks on input vectors? Default TRUE
#' @param ... other optional arguments used by particular methods.
#' @importFrom Rdpack reprompt
#' @export

ineq <- function(age, 
                 dx, 
                 lx, 
                 ex, 
                 ax, 
                 method = c("var","sd","cov","iqr","aid",
                            "gini","mld","edag","cp","theil","H","mad"), 
                 check = TRUE, 
                 ...){
  
  # make sure just one
  method         <- match.arg(method)
  # fun is now the function we need
  fun_name       <- paste0("ineq_", method)
  fun            <- match.fun(paste0("ineq_", method))
  
  # what do we need and what do we have?
  # mget(names(formals(sys.function(sys.parent()))), sys.frame(sys.nframe() - 1L))
  given_args <- mget(names(formals()),sys.frame(sys.nframe()))
  given_args[["..."]] <- NULL
  
  
  my_args     <- lapply(as.list(match.call())[-1], eval)
  # cl        <- sys.call(0)
  # f         <- get(as.character(cl[[1]]), mode="function", sys.frame(-1))
  # cl        <- match.call(definition=f, call=cl)
  # my_args <- as.list(cl)[-1]
  have_args <- c(as.list(environment()), list(...)) 
  # remove objects created in this function prior to this line..
  have_args <- have_args[!names(have_args)
                         %in%
                          c("need_args","fun","have_args","fun_name",
                            "given_args","my_args","cl","f")]
  have_args <- c(have_args,given_args[!names(given_args) %in% names(have_args)])
  # names_have_arg <- names(have_args)

  # remove unneeded args
  need_args      <- names(formals(fun))
  use_args       <- have_args[names(have_args) %in% need_args]
  # remove NULL entries
  use_args       <- use_args[!is.na(use_args)]
  use_args       <- use_args[!lapply(use_args,is.symbol) |> unlist()]
  # potentially warn about unused arguments
  superfluous_args <- 
    names(have_args[!names(have_args) %in% 
                     c(need_args, "need_args","fun","method")])
  
  # throw error if arg missing
  if (!all(need_args %in% names(use_args))){
    missing_arg <- need_args[!need_args %in% names(use_args)] |>
      paste(collapse = ", ")
    stop(paste(method,"method requires missing argument(s)",missing_arg))
  }
 
  
  if (length(superfluous_args) > 0 & check){
    superfluous_args <- paste(superfluous_args,collapse = ", ")
    message("following arguments not used: ",superfluous_args)
  }
  
  # pass in filtered-down args as list
  do.call(fun, as.list(use_args))
}



