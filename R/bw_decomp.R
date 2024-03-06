#' @title between-within decomposition of lifespan inequality measures

#' @description 
#' Partition a lifespan inequality index into additive components 
#' of between-group inequality and within-group inequality. Presently 
#' implemented for Theil's index, e-dagger, variance, mean log deviation, 
#' the gini coeficient, mean absolute deviatation, the absolute inter-individual 
#' difference, and life table entropy (H). 
#' 
#' @param age numeric vector of lower age bounds.
#' @param dx numeric matrix of the lifetable death distribution with age in rows and subgroups in columns.
#' @param lx numeric matrix of the lifetable survivorship with age in rows and subgroups in columns.
#' @param ex numeric matrix of remaining life expectancy with age in rows and subgroups in columns.
#' @param ax numeric matrix of the average time spent in the age interval of those dying within the interval with age in rows and subgroups in columns.
#' @param distribution_type \code{"aad"} for age-at death, \code{"rl"} for remaining life
#' @param prop numeric vector of starting fractions for each of the subgroups.
#' @param method character one of \code{"theil", "edag","var","mld","gini","mad","aid","H"}
#' @import LifeIneq
#' @export
#' @examples
#' # example code
#' 
bw_decomp <- function(age, 
                      ax, 
                      dx, 
                      lx, 
                      ex,  
                      distribution_type = "aad",
                      prop = rep(1/ncol(dx), ncol(dx)),
                      method = c("theil", "edag","var","mld","gini","mad","aid","H")){
  
  
  # Turn data frames into matrices
  if(is.data.frame(age)) age <- as.numeric(age[,1])
  if(is.data.frame(ax)) ax <- as.matrix(ax)
  if(is.data.frame(dx)) dx <- as.matrix(dx)
  if(is.data.frame(lx)) lx <- as.matrix(lx)
  if(is.data.frame(ex)) ex <- as.matrix(ex)
  
  # check dimensions: number of groups
  K <- length(prop)
  stopifnot(all.equal(ncol(ax),
                      ncol(dx), 
                      ncol(lx), 
                      ncol(ex), 
                      K))
  
  # check dimensions: age groups
  N <- length(age)
  stopifnot(all.equal(nrow(ax),
                      nrow(dx), 
                      nrow(lx), 
                      nrow(ex),
                      N))
  
  # validate method selection
  method <- match.arg(method)
  
  # 1) standardize inputs
  prop <- prop / sum(prop)
  lx   <- lx %*% diag(1 / lx[1, ])
  dx   <- dx %*% diag(1 / colSums(dx))
  
  # 2) get pop avgs and other goods
  
  # weighted dx and lx
  pdxm <- dx %*% diag(prop)
  plxm <- lx %*% diag(prop)
  
  # pop avg dx and lx
  plx  <- rowSums(plxm)
  pdx  <- rowSums(pdxm)
  
  # need rows to sum to one to weight ax and ex
  plxc <- diag(1 / plx) %*% plxm
  pdxc <- diag(1 / pdx) %*% pdxm
  
  # pop avg ax is dx-weighted ax
  pax  <- rowSums(ax * pdxc)
  # pop avg ex is lx-weighted ex
  pex  <- rowSums(ex * plxc)
  
  # If method %in% c("H","aid") then use edag or gini, rescale later
  method_lower <- method
  if(method=="H")   method_lower <- "edag"
  if(method=="aid") method_lower <- "gini"
  
  # calculate total inequality index
  args_i <- list(age = age, 
                 dx = pdx,
                 lx = plx,
                 ax = pax,
                 ex = pex, 
                 distribution_type = distribution_type,
                 check = FALSE,
                 method = method_lower)
  tot <- suppressMessages(do.call("ineq", args = args_i, quote = TRUE)[1])
  
  
  # again for each of the k subgroups 
  indices <- rep(0, K)
  for (k in 1:K){
    args_i <- list(age = age, 
                   dx = dx[, k],
                   lx = lx[, k],
                   ax = ax[, k],
                   ex = ex[, k], 
                   distribution_type = distribution_type,
                   method = method_lower)
    indices[k] <- suppressMessages(do.call("ineq", args = args_i, quote = TRUE)[1])
    
  }
  
  # within weighting depends on the measure
  if (method_lower %in% c("edag","var","mld","mad")){ 
    weights <- prop
  }
  
  if(method_lower=="aid") {
    weights <- prop * lx[1,]^2 / plx[1]^2
  }
  
  if (method_lower =="theil" ){
    if(distribution_type=="aad") weights <- prop * (ex[1, ]+min(age)) / (pex[1]+min(age))
    if(distribution_type=="rl") weights <- prop * ex[1, ] / pex[1]
  }
  
  if (method_lower =="gini"){
    if(distribution_type=="aad") weights <- prop * (lx[1,]^2 * (ex[1, ]+min(age))) / (plx[1]^2*(pex[1]+min(age)))
    if(distribution_type=="rl") weights <- prop * (lx[1,]^2 * ex[1, ]) / (plx[1]^2*pex[1])
  }
  
  # combine
  W <- sum(weights * indices)
  
  # Adjust if H or aid
  if(method%in%c("H","aid")) {
    rescale <- sum(ex[1,]*prop)
    if(method=="H") rescale <- 1/rescale
    W <- W*rescale
    tot <- tot*rescale
  }
  
  # between part as residual
  B <- tot - W
  
  # gather minimal goods to return
  out <- list(method = method,
              group_ind = indices,
              tot = tot,
              B = B,
              W = W,
              fB = B / tot,
              fW = W / tot)
  
  # output
  return(out)
}
