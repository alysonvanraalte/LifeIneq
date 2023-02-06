
# here we should see what happens if 0s or negatives are in inputs

check_positive <- function(x){
  stopifnot(all(x > 0))
}
check_nas <- function(x){
  stopifnot(all(!is.na(x)))
}
check_lx <- function(lx){
  stopifnot(all(diff(lx) <= 0))
  check_positive(x=lx)
  check_nas(x=lx)
}

check_age <- function(age){
  stopifnot(all(diff(age) > 0))
  stopifnot(all(age >= 0))
  stopifnot(all(diff(age) > 0))
  check_nas(x=age)
}

check_ax <- function(ax,age){
  check_positive(x=ax)
  stopifnot(all(ax <= c(diff(age),30)))
  stopifnot(all(diff(ax + age) > 0))
  check_nas(x=ax)
}

check_dx <- function(dx){
  check_positive(x=dx)
  check_nas(x=dx)
}

check_ex <- function(ex, age){
  check_positive(x=ex)
  stopifnot(diff(age + ex)>= 0)
  check_nas(x=ex)
}

check_vec_arg <- function(x,item = c("age","lx","ax","dx","ex"),...){
  switch(item,
         age = check_age(age),
         lx = check_lx(lx),
         ax = check_ax(ax,age),
         dx = check_dx(dx),
         ex = check_ex(ex,age))
}
check_args <- function(arg_list){
  L <- lapply(arg_list, length) |> unlist()
  arg_list[L == 0] <- NULL
  age_lengths <- c(length(arg_list$age),
              length(arg_list$lx),
              length(arg_list$ex),
              length(arg_list$ax))
  # dx not in some functions...
  if (any(names(arg_list) == "dx")){
    check_vec_arg(x=arg_list$dx, "dx")
    age_lengths <- c(age_lengths, length(arg_list$dx))
  }
  lengths_match <- diff(range(age_lengths)) == 0
  if (!lengths_match){
    stop("vector argument lengths must match")
  }
  
  check_vec_arg(x=arg_list$age, "age")
  check_vec_arg(x=arg_list$lx, "lx")
  check_vec_arg(x=arg_list$ex, "ex", age = arg_list$age)
  check_vec_arg(x=arg_list$ax, "ax", age = arg_list$age)
}

# to remove CMD check warning
globalVariables(names = c("age","ax","dx","lx", "ex","check"))
