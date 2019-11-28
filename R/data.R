

# Code run on AvR's machine to produce package data file
# library(HMDHFDplus)
# library(tidyverse)
# # install.packages("usethis")
# #options(scipen=3)
# # an HMD life table as example data
# LT <- readHMD("C:\\hmd_statistics\\lt_female\\fltper_1x1\\CAN.fltper_1x1.txt") %>%
#   filter(Year==max(Year))
# 
# usethis::use_data(LT, overwrite = TRUE)

# documentation for data objects.

#' Year 2016 female Canadian single age lifetable from HMD
#'
#' @description Data taken from the file \code{fltper_1x1} using the default settings of the \code{HMDHFDplus} package.
#'
#' @format 
#' A data.frame with 111 rows and 11 columns
#' @source 
#' \url{www.mortality.org}
"LT"

# Note: use a bibtex citation in future for HMD.
