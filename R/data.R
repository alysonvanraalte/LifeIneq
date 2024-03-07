

# Code run on AvR's machine to produce package data file
# library(HMDHFDplus)
# library(tidyverse)
# # install.packages("usethis")
# #options(scipen=3)
# # an HMD life table as example data
# LT <- readHMDweb("CAN","fltper_1x1",username =Sys.getenv("us"),password=Sys.getenv("pw")) |>
# filter(Year == 2016)
# LTm <- readHMDweb("CAN","mltper_1x1",
#                   username = Sys.getenv("us"),
#                   password = Sys.getenv("pw")) |>
#    filter(Year == 2016)
# usethis::use_data(LTm, overwrite = TRUE)
# 
# usethis::use_data(LT, overwrite = TRUE)

# documentation for data objects.

#' Year 2016 female Canadian single age lifetable from HMD
#'
#' @description Data taken from the file \code{fltper_1x1} using the default settings of the \code{HMDHFDplus} package.
#' @details Data downloaded 9 Feb 2023.
#'
#' @format 
#' A data.frame with 111 rows and 11 columns
#' @source 
#' \url{https://mortality.org}
"LT"

#' Year 2016 male Canadian single age lifetable from HMD
#'
#' @description Data taken from the file \code{mltper_1x1} using the default settings of the \code{HMDHFDplus} package.
#' @details Data downloaded 6 march 2024.
#'
#' @format 
#' A data.frame with 111 rows and 11 columns
#' @source 
#' \url{https://mortality.org}
"LTm"


# Note: use a bibtex citation in future for HMD.
# LTabr <- readHMDweb("CAN","fltper_5x1",username =Sys.getenv("us"),password=Sys.getenv("pw")) |>
# filter(Year == 2016)
# usethis::use_data(LTabr, overwrite = TRUE)

#' Year 2016 female Canadian abridged age lifetable from HMD
#'
#' @description Data taken from the file \code{fltper_1x1} using the default settings of the \code{HMDHFDplus} package.
#' @details Data downloaded 9 Feb 2023.
#' 
#' @format 
#' A data.frame with 24 rows and 11 columns
#' @source 
#' \url{https://mortality.org}
"LTabr"
