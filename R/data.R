#' Colorado daily precipitation accumulations
#'
#' Three objects: 1) \code{COprcp}, a 404,326-row
#' data frame with columns \code{date}, \code{prcp} and \code{meta_row};
#' 2) \code{COprcp_meta}, a 64-row data frame, with meta data for 64 stations.
#' 3) \code{COelev}, a list of elevation for the domain at 0.02 x 0.02
#' degree resolution. Precipitation amounts are only given for April
#' to October in the years 1990 - 2019. The domain has a longitude range
#' of [-106, -104] and a latitude range [37, 41]. These choices reflect
#' the analysis of Cooley et al. (2007).
#'
#' @format A data frame with 2383452 rows and 8 variables
#' 
#' The variables are as follows:
#'
#' \describe{
#'   \item{date}{date of observation}
#'   \item{prcp}{daily rainfall accumulation in mm}
#'   \item{meta_row}{an identifier for the row in COprcp_meta; see `Examples'}
#'   \item{lon}{longitude of station}
#'   \item{lat}{latitude of station}
#'   \item{elev}{elevation of station in metres}
#'   \item{id}{GHCDN identifier}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name COprcp
#' @usage data(COprcp) # loads all three objects
#'
#' @references
#'
#' Cooley, D., Nychka, D., & Naveau, P. (2007). Bayesian spatial modeling of 
#' extreme precipitation return levels. Journal of the American Statistical 
#' Association, 102(479), 824-840.
#'
#' @examples
#' 
#' library(evgam)
#' data(COprcp)
#'
#' brks <- pretty(COelev$z, 50)
#' image(COelev, breaks=brks, col=rev(heat.colors(length(brks[-1]))))
#' colplot(COprcp_meta$lon, COprcp_meta$lat, COprcp_meta$elev, breaks=brks, add=TRUE)
#'
NULL

#' @rdname COprcp
#' @name COprcp_meta
NULL

#' @rdname COprcp
#' @name COelev
NULL

#' Fort Collins, Colorado, US daily max. temperatures
#'
#' Daily maximum temperatures at Fort Collins, Colorado, US 
#' from 1st January 1970 to 31st December 2019
#'
#' @format A data frame with 18156 rows and 2 variables
#'
#' The variables are as follows:
#'
#' \describe{
#'   \item{date}{date of observation}
#'   \item{tmax}{daily maximum temperature in degrees Celcius}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name FCtmax
#' @usage data(FCtmax)
#'
#' @examples
#' 
#' library(evgam)
#' data(FCtmax)
NULL
