#' Interval-censored drug users data
#'
#' @docType data
#'
#' @usage data(drugusers)
#'
#' @format Data from a cohort of 940 injecting drug users attending a hospital
#' detoxification unit in Barcelona, Spain. Time is months between initiation
#' of intravenous drug use and HIV seroconversion.
#' A dataframe with five columns:
#'
#' \describe{
#'   \item{\code{left}}{Last negative HIV test (0 if first HIV test was positive)}
#'   \item{\code{right}}{First positive HIV test (\code{Inf} if last HIV test was negative)}
#'   \item{\code{period}}{Period of initiation of drug use,
#'   factor with levels "1972-1980", "1981-1985", "1986-1991", "1992-1997"}
#'   \item{\code{gender}}{Gender, factor with levels "male" and "female"}
#'   \item{\code{age}}{Age at initiation of drug use (years)}
#' }
#'
#' @keywords datasets
#'
#' @references Gomez G, Calle ML, Egea JM & Muga R (2000). Risk of HIV infection as a function of the duration
#' of intravenous drug use: a non-parametric Bayesian approach. Stat Med; 19:2641–2656.
#' 
#' 
"drugusers"
