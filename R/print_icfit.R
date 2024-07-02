#' Print method for an object of class `icfit`
#'
#' @method print icfit
#  Solution found on https://stackoverflow.com/questions/23724815/roxygen2-issue-with-exporting-print-method
#' @export
#' 
#' @param x The object of class 'icfit' to be printed
#' @param digits Number of digits to be printed
#' @param alpha Alpha level to be used of confidence interval ((1-alpha) * 100 percent)
#' @param \dots Further arguments to print
#'
#' @return No return value
#'
print.icfit <- function(x, 
                        digits = max(1L, getOption("digits") - 3L),
                        alpha=0.05, ...) {

  if (!inherits(x, "icfit"))
    stop("'x' must be an 'icfit' object")

  summary(x, lvl=0, digits=digits, alpha=alpha, ...)
  return(invisible())
  
}

