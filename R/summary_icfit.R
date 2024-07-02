#' Summary method for an object of class `icfit`
#'
#' @param object Object of class 'icfit'
#' @param lvl Describes the level of output
#' @param digits Number of digits to be printed
#' @param alpha Alpha level to be used of confidence interval ((1-alpha) * 100 percent)
#' @param \dots Other arguments to summary
#'
#' @return None (invisible \code{NULL})
#'
#' @examples
#' \donttest{
#' icf <- icfit(Surv(left, right, type='interval2') ~ period + gender + age, data=drugusers)
#' summary(icf)
#' summary(icf, lvl=0) # same as print(icf)
#' summary(icf, lvl=1) # extra information on iterations and computation time
#' }
#' 
#' @method summary icfit
#' @export
summary.icfit <- function(object, lvl=1,
                          digits = max(1L, getOption("digits") - 3L),
                          alpha=0.05, ...) {
  
  if (!inherits(object, "icfit"))
    stop("'object' must be a 'icfit' object")
  
  # Call
  cat("Object of class 'icfit'\n")
  if (!is.null(cl <- object$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  # Intervals
  tmin <- 0
  bins <- object$bins
  cat("Bins summary: tmin = ", tmin, ", tmax = ", bins$tmax,
      ", number of bins = ", bins$nt,
      ", bin width = ", bins$dt, "\n", sep="")
  
  # Splines
  splines <- object$control
  cat("P-splines summary: number of segments = ", splines$nseg,
      ", degree = ", splines$bdeg,
      ", penalty order = ", splines$pord, "\n\n", sep="")
  
  # Coefficient table
  cat("Parameter estimates:\n")
  X <- as.matrix(object$data$X, drop = FALSE)
  B <- object$basis$B
  
  cbx <- object$cbx
  n <- nrow(X)
  nt <- nrow(B)
  nb <- ncol(B)
  nx <- ncol(X)
  
  if (nx>2) {
    # Mpen <- object$Mpen
    # Hinv <- solve(object$Mpen)
    Hinv <- solve(object$Itot)
    Covx <- Hinv[nb + (2:nx), nb + (2:nx), drop = FALSE]
    sex <- sqrt(diag(Covx))
    coef <- cbx[nb + (2:nx)]
    z <- coef / sex
    HR <- exp(coef)
    lower <- exp(coef - qnorm(1-alpha/2) * sex)
    upper <- exp(coef + qnorm(1-alpha/2) * sex)
    pvals <- 2 * (1 - pnorm(abs(z)))
    cmat <- cbind(coef, sex, HR, lower, upper, pvals)
    colnames(cmat) <- c("coef", "SE", "HR", "lower", "upper", "pvalue")
    printCoefmat(cmat, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    if (alpha != 0.05) cat("Note: ", 100*(1-alpha), "% confidence interval\n")
  }
  
  # Likelihood, lambda, effective dimension
  cat("\nEffective baseline dimension: ", object$ed,
      ", log-likelihood: ", object$ll,
      ", AIC: ", object$aic, "\n", sep="")
  cat("Smoothness parameter lambda: ", object$lambda, "\n", sep="")
  cat("Number of iterations: ", object$nit, "\n", sep="")
  
  # n, number of events
  cat("n = ", n, ", number of events = ", sum(object$dead), "\n", sep="")
  
  if (lvl >= 1) {
    # Lambda tolerance
    cat("\nNumber of iterations: ", object$nit1, " before updating lambda, ",
        object$nit, " in total\n", sep="")
    cat("Lambda tolerance (used to decide when to switch to updating lambda): ",
        object$tollam)
    cat("\nComputation time:\n")
    print(object$elapsed)
  }
  
  return(invisible())
}