#' Summary method for an object of class `predict.icfit`
#'
#' @param object Object of class 'predict.icfit'
#' @param times The time points at which to summarize the predicted hazards, cumulative hazards
#' and survival probabilities, with associated standard errors and confidence intervals
#' @param \dots Other arguments to plot
#'
#' @return A data frame (if object was a data frame) or a list of data frames (if object was
#' a list of data frames) with hazards etc linearly interpolated between the time points used
#' in the predict function
#'
#' @examples
#' \donttest{
#' icf <- icfit(Surv(left, right, type='interval2') ~ period + gender + age,
#'       data=drugusers)
#' pred_icf <- predict(icf)
#' summary(pred_icf, times=c(0, 30, 183, 365))
#' }
#' 
#' @method summary predict.icfit
#' @export
summary.predict.icfit <- function(object, times, ...) {
  
  if (!inherits(object, "predict.icfit"))
    stop("'object' must be a 'predict.icfit' object")
  
  times <- sort(times)
  maxt <- 0
  if (inherits(object, "data.frame")) maxt <- max(object$time)
  else if (inherits(object, "list")) maxt <- max(object[[1]]$time)
  if (max(times) > maxt) {
    warning("Max of 'times' larger than time domain of 'object', replacing by max time domain")
    times <- times[times < maxt]
    times <- c(times, maxt)
  }
  
  approxmat <- function(x, ymat, xout) {
    nms <- names(ymat)
    nc <- ncol(ymat)
    res <- matrix(NA, length(xout), nc+1)
    res[, 1] <- xout
    for (k in 1:nc) res[, k+1] <- stats::approx(x, ymat[, k], xout)$y
    res <- as.data.frame(res)
    names(res) <- c("time", nms)
    res
  }

  if (inherits(object, "data.frame")) {
    # This is when no newdata object has been used
    res <- approxmat(object$time, object[, -1], times)
  } else if (inherits(object, "list")) {
    # This is when a newdata object has been used
    res <- lapply(object, function(x) approxmat(x$time, x[, -1], times))
  }
  
  return(res)
}
