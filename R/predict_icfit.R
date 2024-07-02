#' Predict method for an object of class `icfit`
#'
#' @importFrom stats as.formula
#'
#' @param object The object of class 'icfit' for which a prediction is to be made
#' @param newdata A data frame containing covariate information for a new subject
#' @param nstep Number of time steps used for calculating cumulative hazards (default is 500)
#' @param alpha The alpha level for the (1-alpha)*100 percent confidence interval
#' @param ... Any other arguments
#'
#' @return An object of class `predict.icfit`, which is a data frame with time points and hazard,
#' cumulative hazard and survival at those time points, along with standard errors
#' and pointwise lower and upper confidence bounds, or a list of such data frames
#' for each subject represented in `newdata`
#'
#' @examples
#' \donttest{
#' icf <- icfit(Surv(left, right, type='interval2') ~ period + gender + age, data=drugusers)
#' pred_icf <- predict(icf)
#' head(pred_icf)
#' ndata <- drugusers[1:4, ]
#' pred_nd_icf <- predict(icf, newdata=ndata)
#' lapply(pred_nd_icf, head)
#' }
#' 
#' @export
#' 
predict.icfit <- function(object, newdata, nstep = 500, alpha = 0.05, ...) {
  if (!inherits(object, "icfit"))
    stop("'object' must be an 'icfit' object")
  
  cov_diag = function(Cov, U) {
    # Compute diagonal of U %*% solve(Cov) %*% t(U) efficiently
    V = U %*% Cov
    d = rowSums(U * V)
  }
  
  # Extract variables
  B <- object$basis$B
  dt <- object$bins$dt
  cbx <- object$cbx
  nt <- nrow(B)
  nb <- ncol(B)
  nb1 <- nb + 1

  X <- as.matrix(object$data$X, drop = FALSE)
  n <- nrow(X)
  nx <- ncol(X)

  # Compute standard errors of curve and confidence intervals 
  Hinv <- solve(object$Itot)
  sd <- 1  # Poisson assumption (may change later)
  cb <- cbx[1:nb1]
  Covcb <- sd * Hinv[1:nb1, 1:nb1]
  z <- qnorm(1-alpha/2)

  steptime <- seq(0, object$bins$tmax, length = nstep+1) # running time for cumulative hazard
  BB <- bbase(x = steptime / dt, xl = 0, xr = nt, nseg = object$control$nseg, deg = object$control$bdeg)
  BB <- cbind(BB, 1) # Extra column of ones needed for proper standard errors
  
  # Cumulative (baseline) hazard, computed using midpoints of intervals
  # (hence recalculation of steptime and BB)
  steptimemid <- (steptime[-1] + steptime[-(nstep+1)]) / 2 # midpoints
  dtstep <- object$bins$tmax / nstep
  # Recompute baseline hazard at midpoints for cumulative hazards
  BBmid <- bbase(x = steptimemid / dt, xl = 0, xr = nt, nseg = object$control$nseg, deg = object$control$bdeg)
  BBmid <- cbind(BBmid, 1) # Extra column of ones needed for proper standard errors
  
  #
  # No newdata
  #
  if (missing(newdata)) {
    # Baseline hazard
    se <- sqrt(cov_diag(Covcb, BB)) # = sqrt(diag(BB %*% Covcb %*% t(BB)))
    eta <- BB %*% cb
    h0 <- exp(eta) / dt
    width <- z * se
    lowerh0 <- exp(eta - width) / dt
    upperh0 <- exp(eta + width) / dt
    seh0 <- se * h0 # Delta method
    
    # Baseline cumulative hazard
    eta <- BBmid %*% cb
    h0mid <- exp(eta) / dt
    Coveta <- BBmid %*% Covcb %*% t(BBmid) # covariance matrix of vector log(h0) at midpoints
    tmp <- h0mid %*% t(h0mid)
    Covh0mid <- tmp * Coveta # covariance matrix of vector h0 at midpoints
    cummat <- diag(dtstep, nstep)
    cummat[lower.tri(cummat)] <- dtstep
    H0 <- c(cummat %*% h0mid) # or cumsum(h0mid) * dtstep
    seH0 <- sqrt(cov_diag(Covh0mid, cummat)) # = sqrt(diag(cummat %*% Covh0mid %*% t(cummat)))

    lowerH0 <- H0 - z * seH0
    lowerH0[lowerH0 < 0] <- 0
    upperH0 <- H0 + z * seH0
    # Add time=0 to H0 and seH0
    H0 <- c(0, H0)
    seH0 <- c(0, seH0)
    lowerH0 <- c(0, lowerH0)
    upperH0 <- c(0, upperH0)
    
    # Baseline survival
    S0 <- exp(-H0)
    seS0 <- S0 * seH0
    lowerS0 <- exp(-upperH0)
    upperS0 <- exp(-lowerH0)
    
    res <- data.frame(time = seq(0, object$bins$tmax, length = nstep+1),
                      haz = h0, sehaz = seh0, lowerhaz = lowerh0, upperhaz = upperh0,
                      Haz = H0, seHaz = seH0, lowerHaz = lowerH0, upperHaz = upperH0,
                      Surv = S0, seSurv = seS0, lowerSurv = lowerS0, upperSurv = upperS0)
    class(res) <- c("predict.icfit", "data.frame")
  } else { # with newdata
    Covcb <- sd * Hinv
    npats <- nrow(newdata)

    # Get design matrix, using make_xy_pred
    # Response part is stripped from original formula
    frml <- stats::as.formula(paste("~", strsplit(as.character(object$formula), "~")[[3]]))
    reg_items <- make_xy_pred(frml, newdata) # adapted version of make_xy for predictions used here
    mm <- reg_items$x
    p <- ncol(mm)
    res <- list()
    
    for (i in 1:npats) {

      # Hazard
      XBB <- matrix(as.numeric(cbind(BB, matrix(mm[i, , drop = FALSE], 
                                                nstep+1, p, byrow = TRUE))),
                    nstep+1, nb1 + p)
      Coveta <- XBB %*% Covcb %*% t(XBB) # covariance matrix of vector log(h0)
      se <- sqrt(diag(Coveta))
      eta <- XBB %*% cbx
      h <- exp(eta) / dt
      width = z * se
      lowerh = exp(eta - width) / dt
      upperh = exp(eta + width) / dt
      seh <- se * h # Delta method
      
      # Cumulative hazard
      XBBmid <- matrix(as.numeric(cbind(BBmid, matrix(mm[i, , drop = FALSE],
                                                      nstep, p, byrow = TRUE))),
                       nstep, nb1 + p)
      eta <- XBBmid %*% cbx
      hmid <- exp(eta) / dt
      Coveta <- XBBmid %*% Covcb %*% t(XBBmid) # covariance matrix of vector log(h) at midpoints
      tmp <- hmid %*% t(hmid)
      Covhmid <- tmp * Coveta # covariance matrix of vector h
      cummat <- diag(dtstep, nstep)
      cummat[lower.tri(cummat)] <- dtstep
      H <- cummat %*% hmid # or cumsum(hmid) * dtstep
      CovH <- cummat %*% Covhmid %*% t(cummat)
      seH <- sqrt(diag(CovH))
      lowerH <- H - z * seH
      lowerH[lowerH < 0] <- 0
      upperH <- H + z * seH
      # Add time=0 to H and seH
      H <- c(0, H)
      seH <- c(0, seH)
      lowerH <- c(0, lowerH)
      upperH <- c(0, upperH)
      
      # Survival
      S <- exp(-H)
      seS <- S * seH
      lowerS <- exp(-upperH)
      upperS <- exp(-lowerH)
      
      res[[i]] <- data.frame(time = seq(0, object$bins$tmax, length = nstep+1),
                        haz = h, sehaz = seh, lowerhaz = lowerh, upperhaz = upperh,
                        Haz = H, seHaz = seH, lowerHaz = lowerH, upperHaz = upperH,
                        Surv = S, seSurv = seS, lowerSurv = lowerS, upperSurv = upperS)
      class(res[[i]]) <- c("predict.icfit", "data.frame")
    }
    
    class(res) <- c("predict.icfit", "list")
  }
  attr(res, "alpha") <- alpha
  return(res)
}
