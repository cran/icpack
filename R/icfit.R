#' Fit a proportional hazards model with baseline hazard modeled by P-splines
#'
#' @export
#' @importFrom methods is
#'
#' @param formula A formula object with response of the left of a ~ operator
#' and covariate terms on the right. The response must be a survival object
#' as returned by the `Surv` function, with type either right', 'counting' or
#' 'interval2'
#' @param data A data frame in which to interpret the variable names in the
#' `formula`
#' @param entry When appropriate, a vector of entry (left truncation) times,
#' or a string indicating the column name in `data` containing entry times;
#' only used if Surv object is of type 'interval2'
#' @param lambda Starting value of penalty tuning parameter
#' @param nt The number of time bins
#' @param tmax The end of time domain (default 1.01 times largest observation)
#' @param nseg The number of B-spline segments
#' @param bdeg The degree of the B-splines
#' @param pord The order of the differences used in the penalty
#' @param nit Maximum number of iterations (integer)
#' @param tol Tolerance for final fit
#' @param tollam Tolerance for switching to lambda update
#' @param kappa Ridge parameter (number)
#' @param update_lambda Automatic update of lambda (Boolean)
#' @param ic_update Update risk and event probabilities (Boolean)
#' @param monitor Monitor convergence (Boolean)
#'
#' @return An object of class `icfit`
#' 
#' @examples
#' # Fit proportional hazards model to interval-censored data
#' icfit(Surv(left, right, type='interval2') ~ period + gender + age,
#'       data=drugusers)
#' # Fit proportional hazards model to right-censored data
#' icfit(Surv(time, d) ~ Diameter + FIGO + Karnofsky, data = Ova)
#'
icfit <- function(formula, data, entry, lambda = 10, nt = 100, tmax,
  nseg = 20, bdeg = 3, pord = 2, nit = 50, tol = 1e-06, tollam = 0.01, 
  kappa = 1e-06, update_lambda = TRUE, ic_update = TRUE, monitor = FALSE) {
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", 
               "weights", "na.action", "offset"), 
             names(mf), 0L)
  tmp <- get_input_icfit(formula, data, entry)
  # Check arguments after formula, data, entry
  checkmate::assert_number(lambda, lower = 0)
  checkmate::assert_number(nt, lower = 0)
  checkmate::assert_number(nseg, lower = 1)
  checkmate::assert_number(bdeg, lower = 0)
  checkmate::assert_number(pord, lower = 1)
  checkmate::assert_number(nit, lower = 0)
  checkmate::assert_number(tol, lower = 0)
  checkmate::assert_number(tollam, lower = 0)
  checkmate::assert_number(kappa, lower = 0)
  checkmate::assert_logical(update_lambda, len = 1)
  checkmate::assert_logical(ic_update, len = 1)
  checkmate::assert_logical(monitor, len = 1)
  
  # Initialise (data in discrete bins, B-spline basis function,
  # starting values for baseline hazard)
  Times <- tmp$Ymat
  X <- cbind(1, tmp$X) # add column of 1's to X (will be penalized)
  initres <- init(Times = Times, X = X, nt = nt, tmax = tmax, nseg = nseg, 
    bdeg = bdeg, pord = pord, kappa = kappa)

  # Need single E-step to obtain first 'estimate' of Y and R
  B <- initres$basis$B
  nb <- ncol(B)
  nx <- ncol(X)
  n <- initres$data$n
  HR <- matrix(1, n, 1)
  h0 <- initres$start$h0 * initres$bins$dt
  H <- outer(rep(h0, nt), c(HR))
  Ic <- initres$data$Ic
  R1 <- initres$data$R1
  dead <- initres$data$dead
  es <- Estep(H, Ic, R1, dead)
  Y <- es$Y
  R <- es$R
  
  # Test out whether data is "exactly observed", relative to nt
  # ic_update will be set to TRUE if there is anyone with an event in possibly multiple bins
  if (ic_update) # only possibly overwrite to FALSE if the user has already set ic_update to TRUE
    ic_update <- any(apply(Ic, 2, sum) > 1 & dead)
  # deb(ic_update)

  # Starting values of cbx
  cbx <- rep(0, nb + nx)
  cbx[1:nb] <- log(h0)
  
  # Penalty matrix
  Pdiff <- initres$penalty$Pdiff
  Pridge <- initres$penalty$Pridge

  # Fit it
  if (monitor)
    cat("Monitoring it, log-lik, abs diff estimates, lambda's, eff dim\n\n")
  
  ptm <- proc.time()
  res <- fitit(Y = Y, R = R, dead = dead, X = X, B = B, Ic = Ic, R1 = R1,
               cbx = cbx, Pdiff = Pdiff, Pridge = Pridge, lambda = lambda,
               nit = nit, tol = tol, tollam = tollam, 
               update_lambda = update_lambda, ic_update = ic_update, 
               monitor = monitor)
  
  # Calculation of Fisher information matrices
  info <- InfoMatrix(res, initres)
  res$Ifull <- info$Ifull
  res$Iloss <- info$Iloss
  res$Itot <- info$Itot
  
  # Pack elements in result
  elapsed <- proc.time() - ptm
  res$data <- initres$data
  res$formula <- formula
  res$bins <- initres$bins
  res$basis <- initres$basis
  res$control <- initres$control
  res$call <- cl
  res$m <- m
  res$dead <- dead
  res$elapsed <- elapsed
  
  class(res) <- "icfit"
  return(res)
}
