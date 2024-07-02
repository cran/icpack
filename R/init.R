#' Generate a discrete IC object
#'
#' @param Times The (possibly interval censored) survival data, in a matrix
#' @param X The design matrix containing covariates
#' @param nt The number of bins for discretization
#' @param tmax The end of time domain (default 1.01 times largest observation)
#' @param nseg The number of B-spline segments
#' @param bdeg The degree of the B-splines
#' @param pord The order of the differences used in the penalty
#' @param kappa Ridge parameter (number)
#'
#' @return A list with items
#' \item{data}{List containing the original data as well as the binned data}
#' \item{bins}{List with information on bins used}
#' \item{basis}{List containing the B-spline matrix}
#' \item{start}{List containing information on starting values}
#' \item{penalty}{List containing Pdiff and Pridge}
#' \item{control}{List with information on control of B-spline basis}
#'
#' @export
#' 
init <- function(Times, X, nt, tmax, nseg = 20, bdeg = 3, pord = 2, kappa = 1e-6) {
  
  entry <- Times[, 1]
  left <- Times[, 2]
  right <- Times[, 3]
  nx <- ncol(X)  # number of covariates
  n <- length(right)  # sample size
  
  # Discretize and set up bins
  maxt <- max(c(left, right[!is.infinite(right)]))
  if (missing(tmax)) tmax <- maxt * 1.01
  checkmate::assert_numeric(tmax, lower=tmax)
  
  dt <- tmax / nt
  j0 <- floor(entry / dt) + 1
  j1 <- floor(left / dt) + 1
  j2 <- floor(right / dt) + 1
  
  # Discriminate between dead and censored
  dead <- j2 != Inf
  j2 <- pmin(nt, j2)
  
  # Code intervals as 0/1 matrix
  ones <- rep(1, nt)
  II <- outer(1:nt, rep(1, n))
  Ic1 <- outer(ones, j1) <= II
  Ic2 <- outer(ones, j2) >= II
  Ic <- Ic1 * Ic2
  R1 <- outer(ones, j0) <= II

  # Basis functions and penalty matrix
  tt <- (1:nt) - 0.5
  B <- bbase(tt, xl = 0, xr = nt, nseg = nseg, deg = bdeg)
  nb <- ncol(B)

  # Initialise regression coefficients with 0, baseline hazard
  # as constant with value = (number of events) / (exposure time)
  ndead <- sum(dead)
  exptime <- sum(right[dead]) + sum(left[!dead])
  h0 <- ndead / sum(exptime)

  # Penalty matrix
  E <- diag(nb)
  D <- diff(E, diff = pord)
  Pdiff = Pridge = 0 * diag(nb + nx)
  Pdiff[1:nb, 1:nb] = t(D) %*% D
#  diag(Pridge)[nb + (1:nx)] = kappa
  diag(Pridge) = kappa
  
  return(list(data = list(Times = Times, Ic = Ic, R1 = R1, 
    X = X, dead = dead, n = n), bins = list(tmax = tmax, 
    nt = nt, dt = tmax / nt), basis = list(B = B), start = list(h0 = h0, 
    cb = rep(0, nx)), penalty = list(Pdiff = Pdiff, Pridge = Pridge),
    control = list(nseg = nseg, bdeg = bdeg, pord = pord)))
}
