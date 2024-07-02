#' Compute a B-spline basis
#'
#' @export
#' @param x The vector of values for which the basis is to be evaluated
#' @param xl The left boundary of the domain
#' @param xr The right boundary of the domain
#' @param nseg The number of inter-knot segments on the domain
#' @param deg The degree of the B-splines (2 means quadratic, 3 means cubic, and so on)
#' @return A matrix containing the basis
#'
#' @examples
#' x = runif(100)
#' B = bbase(x, 0, 1, 20, 3)

bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3) {
  # Construct B-spline basis
  
  tpower <- function(x, t, p) {
    # Truncated p-th power function
    (x - t) ^ p * (x >= t)
  }
  
  dx <- (xr - xl) / nseg
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
}
