#' Fills space between two lines in a graph
#'
#' Taken from mstate
#'
#' @export
#' @param x Points on the x-axis
#' @param y1 First set of points on y-axis
#' @param y2 Second set of points on y-axis
#' @param col The color to fill space with
#'
#' @return Nothing

fillplot <- function(x, y1, y2, col) {
  nx <- length(x)
  # add mini-bit of space, this is to incorporate the
  # possibility of a jump at the end
  x <- c(x, x[nx])
  xx <- c(rep(x, c(1, rep(2, nx - 1), 1)), rep(rev(x), c(1, 
    rep(2, nx - 1), 1)))
  yy <- c(rep(y1, rep(2, nx)), rep(rev(y2), rep(2, nx)))
  polygon(xx, yy, col = col, border = NA)
}
