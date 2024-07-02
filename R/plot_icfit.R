#' Plot method for an object of class `icfit`
#'
#' @import ggplot2
#' @importFrom graphics lines polygon
#' @importFrom stats contrasts integrate model.matrix model.response pnorm printCoefmat qnorm
#' 
#' @param x The object of class 'icfit' to be plotted
#' @param type Type of plot. Accepted choices: 'hazard' (default), 'cumhazard', 'survival' or 'cumprob'
#' @param conf.int If `TRUE` a 100*(1 - alpha) percent confidence interval is plotted
#' @param fill Fill area between lower and upper
#' @param fillcol The color for filling (default 'lightgrey')
#' @param ylim The y-limits for the plot
#' @param title Optional title string
#' @param xlab Text for x-label
#' @param ylab Text for y-label
#' @param \dots Other arguments to plot (except `type`, which is set to 'l')
#'
#' @return A ggplot grob, containing the plot. Use \code{print()} or \code{plot()} to show it
#' Multiple objects can be combined by using functions in the package \code{gridExtra}.
#'
#' @examples
#' \donttest{
#' icf <- icfit(Surv(left, right, type='interval2') ~ period + gender + age, 
#'              data = drugusers)
#' plot(icf)
#' }
#' 
#' @export

plot.icfit <- function(x, 
                       type = c('hazard', 'cumhazard', 'survival', 'probability'),
                       conf.int = TRUE, ylim = NULL, title = NULL, xlab = NULL, ylab = NULL,
                       fill = TRUE, fillcol = 'lightgrey', ...) {
  
   if (!inherits(x, "icfit"))
      stop("'x' must be an 'icfit' object")
   
   res <- predict.icfit(x)
   type <- match.arg(type)
   plt <- plot(res, type=type, conf.int=conf.int, fill=fill, fillcol=fillcol,
              ylim=ylim, title=title, xlab=xlab, ylab=ylab, ...)
   return(plt)
}      
