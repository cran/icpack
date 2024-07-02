#' Plot method for an object of class `predict.icfit`
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon xlab ylab ggtitle ylim xlim theme theme_bw
#' @importFrom graphics lines polygon
#' @importFrom stats contrasts integrate model.matrix model.response pnorm printCoefmat qnorm
#' @importFrom rlang .data
#' 
#' 
#' @param x The object of class 'predict.icfit' to be plotted
#' @param type Type of plot. Accepted choices: 'hazard' (default), 'cumhazard', 'survival' or 'probability'
#' @param conf.int If `TRUE` a 100*(1 - alpha) percent confidence interval is plotted
#' @param fill Fill area between lower and upper
#' @param fillcol The color for filling (default 'lightgrey')
#' @param ylim The y-limits for the plot
#' @param title Optional title string, or, if x is a list, obtained from 'predict.icfit' using `newdata`,
#' a vector of title strings
#' @param xlab Text for x-label
#' @param ylab Text for y-label
#' @param selection If x is a list, obtained from 'predict.icfit' using `newdata`, then a
#' vector containing the subset of list elements to be plotted, default is to plot
#' all elements of the list
#' @param nrow If x is a list, obtained from 'predict.icfit' using `newdata`, then a number
#' specifying the number of rows to plot; default the square root of the number of list
#' elements to be plotted
#' @param ncol If x is a list, obtained from 'predict.icfit' using `newdata`, then a number
#' specifying the number of columns to plot; default the square root of the number of list
#' elements to be plotted
#' @param do_plot Boolean indicating whether or not to actually plot (default is TRUE)
#' @param ... other graphical parameters to be passed on
#'
#' @return A ggplot grob, containing the plot. Use \code{print()} or \code{plot()} to show it
#' Multiple objects can be combined by using functions in the package \code{gridExtra}.
#'
#' @examples
#' \donttest{
#' icf <- icfit(Surv(left, right, type='interval2') ~ period + gender + age, data=drugusers)
#' pred_icf <- predict(icf)
#' plot(pred_icf)
#' library(ggplot2)
#' plot(icf) + xlim(0, 200) + ylim(0, 0.05)
#' ndata <- drugusers[1:4, ]
#' pred_nd_icf <- predict(icf, newdata=ndata)
#' plot(pred_nd_icf) # plot all four
#' plot(pred_nd_icf[[2]]) # plot only the second
#' plot(pred_nd_icf, type = "cumhazard") # plot four cumulative hazard curves
#' plot(pred_nd_icf[[3]], type = "prob", ylim = c(0, 1)) # plot probability curve for nr 3
#' plot(pred_nd_icf[[4]], type = "surv", ylim = c(0, 1)) # plot survival curve for nr 4
#' }
#' 
#' @export
#' 
plot.predict.icfit <-
  function(x, 
           type = c('hazard', 'cumhazard', 'survival', 'probability'),
           conf.int = TRUE, 
           fill = TRUE, fillcol = 'lightgrey', ylim = NULL,
           title = NULL, xlab = NULL, ylab = NULL, 
           selection = NULL, nrow = NULL, ncol = NULL,
           do_plot = TRUE, ...)
{
  
  if (!inherits(x, "predict.icfit")) stop("object not of class 'predict.icfit'")
  
  type <- match.arg(type)
  
  if (inherits(x, "data.frame"))
  {
    sp <- singleplot(
      x, type = type, conf.int = conf.int,
      fill = fill, fillcol= fillcol, ylim = ylim,
      title = title, xlab = xlab, ylab = ylab)
    if (do_plot) plot(sp)
    return(sp)
  }
  else if (inherits(x, "list")) {
      plts <- listplot(
        x, type = type, conf.int = conf.int,
        fill = fill, fillcol= fillcol, ylim = ylim,
        title = title, xlab = xlab, ylab = ylab,
        selection = selection)
      
      # Put in agrob
      agrob <- gridExtra::arrangeGrob(grobs = plts, nrow = nrow, ncol = ncol)
      
      if (do_plot) plot(agrob)
      
      return(agrob)
  }
}

singleplot <-
  function(x, 
           type = c('hazard', 'cumhazard', 'survival', 'probability'),
           conf.int = TRUE, 
           fill = TRUE, fillcol = 'lightgrey', ylim = NULL,
           title = NULL, xlab = NULL, ylab = NULL, ...)
{
  if (!inherits(x, "predict.icfit")) stop("object not of class 'predict.icfit'")

  res <- as.data.frame(x)
  type <- match.arg(type)
  if (is.null(xlab)) xlab <- 'Time'
  if (is.null(title)) title <- ''
  
  # Make plot with options
  if (type == 'hazard') {
    if (is.null(ylab)) ylab <- 'Hazard'
    plt = ggplot(res, ...)
    if (conf.int & fill) plt = plt +
      geom_ribbon(aes(x = .data$time, ymin = .data$lowerhaz, ymax = .data$upperhaz),
                  fill = fillcol)
    if (conf.int) {
      plt = plt + 
        geom_line(aes(x = .data$time, y = .data$lowerhaz), colour = 'blue') +
        geom_line(aes(x = .data$time, y = .data$upperhaz), colour = 'blue')
    }
    plt = plt + geom_line(aes(x = .data$time, y = .data$haz), size = 1,
                          colour = "blue", lty = 1)
  }
  
  if (type == 'cumhazard') {
    if (is.null(ylab)) ylab <- 'Cumulative hazard'
    plt = ggplot(res, ...)
    if (conf.int & fill) plt = plt +
      geom_ribbon(aes(x = .data$time, ymin = .data$lowerHaz, ymax = .data$upperHaz),
                  fill = fillcol)
    if (conf.int) {
      plt = plt + geom_line(aes(x = .data$time, y = .data$lowerHaz), colour = 'blue') +
        geom_line(aes(x = .data$time, y = .data$upperHaz), colour = 'blue')
    }
    plt = plt + geom_line(aes(x = .data$time, y = .data$Haz), size = 1,
                          colour = "blue", lty = 1)
  }
  
  if (type == 'survival') {
    if (is.null(ylab)) ylab <- 'Survival'
    plt = ggplot(res, ...)
    if (conf.int & fill) plt = plt +
      geom_ribbon(aes(x = .data$time, ymin = .data$lowerSurv, ymax = .data$upperSurv),
                  fill = fillcol)
    if (conf.int) {
      plt = plt + geom_line(aes(x = .data$time, y = .data$lowerSurv), colour = 'blue') +
        geom_line(aes(x = .data$time, y = .data$upperSurv), colour = 'blue')
    }
    plt = plt + geom_line(aes(x = .data$time, y = .data$Surv), size = 1,
                          colour = "blue", lty = 1)
  }
  
  if (type == 'probability') {
    if (is.null(ylab)) ylab <- 'Probability'
    plt = ggplot(res, ...)
    if (conf.int & fill) plt = plt +
      geom_ribbon(aes(x = .data$time, ymin = 1- .data$upperSurv, ymax = 1 - .data$lowerSurv),
                  fill = fillcol)
    if (conf.int) {
      plt = plt + 
        geom_line(aes(x = .data$time, y =  1- .data$upperSurv), colour = 'blue') +
        geom_line(aes(x = .data$time, y = 1 - .data$lowerSurv), colour = 'blue')
    }
    plt = plt + geom_line(aes(x = .data$time, y = 1 - .data$Surv), size = 1,
                          colour = "blue", lty = 1)
  }
  
  plt = plt + ggtitle(title) + xlab(xlab) + ylab(ylab) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  if (!is.null(ylim))
    plt <- plt + ylim(ylim)
  
  return(plt)
}

listplot <-
  function(x, 
           type = c('hazard', 'cumhazard', 'survival', 'probability'),
           conf.int = TRUE, 
           fill = TRUE, fillcol = 'lightgrey', ylim = NULL,
           title = NULL, xlab = NULL, ylab = NULL, 
           selection = NULL, ...)
{
  
  if (!inherits(x, "predict.icfit")) stop("object not of class 'predict.icfit'")
    
  nsubjs <- length(x)
  
  # If selection present, check and select
  if (!is.null(selection)) {
    if (!all(selection %in% 1:nsubjs)) stop("too many selected")
    x <- x[selection]
  }
  nsubjs <- length(x)
  
  type <- match.arg(type)
  if (is.null(xlab)) xlab <- 'Time'
  if (is.null(title)) title <- ""

  # Set up empty list and fill list with plots calling singleplot()
  plts <- list()
  
  for (i in 1:nsubjs) {
    res <- x[[i]]
    class(res) <- c("predict.icfit", "data.frame")
    plts[[i]] <- singleplot(
      res, type = type, conf.int = conf.int,
      fill = fill, fillcol= fillcol, ylim = ylim,
      title = title, xlab = xlab, ylab = ylab)

  }
  return(plts)
}
