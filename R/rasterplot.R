#' Plot probabilities as a raster`
#'
#' @importFrom ggplot2 geom_raster scale_colour_viridis_c labs theme_minimal
#' @importFrom dplyr mutate
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' 
#' @param icf an object of class 'icfit'
#' @param type a string giving the type of the plot. Accepted choices: 'R' (risk probabilities) and 'Y' (event probabilities) 
#' @param sel a vector of integers for selection of subject (rows of the matrix)
#' @param label character vector containing labels for the individuals to be plotted in selection
#' @param show_label Boolean, whether or not to show the labels
#' @param pow a number, giving he power to which the probabilities will be raised, 
#' to improve the clarity of the plot
#' @param order Boolean, default (TRUE) is to order according to first positive in Y,
#' then first zero in Y, then first zero in R; if FALSE order of occurrence in data is used
#' @param do_plot Boolean, default (TRUE) shows the plot, if FALSE object is returned but not plotted
#' 
#' @returns a ggplot object (Grob)
#' 
#' @examples
#' \donttest{
#' icf <- icfit(Surv(left, right, type='interval2') ~ period + gender + age,
#'   data=drugusers)
#' rasterplot(icf)
#' rasterplot(icf, type = 'R')
#' rasterplot(icf, type = 'Y')
#' rasterplot(icf, pow = 0.05) # very small power basically shows 0/1
#' sel <- c(
#'   11, 18,  # right-censored, event in (L, \infty)
#'   1:2,     # event in (0, R)
#'   115, 133 # event in (L, R)
#' )
#' rasterplot(icf, sel = sel)
#' rasterplot(icf, sel = sel, label = c("e", "p", "g", "c", "m", "n"), show_label = TRUE)
#' rasterplot(icf, sel = sel, label = c("e", "p", "g", "c", "m", "n"), show_label = TRUE,
#'   type = 'Y')
#' }
#' 
#' @export
#' 
rasterplot <- function(icf,
                       type = c('both', 'R', 'Y'),
                       sel = NULL,
                       label = NULL,
                       show_label = FALSE,
                       pow = 0.2,
                       order = TRUE,
                       do_plot = TRUE) {
  # Apply selection and match label
  if (is.null(sel)) sel <- 1:ncol(icf$Y)
  nplot <- length(sel)
  if (is.null(label)) label <- sel
  if (length(label) != nplot) stop("'label' does not have same length as 'sel'")

  # Extract Y and R from icf object
  Y <- icf$Y[, sel] ^ pow
  R <- icf$R[, sel] ^ pow
  colnames(R) <- NULL

  # Order according to right-censored, then interval-censored, according to first left, then right
  if (order) {
    Yfirstpos <- apply(Y, 2, function(x) {u = which(x > 0); ifelse(length(u) > 0, min(u), 0)})
    Yfirstzero <- apply(Y, 2, function(x) {u = which(x == 0); ifelse(length(u) > 0, min(u), 0)})
    Rfirstzero <- apply(R, 2, function(x) {u = which(x == 0); ifelse(length(u) > 0, min(u), 0)})
    ord <- order(Yfirstpos, Yfirstzero, Rfirstzero)
    Y <- Y[, ord]
    R <- R[, ord]
  }
  
  type <- match.arg(type)
  if (!(type == 'R')) { # build plot for Y unless only R requested
    dfY <- melt(Y) |> 
      mutate(Time = .data$Var1 * icf$bins$dt)
    if (show_label)
      dfY$Var2 <- factor(dfY$Var2, labels = label)
    pY <- ggplot(dfY, aes(.data$Time, .data$Var2)) +
      geom_raster(aes(fill = .data$value)) +
      scale_colour_viridis_c() +
      labs(x = "Time",
           y = "",
           title = "Y") + 
      theme_minimal()
  }
  if (!(type == 'Y')) { # build plot for R unless only Y requested
    dfR <- melt(R) |> 
      mutate(Time = .data$Var1 * icf$bins$dt)
    if (show_label)
      dfR$Var2 <- factor(dfR$Var2, labels = label)
    pR <- ggplot(dfR, aes(.data$Time, .data$Var2)) +
      geom_raster(aes(fill = .data$value)) +
      scale_colour_viridis_c(direction = -1) + 
      labs(x = "Time",
           y = "",
           title = "R") + 
      theme_minimal()
  }
  if (type == 'both')
    res <- lemon::grid_arrange_shared_legend(pY, pR, nrow=1, ncol=2)
  if (type == 'Y')
    res <- pY
  if (type == 'R')
    res <- pR
  if (do_plot) plot(res)
  return(res)
}
