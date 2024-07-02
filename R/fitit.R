#' Fit proportional hazard model with smooth baseline hazard and (optional) interval censoring
#'
#' @param Y Events (matrix, number of bins by subjects)
#' @param R Risk sets (matrix, number of bins by subjects)
#' @param dead (Boolean vector, TRUE if event, FALSE if right censored)
#' @param X Covariates (matrix, number of covariates (+1) by subjects)
#' @param B B-spline basis matrix
#' @param Ic Censoring interval per individual, coded as 0/1 (in columns)
#' @param R1 Left truncation interval per individual, coded as 0/1 (in columns)
#' @param cbx Vector of starting values
#' @param Pdiff B-spline part of penalty matrix
#' @param Pridge Ridge part of penalty matrix (for intercept)
#' @param lambda Smoothing parameter (number)
#' @param nit Maximum number of iterations (integer)
#' @param tol Tolerance for final fit
#' @param tollam Tolerance for switching to lambda update
#' @param update_lambda Automatic update of lambda (Boolean)
#' @param ic_update Update risk and event probabilities (Boolean)
#' @param monitor Monitor convergence (Boolean)
#'
#' @return A list with items
#' \item{cbx}{Vector of }
#' \item{ll}{Poisson GLM log-likelihood}
#' \item{lambda}{Final tuning parameter}
#' \item{pen}{Penalty part of penalized log-likelihood}
#' \item{ed}{Effetive dimension of the baseline hazard}
#' \item{nit1}{Number of iterations used in first phase}
#' \item{nit}{Total number of iterations used (first plus second phase)}
#' \item{tollam}{Tolerance used for switching to lambda update}
#' 
fitit <- function(Y, R, dead, X, B, Ic, R1, cbx,
                  Pdiff, Pridge, lambda,
                  nit = 50, tol = 1e-06, tollam = 0.01,
                  update_lambda = FALSE, ic_update = TRUE,
                  monitor = FALSE) {
  
  # Some matrix dimensions
  nt <- nrow(Y)
  nx <- ncol(X)
  nb <- ncol(B)

  # Fit the models
  phase <- 1
  it1 <- 1
  
  for (it in 1:nit) {
    
    # Penalty matrix
    Pen <- lambda * Pdiff + Pridge
    # M-step
    ms <- Mstep(Y, R, X, B, Pen, lambda, cbx)
    # Keep track of difference between successive estimates    
    dc <- max(abs(cbx - ms$cbx))
    if (monitor) 
      cat(it, ms$ll, dc, lambda, ms$lambda, ms$ed, "\n")
    
    cbx <- ms$cbx
    # Go to phase 2 (lambda update)
    if (dc < tollam & update_lambda) phase <- 2
    # Stop iteration when converged
    if (dc < tol) break

    # Update lambda
    if (phase == 2 & update_lambda) {
      dla <- abs((lambda - ms$lambda)/lambda)
      lambda <- ms$lambda
    } else it1 <- it
    
    # One (optional) E-step to update Y and R
    if (!is.null(Ic) & ic_update) {
      es <- Estep(ms$H, Ic, R1, dead)
      Y <- es$Y
      R <- es$R
    }
  }
  
  # Collect and return output
  h0 <- exp(B %*% cbx[1:nb])
  output <- list(cbx = cbx, ll = ms$ll, lambda = lambda, pen = ms$pen, 
    ed = ms$ed, aic = -2 * ms$ll + 2 * ms$ed, h0 = h0, Y = Y, 
    R = R, Mpen = ms$Mpen, nit1 = it1, nit = it, tollam = tollam)
  return(output)
}
