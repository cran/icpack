#' Perform the E-step in the EM algorithm
#'
#' @param H Hazards per individual (in columns)
#' @param Ic Censoring interval per individual, coded as 0/1 (in columns)
#' @param R1 Left truncation interval per individual, coded as 0/1 (in columns)
#' @param dead Boolean vector (TRUE is event, FALSE is right censored)
#'
#' @return A list with two matrices
#' \item{Y}{Expected probability of event per bin per subject}
#' \item{R}{Expected probability of at risk per bin per subject}

Estep <- function(H, Ic, R1, dead) {

  nt <- nrow(H)
  n <- ncol(H)
  
  #!!!! Problems when H > 1 !!!!
  # Will be "dealt with" by setting H[H>1] to 0.9
  H[H > 1] <- 0.9
  
  Su <- apply(1 - H, 2, cumprod) # Survival at end of bin
  Sl <- rbind(1, Su)             # Survival at start of bin
  F <- -apply(Sl, 2, diff)       # Density
  # F <- -diff(Sl)                 # Density
  
  # cat("Estep():\n")
  # deb(t(F[, c(1, 2, 3, 143, 361, 760)]))
  
  # Conditional components
  Y <- R <- matrix(0, nt, n)
  F_sel <- F * Ic

  # deb(t(F_sel[, c(1, 2, 3, 143, 361, 760)]))
  
  # Rescale F for subjects with event
  ds <- colSums(F_sel[, dead]) # rescaling of F_sel
  
  # deb(ds[c(1, 2, 3, 141, 356, 751)])
  
  Y[, dead] <- scale(F_sel[, dead], center=FALSE, scale=ds)
  
  # deb(t(Y[, c(1, 2, 3, 143, 361, 760)]))
  # deb(sum(Y[, 361]))
  # deb(c(0, Y[-nt, 361]))
  # deb(cumsum(c(0, Y[-nt, 361])))
  # deb(1 - cumsum(c(0, Y[-nt, 361])))
  
  R[, dead] <- 1 - apply(rbind(0, Y[-nt, dead]), 2, cumsum)
  # For right-censored subjects Y=0 and R=1 until and including censoring bin
  R[, !dead] <- 1 - rbind(0, Ic[-nt, ])[, !dead]
  # Avoid computational negative values
  R[R < 0] <- 0
  
  # Incorporate left-truncation
  R <- R * R1

  # deb(any(Y<0))
  # deb(any(R<0))
    
  output <- list(Y = Y, R = R)
  output
}
