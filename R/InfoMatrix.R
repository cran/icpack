#' Compute the information matrix
#'
#' @importFrom matrixStats colCumsums
#'
#' @param object Fit obtained from \code{fitit}
#' @param initres Result from \code{init}
#' 
#' @return A list with three items
#' \item{Itrue}{Total of \code{Ifull} and \code{Iloss}, true Fisher information
#' matrix}
#' \item{Ifull}{Full Fisher information matrix}
#' \item{Iloss}{Loss of information due to intervals (missing event times)}
#'
#' @details Three information matrices are computed. One is \code{Ifull} which interprets
#' the imputed \code{R} and \code{Y} data from \code{object} as actual observations.
#' \code{Iloss} gives the loss of information due to imputation. The sum of both matrices
#' is the true information matrix.

InfoMatrix <- function(object, initres) {
  
  # Retrieve information from object and initres
  B <- initres$basis$B
  nb <- ncol(B)
  X <- initres$data$X
  nx <- ncol(X)
  Ic <- initres$data$Ic
  dead <- initres$data$dead
  
  # Last values of event and at-risk matrices
  Y <- object$Y
  R <- object$R

  # Covariates and subject-specific hazards
  cbx <- object$cbx
  cb <- c(cbx[1:nb]) # B-spline part
  cx <- c(cbx[nb + (1:nx)]) # covariate part
  eta0 <- c(B %*% cb)  # log baseline hazard
  f <- c(X %*% cx) # linear predictors
  H <- outer(exp(eta0), exp(f)) # subject-specific hazards

  # Penalty matrices
  Pdiff <- initres$penalty$Pdiff
  Pridge <- initres$penalty$Pridge
  Pen <- Pridge + object$lambda * Pdiff
  
  # Preparations
  nt <- nrow(H)
  n <- ncol(H)
  Su <- apply(1 - H, 2, cumprod) # Survival at end of bin
  Sl <- rbind(1, Su)             # Survival at start of bin
  F <- -diff(Sl)                 # Density
  
  # Vector s0 will contain (for those with event) S_j(l_j) - S_j(r_j)
  s0 <- rep(0, n)
  F_sel <- F * Ic
  s0[dead] <- colSums(F_sel[, dead])

  # Array S1[j, i, ] will contain the derivatives of S_i(tau_{j-1}) - S_i(tau_j)
  S1 <- array(0, c(nt, n, nx + nb))
  H0 <- rbind(rep(0, n), H[-nt, ])
  # Compute S1 for B and X separately
  S0 <- Sl[1:nt, ]
  C <- apply(H0, 2, cumsum)
  
  # Computations for term1 
  for (i in 1:n) {
    if (dead[i]) {
      HB <- c(H[, i]) * B
      HB <- rbind(rep(0, nb), HB[-nt,])
      M1 <- c(H[, i] * S0[, i]) * B 
      M2 <- c(F[, i]) * apply(HB, 2, cumsum)
      S1[, i, 1:nb] <- M1 - M2
      M1 <- outer(c(H[, i] * S0[, i]), X[i, ])
      M2 <- outer(c(F[, i]), X[i, ]) * C[, i]
      S1[, i, nb + (1:nx)] <- M1 - M2
    }
  }
  Ica <- array(Ic, c(nt, n, nx + nb)) # array version of Ic, for multiplication
  s1 <- apply(S1 * Ica, c(2, 3), sum)

  # Derivatives of R and Y
  derR <- derY <- array(0, c(nt, n, nx + nb))
  for (i in 1:n) {
    if (dead[i]) {
      U <- (s0[i] * S1[, i, ] - outer(F[, i], s1[i, ])) / s0[i] ^ 2
      U <- c(Y[, i] != 0) * U
      derY[, i, ] <- U
      # U0 <- apply(rbind(0, U[-nt, ]), 2, cumsum)
      U0 <- matrixStats::colCumsums(rbind(0, U[-nt, ]))
      derR[, i,] <- c(Y[, i] != 0) * U0
    }
  }
  # Derivatives of R and Y, term1
  V <- H * R
  # Tbb <- t(B) %*% (diag(c(rowSums(V))) %*% B)
  # Txx <- t(X) %*% (diag(c(colSums(V))) %*% X)
  Tbb <- t(B) %*% (c(rowSums(V)) * B)
  Txx <- t(X) %*% (c(colSums(V)) * X)
  Tbx <- t(B) %*% V %*% X
  term1 <- rbind(cbind(Tbb, Tbx), cbind(t(Tbx), Txx))
  
  # Computations of term2 with index trick to simulate Kronecker products
  DYb <- derY[, , 1:nb]
  dim(DYb) <- c(nt * n, nb)
  DYx <- derY[, , nb + (1:nx)]
  dim(DYx) <- c(nt * n, nx)
  DRb <- derR[, , 1:nb]
  dim(DRb) <- c(nt * n, nb)
  DRx <- derR[, , nb + (1:nx)]
  dim(DRx) <- c(nt * n, nx)
  DTb <- DYb - c(H) * DRb
  DTx <- DYx - c(H) * DRx
  ib <- rep(1:nt, n)
  ix <- rep(1:n, nt)
  Qbb <- t(B[ib, ]) %*% DTb
  Qbx <- t(B[ib, ]) %*% DTx
  Qxb <- t(X[ix, ]) %*% DTb
  Qxx <- t(X[ix, ]) %*% DTx
  term2 <- rbind(cbind(Qbb, Qbx), cbind(Qxb, Qxx))
  
  nh <- nt * n
  hh <- 1:nh
  jj <- rep(1:nt, n)
  ii <- rep(1:n, each = nt)
  BB <- kronecker(rep(1, n), B)
  XX <- kronecker(X, rep(1, nt))
  Qbb <- t(BB) %*% DTb
  Qbx <- t(BB) %*% DTx
  Qxb <- t(XX) %*% DTb
  Qxx <- t(XX) %*% DTx
  Q <- rbind(cbind(Qbb, Qbx), cbind(Qxb, Qxx))

  # Deliver output
  Ifull <- term1 + Pen
  Ifull <- (Ifull + t(Ifull)) / 2 # Ifull has to be symmetric
  Iloss <- -Q
  Iloss <- (-Q - t(Q)) / 2 # Iloss has to be symmetric
  Itot <- Ifull - Iloss
  info <- list(Ifull = Ifull, Iloss = Iloss, Itot = Itot)
  return(info)
}