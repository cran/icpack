#'  Function for fitting proportioal hazard model with baseline hazard
#'
#' @param Y Expected events (matrix)
#' @param R Expected risk sets (matrix)
#' @param X Covariates (matrix)
#' @param B B-spline basis
#' @param Pen Penalty matrix
#' @param lambda Smoothing parameter (number)
#' @param cbx Current coefficient estimates
#'
#' @return An object with fields:
#' \code{H} = hazards (matrix),
#' \code{cbx} = coefficient estimates (vector),
#' \code{lambda} = proposal for new lambda,
#' \code{ed} = effective dimension,
#' \code{G} = G matrix,
#' \code{ll} = log-likelihood,
#' \code{pen} = penalized part of log-likelihood,
#' \code{Mpen} = penalized M matrix
#'
Mstep <- function(Y, R, X, B, Pen, lambda, cbx) {
  
  nb <- ncol(B)
  nx <- ncol(X)

  # Compute hazard, expected values and residuals
  cb <- c(cbx[1:nb])
  cx <- c(cbx[nb + (1:nx)])
  eta0 <- c(B %*% cb)  # log baseline hazard
  f <- c(X %*% cx)  # covariate contribution to log hazard
  Eta <- outer(eta0, f, "+")
  H <- exp(Eta)
  Mu <- R * H
  Res <- Y - Mu

  # cat("Mstep():\n")
  # deb(eta0)
  # deb(exp(eta0))
  
  # Likelihood
  ll <- sum(Y * Eta - Mu)
  pen <- lambda * t(cb) %*% Pen[1:nb, 1:nb] %*% cb/2
  
  # Construct equation system
  wb <- rowSums(Mu)
  wx <- colSums(Mu)
  Bt <- t(B)
  # qb <- Bt %*% rowSums(Res)
  qb <- rowSums(Res) %*% B
  qx <- t(X) %*% colSums(Res)
  Mbb <- Bt %*% (wb * B)
  Mxx <- t(X) %*% (wx * X)
  Mbx <- Bt %*% Mu %*% X
  M <- rbind(cbind(Mbb, Mbx), cbind(t(Mbx), Mxx))
  
  # Solve for new coeffcients and check convergence
  Mpen <- M + Pen
  cnew <- solve(Mpen, c(qb, qx) + M %*% cbx)
  G <- solve(Mpen, M)
  ed <- sum(diag(G)[1:nb])
  cbx <- cnew
  u <- cbx[1:nb]
  v <- t(u) %*% Pen[1:nb, 1:nb] %*% u / lambda
  lambdanew <- max(1e-4, min(1e+06, ed / v)) # Because of error that Mpen was computationally singular
  
  output <- list(H = H, cbx = cbx, lambda = lambdanew, ed = ed, 
                 G = G, ll = ll, pen = pen, Mpen = Mpen)
}
