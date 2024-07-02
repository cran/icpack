#' This function prepares response + feature matrices, borrowed from icenReg
#' @importFrom stats model.frame na.pass
#' @param frml Formula object
#' @param data Data frame
#' @param gety Boolean, do we need to get y as well as x
#' @noRd
make_xy = function(frml, df, gety=TRUE) {
  ans <- list()
  # Making a model.frame
  mod_frame <-  model.frame(formula = eval(frml), 
                            data = as.data.frame(df), 
                            # NA's allowed in y (right censoring) but not x
                            # We check x for nas manually
                            na.action = na.pass)
  
  # Get x from model frame
  # trapping weird problem happens with diag_covar if rightside is ~0 
  x <- try(model.matrix(eval(frml), df), silent = TRUE)
  if(inherits(x, "try-error")) {
    if(frml[[3]] == "0") {
      # In this case, n x 0 matrix is needed
      x <- matrix(0, nrow = nrow(df), ncol = 0)
    }
  }
  if(any(is.na(x))){ stop("Not allowed to have NAs for predictors") }
  
  # Remove intercepts
  if('(Intercept)' %in% colnames(x)){	
    ind = which(colnames(x) == '(Intercept)')
    x <- x[,-ind, drop = F]
  }
  ans$x <- x
  ans$xNames = colnames(x)
  
  # Get y
  if (gety) {
    base_y = model.response(mod_frame)
    
    yMat = makeIntervals(base_y, mod_frame)
    ans$y = yMat
  }
  return(ans)
}

#' This function prepares response + feature matrices, borrowed from icenReg
#' further adapted (see comments) to be used in predict.icfit
#' @importFrom stats model.frame na.pass
#' @param frml Formula object
#' @param data Data frame
#' @noRd
make_xy_pred = function(frml, df) { # adapted from make_xy, gety set to FALSE
  ans <- list()
  # Making a model.frame not needed
  # mod_frame <-  model.frame(formula = eval(frml), 
  #                           data = as.data.frame(df), 
  #                           # NA's allowed in y (right censoring) but not x
  #                           # We check x for nas manually
  #                           na.action = na.pass)
  
  # Get x from model frame
  # trapping weird problem happens with diag_covar if rightside is ~0 
  x <- try(model.matrix(eval(frml), df), silent = TRUE)
  if(inherits(x, "try-error")) {
    if(frml[[2]] == "0") { # changed from frml[[3]] to frml[[2]] (since response has been stripped from frml)
      # In this case, n x 0 matrix is needed
      x <- matrix(0, nrow = nrow(df), ncol = 0)
    }
  }
  if(any(is.na(x))){ stop("Not allowed to have NAs for predictors") }
  
  # Remove intercepts
  if('(Intercept)' %in% colnames(x)){	
    ind = which(colnames(x) == '(Intercept)')
    x <- x[,-ind, drop = F]
  }
  ans$x <- x
  ans$xNames = colnames(x)
  
  # Get y not needed
  # if (gety) {
  #   base_y = model.response(mod_frame)
  #   
  #   yMat = makeIntervals(base_y, mod_frame)
  #   ans$y = yMat
  # }
  return(ans)
}

#' Function makeIntervals used by makexy
#' @param y Matrix of left and right time endpoints
#' @param mf Model frame
#' @noRd
makeIntervals <- function(y, mf){
    yMat <- as.matrix(y)[, 1:2]
  if(is(y, 'Surv')){
    yy <- mf[,1]
    status <- yy[, ncol(yy)] # last column contains status variable
    rightCens <- status == 0
    yMat[rightCens, 2] <- Inf
    exact <- status == 1
    yMat[exact, 2] = yMat[exact, 1]
  }
  else{
    rightIsNA <- is.na(yMat[,2])
    yMat[rightIsNA,2] <- Inf
  }
  storage.mode(yMat) <- 'double'
  return(yMat)
}
