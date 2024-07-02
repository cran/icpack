#' Get and check input of icfit
#'
#' @importFrom survival Surv
#' @param formula A formula object with response of the left of a ~ operator
#' and terms on the right. The response must be a survival object as returned
#' by the `Surv` function, with type either right', 'counting' or 'interval2'
#' @param data A data frame in which to interpret the variable names in the
#' `formula`
#' @param entry When appropriate, a vector of entry (left truncation) times,
#' or a string indicating the column name in `data` containing entry times;
#' only used if Surv object is of type 'interval2'
#'
#' @return A list with items
#' \item{Ymat}{Matrix (number of subjects x 3) containing entry, left and right
#' hand of intervals}
#' \item{X}{Matrix (number of subjects x number of covariates + 1) with design
#' matrix of covariates}
#'
#' @export
#' 
get_input_icfit <- function(formula, data, entry) {

  # Borrowing make_xy from icenReg package, which is a function that takes a formula and data, and returns a list with the response and the design matrix
  
  # reg_items is a list with the response and the design matrix
  reg_items <- make_xy(formula, data)
  
  # Y is the response and X is the design matrix
  Y <- reg_items$y
  X <- reg_items$x
  
  if (is.matrix(X)) # xnms is the names of the columns of the design matrix, if X is a matrix
    xnms <- colnames(X) 
  else # if X is not a matrix, then xnms is the names of the columns of the design matrix
    xnms <- as.character(formula[[3]])
  # n is the number of rows of the response
  n <- nrow(Y)
  
  # If entry is missing, then entry is a vector of 0's of length n,
  # else if entry is a character, then entry is the data of the column with that name
  if (missing(entry)) 
    entry <- rep(0, n) 
  else if (is.character(entry)) 
    entry <- data[, entry]
  # If entry is a single number, then entry is a vector of that number repeated n times  
  # Check that entry is numeric and has length n
  if (length(entry) == 1) entry <- rep(entry, n)
  checkmate::assert_numeric(entry, len=n)
  
  # If any entry is greater than the left end of the interval, then stop
  if (any(entry > Y[, 1]))
    stop("entry should be before left end of intervals")
  
  # Add vector of entry times to the left of matrix Y
  Ymat <- cbind(entry, as.matrix(Y)[, 1:2])

  return(list(Ymat = Ymat, X = X))
}
