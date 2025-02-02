% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_predict_icfit.R
\name{summary.predict.icfit}
\alias{summary.predict.icfit}
\title{Summary method for an object of class `predict.icfit`}
\usage{
\method{summary}{predict.icfit}(object, times, ...)
}
\arguments{
\item{object}{Object of class 'predict.icfit'}

\item{times}{The time points at which to summarize the predicted hazards, cumulative hazards
and survival probabilities, with associated standard errors and confidence intervals}

\item{\dots}{Other arguments to plot}
}
\value{
A data frame (if object was a data frame) or a list of data frames (if object was
a list of data frames) with hazards etc linearly interpolated between the time points used
in the predict function
}
\description{
Summary method for an object of class `predict.icfit`
}
\examples{
\donttest{
icf <- icfit(Surv(left, right, type='interval2') ~ period + gender + age,
      data=drugusers)
pred_icf <- predict(icf)
summary(pred_icf, times=c(0, 30, 183, 365))
}

}
