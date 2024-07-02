## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 4,
                      echo = TRUE, warning = FALSE, message = FALSE)

## ----ova----------------------------------------------------------------------
library(icpack)

## ----dusersummary-------------------------------------------------------------
table(Ova$Diameter)
table(Ova$FIGO)
table(Ova$Karnofsky)

## ----icf1---------------------------------------------------------------------
ic_fit <- icfit(Surv(time, d) ~ factor(Diameter) + FIGO + factor(Karnofsky),
                data = Ova)
ic_fit

## ----coxph1-------------------------------------------------------------------
cox_fit <- coxph(Surv(time, d) ~ factor(Diameter) + FIGO + factor(Karnofsky), data = Ova)
cox_fit

## ----icfcoxph-----------------------------------------------------------------
ic_fit <- icfit(Surv(time, d) ~ Diameter + FIGO + Karnofsky, data = Ova)
cox_fit <- coxph(Surv(time, d) ~ Diameter + FIGO + Karnofsky, data = Ova)
ic_fit
cox_fit

## ----ploticf------------------------------------------------------------------
plot(ic_fit, type="survival", conf.int = FALSE,
     title = "Baseline survival function from icfit")

## ----baselinecoxph------------------------------------------------------------
bh <- basehaz(cox_fit, centered=FALSE)
bh$surv <- exp(-bh$hazard)
plot(c(0, bh$time), c(1, bh$surv), type="s",
     xlab="Time", ylab="Survival", main="Baseline survival functions from coxph and icfit")
nd_baseline <- data.frame(Diameter = 0, FIGO = 0, Karnofsky = 0)
pred_baseline <- predict(ic_fit, newdata = nd_baseline)[[1]]
lines(pred_baseline$time, pred_baseline$Surv, col = "blue")
legend("topright", c("coxph", "icfit"), lwd=1, col=c("black", "blue"), bty="n")

## ----coxbasehaz---------------------------------------------------------------
plot(bh$time, diff(c(0, bh$hazard)),
     xlab="Time", ylab="Hazard", main="Baseline hazard function from coxph")

## ----ploticfit2---------------------------------------------------------------
plot(ic_fit, type = 'hazard', title = 'Baseline hazard function from icfit')

## ----newd---------------------------------------------------------------------
newd <- data.frame(subject = c("A", "B"),
                   Diameter = c(0, 4),
                   FIGO = c(0, 0),
                   Karnofsky = c(10, 6))
newd

## ----predict------------------------------------------------------------------
picf <- predict(ic_fit, newdata=newd)

## ----summarypredicticfit------------------------------------------------------
summary(picf, times=(0:7) * 365.25)

## ----icfitAB------------------------------------------------------------------
library(gridExtra)
plot(picf, type="surv", ylim=c(0, 1), xlab="Days since randomisation",
     title = c("Subject A", "Subject B"), nrow=1)

## ----coxAB--------------------------------------------------------------------
sfAB <- survfit(cox_fit, newdata=newd)
sfA <- sfAB[1]
sfB <- sfAB[2]
# Retrieve separate plots from icfit for comparison
picf1 <- picf[[1]]
class(picf1) <- c("predict.icfit", "data.frame")
plt1 <- plot(picf1, type='surv', ylim=c(0, 1),
             xlab = "Days since randomisation",
             title = "Subject A")
picf2 <- picf[[2]]
class(picf2) <- c("predict.icfit", "data.frame")
plt2 <- plot(picf2, type='surv', ylim=c(0, 1),
             xlab = "Days since randomisation",
             title = "Subject B")

oldpar <- par(mfrow=c(1, 2))
plot(sfA, xlab="Days since randomisation", ylab="Survival", main="Subject A")
lines(plt1$data$time, plt1$data$Surv, col = "blue")
lines(plt1$data$time, plt1$data$lowerSurv, col = "blue")
lines(plt1$data$time, plt1$data$upperSurv, col = "blue")
legend("topright", c("coxph", "icfit"), lwd=1, col=c("black", "blue"), bty="n")
plot(sfB, xlab="Days since randomisation", ylab="Survival", main="Subject B")
lines(plt2$data$time, plt2$data$Surv, col = "blue")
lines(plt2$data$time, plt2$data$lowerSurv, col = "blue")
lines(plt2$data$time, plt2$data$upperSurv, col = "blue")
legend("topright", c("coxph", "icfit"), lwd=1, col=c("black", "blue"), bty="n")

par(oldpar)

