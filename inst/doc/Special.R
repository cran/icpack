## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 4,
                      echo = TRUE, warning = FALSE, message = FALSE)

## ----icpack-------------------------------------------------------------------
library(icpack)
library(ggplot2)

## ----icfit_cst----------------------------------------------------------------
ic_fit_cst <- icfit(Surv(time, d) ~ 1, data = Ova,
                    bdeg=1, pord=1, lambda=1e+06, update_lambda=FALSE)
ic_fit_cst

## ----plot_icfit_cst-----------------------------------------------------------
plot(ic_fit_cst)

## ----predict_icfit_cst--------------------------------------------------------
pic_fit_cst <- predict(ic_fit_cst)[, c("time", "haz", "sehaz", "lowerhaz", "upperhaz")]
head(pic_fit_cst)
tail(pic_fit_cst)

## ----DTrate-------------------------------------------------------------------
D <- sum(Ova$d)
T <- sum(Ova$time)
rate <- D / T
lograte <- log(rate)
selograte <- 1 / sqrt(D)
lowerlograte <- lograte - qnorm(0.975) * selograte
upperlograte <- lograte + qnorm(0.975) * selograte
data.frame(D=D, T=T, rate=rate, lower=exp(lowerlograte), upper=exp(upperlograte))

## ----icfit_cst_covs-----------------------------------------------------------
ic_fit_cst <- icfit(Surv(time, d) ~ Diameter + FIGO + Karnofsky, data = Ova,
                    bdeg=1, pord=1, lambda=1e+06, update_lambda=FALSE)
ic_fit_cst

## ----Poissonreg---------------------------------------------------------------
Ova$logtime <- log(Ova$time)
pois_fit <- glm(d ~ Diameter + FIGO + Karnofsky + offset(logtime), family=poisson(), data = Ova)
spois_fit <- summary(pois_fit)
spois_fit

## ----baselineCI---------------------------------------------------------------
lograte <- spois_fit$coefficients[1, 1]
selograte <- spois_fit$coefficients[1, 2]
lowerlograte <- lograte - qnorm(0.975) * selograte
upperlograte <- lograte + qnorm(0.975) * selograte
data.frame(rate=exp(lograte), lower=exp(lowerlograte), upper=exp(upperlograte))

## ----plot_icfit_cst_covs------------------------------------------------------
plot(ic_fit_cst, xlab="Days since randomisation", ylab="Hazard",
     title="Constant baseline hazard") + ylim(0, 0.0052)

## ----icfit_nt500--------------------------------------------------------------
ic_fit_cst <- icfit(Surv(time, d) ~ 1, data = Ova, nt=500,
                    bdeg=1, pord=1, lambda=1e+06, update_lambda=FALSE)
pic_fit_cst <- predict(ic_fit_cst)[, c("time", "haz", "sehaz", "lowerhaz", "upperhaz")]
head(pic_fit_cst)
tail(pic_fit_cst)

## ----icfit_gompertz-----------------------------------------------------------
ic_fit_gompertz <- icfit(Surv(time, d) ~ 1, data = Ova,
                         bdeg=1, pord=2, lambda=1e+06, update_lambda=FALSE)
ic_fit_gompertz

## ----plot_gompertz------------------------------------------------------------
plot(ic_fit_gompertz) + ylim(0, 0.00125)

## ----plot_gompertz_logy-------------------------------------------------------
plot(ic_fit_gompertz) + scale_y_continuous(trans='log10')

## ----icfit_weibull------------------------------------------------------------
Ova$logtime <- log(Ova$time)
ic_fit_weibull <- icfit(Surv(logtime, d) ~ 1, data = Ova,
                         bdeg=1, pord=2, lambda=1e+06, update_lambda=FALSE)
ic_fit_weibull

## ----plot_icfit_weibull-------------------------------------------------------
plot(ic_fit_weibull) + scale_y_continuous(trans='log')

## ----pweib--------------------------------------------------------------------
pw <- plot(ic_fit_weibull)$data
head(pw)[, 1:5]
tail(pw)[, 1:5]
pw$u <- pw$time
pw$time <- exp(pw$u)
pw$haz <- pw$haz / pw$time
pw$lowerhaz <- pw$lowerhaz / pw$time
pw$upperhaz <- pw$upperhaz / pw$time
head(pw)[, 1:5]
tail(pw)[, 1:5]

## ----plot_pweib---------------------------------------------------------------
plt <- ggplot(pw) +
  geom_ribbon(aes(x = pw$time, ymin = pw$lowerhaz, ymax = pw$upperhaz),
              fill = 'lightgrey') +
  geom_line(aes(x = pw$time, y = pw$lowerhaz), colour = 'blue') +
  geom_line(aes(x = pw$time, y = pw$upperhaz), colour = 'blue') +
  geom_line(aes(x = pw$time, y = pw$haz), size = 1, colour = 'blue', lty = 1) +
  ylim(0, 0.002) + 
  ggtitle("Weibull hazard fit") + xlab("Days since randomisation") +
  ylab("Hazard") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
plt

## ----flexsurvreg--------------------------------------------------------------
library(flexsurv)
fsr <- flexsurvreg(Surv(time, d) ~ 1, data = Ova, dist = "weibull")
fsr

alp <- exp(fsr$coef[1])
mu <- exp(fsr$coef[2])
lam <- 1 / mu^alp
hw <- function(t, alp, lam) alp * lam * t^(alp - 1)
maxt <- max(Ova$time)
pw$haz_survreg <- hw(pw$time, alp, lam)

plt <- plt +
  geom_line(aes(x = pw$time, y = pw$haz_survreg), size = 1, lty = 2)
plt

