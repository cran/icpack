## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 4,
                      echo = TRUE, warning = FALSE, message = FALSE)

## ----icpack-------------------------------------------------------------------
library(icpack)
library(survival)

## ----dusersummary-------------------------------------------------------------
table(drugusers$period)
table(drugusers$gender)
summary(drugusers$age)

## ----icf1---------------------------------------------------------------------
icf <- icfit(Surv(left, right, type = "interval2") ~ period + gender + age, data = drugusers)

## ----printicf-----------------------------------------------------------------
print(icf)

## ----ploticf------------------------------------------------------------------
plot(icf, title = 'Drug users')

## ----ploticfit2---------------------------------------------------------------
plot(icf, type = 'cumhazard', title = 'Drug users')
plot(icf, type = 'survival', title = 'Drug users')

## ----predict------------------------------------------------------------------
newd <- data.frame(period=1:4, gender=1, age=25)
newd$period <- factor(newd$period, levels=1:4, labels=levels(drugusers$period))
newd$gender <- factor(newd$gender, levels=1:2, labels=levels(drugusers$gender))
picf <- predict(icf, newdata=newd)

## ----summarypredicticfit------------------------------------------------------
summary(picf, times=seq(0, 60, by=12))

## ----plotpredicticfit---------------------------------------------------------
library(gridExtra)
plot(picf, type="surv", ylim=c(0, 1),
     xlab = "Months since start intravenous drug use",
     title = levels(drugusers$period))

