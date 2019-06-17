# install the packages before loading them
library(survival)
library(rio)
library(survminer)
library(bshazard)
library(tidyverse)
library(eha)
library(flexsurv)

# import stata dataset
agency <- import("./agency.dta")

# create censorship indicator
agency <- cbind(agency, agency$enddate != "1997-12-31")
colnames(agency)[100] <- "terminated"
agency$terminated <- as.numeric(agency$terminated)

sum(agency$terminated[!is.na(agency$terminated)])

length(agency$terminated) - sum(agency$terminated)

# creating a Surv object and fitting it
agencysurv <- Surv(agency$enddate - agency$startdat, event = agency$terminated)
fit <- survfit(agencysurv ~ 1)
fit; summary(fit)

# plots (KM, hazard function, smoothed hazard, legislative survival and hazard function)

ggsurvplot(fit,
           data = agency,
           legend = "none",
           conf.int = TRUE,
           censor = FALSE,
           palette = "darkturquoise",
           title = "Survival function of US government agencies (KM)",
           xlab = "Time (number of days)",
           ggtheme = theme_bw()) -> KM

ggsurvplot(fit,
           data = agency,
           legend = "none",
           fun = "cumhaz",
           conf.int = TRUE,
           censor = FALSE,
           palette = "darkturquoise",
           title = "Cumulative hazard function of US government agencies",
           xlab = "Time (number of days)",
           ggtheme = theme_bw()) -> H

fitleg <- survfit(agencysurv ~ leg, data = agency)
fitleg

ggsurvplot(fitleg, 
           data = agency,
           linetype = "solid",
           conf.int = TRUE,
           censor = FALSE,
           palette = c("gold", "darkturquoise"),
           title = "Survival function of US government agencies (KM)",
           xlab = "Time (number of days)",
           ggtheme = theme_bw()) -> leg

ggsurvplot(fitleg, 
           data = agency,
           fun = "cumhaz",
           linetype = "solid",
           conf.int = TRUE,
           censor = FALSE,
           palette = c("gold", "darkturquoise"),
           title = "Cumulative hazard function of US government agencies",
           xlab = "Time (number of days)",
           ggtheme = theme_bw()) -> hazardleg

png("KM.png", width = 600, height = 400)
KM
dev.off()

png("H.png", width = 600, height = 400)
H
dev.off()

png("leg.png", width = 600, height = 400)
leg
dev.off()

png("hazardleg.png", width = 600, height = 400)
hazardleg
dev.off()

png("all.png", width = 750, height = 500)
arrange_ggsurvplots(list(KM, leg, H, hazardleg), 
                    nrow = 2, 
                    ncol = 2) 
dev.off()

# log rank test (sts test varname in stata)
survdiff(agencysurv ~ leg, 
         data = agency)

# hazard graph
png("hazard.png", width = 600, height = 500)
plot(bshazard(agencysurv ~ 1, 
              data = agency, 
              lambda = 30), 
     col = "blue", 
     col.fill = "khaki1",
     main = "Smoothed hazard",
     xlab = "Time in days")
dev.off()

# hazard graph for two given variables (I have to separate leg == 1 and leg == 0)

agency %>% filter(leg == 1) -> agency1 
agency %>% filter(leg == 0) -> agency0

agencysurv1 <- Surv(agency1$enddate - agency1$startdat, 
                    event = agency1$terminated)
agencysurv0 <- Surv(agency0$enddate - agency0$startdat, 
                    event = agency0$terminated)

png("Hcompare.png", width = 600, height = 400)
plot(bshazard(agencysurv0 ~ 1, data = agency0, lambda = 30), 
     conf.int = FALSE,
     col = "blue",
     main = "PH test for leg == 1 and leg == 0")
lines(bshazard(agencysurv0 ~ 1, data = agency0, lambda = 30), 
      col = "gold", conf.int = FALSE, lwd = 2)
lines(bshazard(agencysurv1 ~ 1, data = agency1, lambda = 30), 
      col = "darkturquoise", conf.int = FALSE, lwd = 2)
text(x = 3000, y = 0.00005, "leg = 1")
text(x = 6000, y = 0.00015, "leg = 0")
dev.off()

# parametric survival regression (AFT)
survreg(agencysurv ~ leg, 
        data = agency, 
        dist = "weibull") -> s
summary(s)

# survreg only has AFT
# flexsurvreg has AFT for Weibull and PH for exponential
# eha has PH for both

# parametric survival regression with eha (PH)
exp <- phreg(agencysurv ~ leg + num + exec, 
             data = agency, 
             shape = 1) # exponential as weibull with shape 1
weib <- phreg(agencysurv ~ leg + num + exec, 
              data = agency) #weibull
exp; weib

# compute AIC from maximum log-likelihood in model results
# loglik doc: Vector of length 2. The first component is the maximized loglihood with only...
# ...scale and shape in the model, the second the final maximum, which one is most relevant?
ehaAIC <- function (fit) return(2 * fit$df - 2 * fit$loglik[2])
ehaAIC(exp)
ehaAIC(weib)

# post-estimation graphs
png("exp.png", height = 500, width = 500)
plot(exp)
dev.off()

png("weib.png", height = 500, width = 500)
plot(weib)
dev.off()
?coxph

# Cox PH model
coxph(agencysurv ~ leg + num + com + mem, 
      data = agency) -> cox
cox
coxph(agencysurv ~ leg + num + com + mem + exec, 
      data = agency)
cox.zph(cox) -> zph
zph

par(mfrow = c(1,1))
plot(cox.zph(cox))
plot(zph[2])

png("Schoenfeld.png", width = 700, height = 700)
par(mfrow = c(2,2))
ggcoxzph(cox.zph(cox)) # equivalent with survminer
dev.off()

ggcoxdiagnostics(cox) # martingale residuals

# time-varying covariates
as.numeric(unique(agency$enddate[agency$terminated == 1])) # natural cut points

# wide to long
survSplit(formula = agencysurv ~ leg + num + com + mem + enddate + exec, 
          data = agency,
          id = "agencyid",
          cut = seq(min(as.numeric(agency$enddate)), 
                    max(as.numeric(agency$enddate)), 
                    by = 50),
          end = "enddate",
          start = "startdat",
          event = "terminated") -> long

long$enddate <- as.numeric(long$enddate) - min(as.numeric(long$enddate))

# mem tvc
long$lmem <- long$mem * log(as.numeric(long$enddate))
long$lmem

# fitting the model with the tvc / cluster 
coxph(agencysurv ~ leg + num + com + mem + exec, 
      data = agency) -> notvc
coxph(agencysurv ~ leg + num + com + mem + lmem + exec, 
      data = long, 
      cluster(agencyid)) -> coxtvc
cox.zph(notvc)
cox.zph(coxtvc)

# specify the id when creating Surv object (+ cluster)

agencysurv2 <- Surv(time = as.numeric(agency$startdat),
                    time2 = as.numeric(agency$enddate), 
                    event = agency$terminated)


coxph(agencysurv ~ leg + reorg + bdivided, 
      data = agency) -> test

coxph(agencysurv ~ leg + num + com + mem, 
      data = agency)

# compare Cox PH model with Weibull and exponential
phreg(agencysurv ~ leg + num + com + mem, 
      data = agency, 
      shape = 1)
phreg(agencysurv ~ leg + num + com + mem, 
      data = agency)
