# install the packages before loading them
library(survival)
library(rio)
library(survminer) # make sure you're using the latest version of R
library(bshazard)
library(eha)

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
fit <- survfit(agencysurv~1)
fit
summary(fit)

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

plot(bshazard(agencysurv~1, data = agency), 
     col = "blue", 
     col.fill = "gold",
     main = "Smoothed hazard",
     xlab = "Time in days")

fitleg <- survfit(agencysurv~leg, data = agency)
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
arrange_ggsurvplots(list(KM, leg, H, hazardleg), nrow = 2, ncol = 2) 
dev.off()

# log rank test (sts test varname in stata)
survdiff(agencysurv ~ leg, data = agency)

# parametric survival regression (AFT)
survreg(agencysurv ~ leg, data = agency, dist = "weibull") -> s
summary(s)

# survreg only has AFT
# flexsurvreg has AFT for Weibull and PH for exponential
# eha has PH for both

# parametric survival regression with eha (PH)
phreg(agencysurv ~ leg, data = agency)
