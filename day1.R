# install the packages before loading them
library(survival)
library(KMsurv)
library(rio)
library(survminer) # make sure you're using the latest version of R

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

# plots (KM, hazard function, legislative variable survival and hazard function)
png("KM.png", width = 800, height = 600)
plot(fit, 
     main = "Survival of US government agencies",
     xlab = "Time in days",
     ylab = "Probability of surviving (Kaplan-Meier)")
dev.off()

png("H.png", width = 800, height = 600)
ggsurvplot(fit,
           data = agency,
           legend = "none",
           fun = "cumhaz",
           conf.int = TRUE,
           censor = FALSE,
           palette = "darkturquoise",
           title = "Cumulative hazard function of US government agencies",
           xlab = "Time (number of days)",
           ggtheme = theme_bw())
dev.off()

fitleg <- survfit(agencysurv~leg, data = agency)
png("leg.png", width = 800, height = 600)
ggsurvplot(fitleg, 
           data = agency,
           linetype = "solid",
           conf.int = TRUE,
           censor = FALSE,
           palette = c("gold", "darkturquoise"),
           title = "Survival function of US government agencies",
           xlab = "Time (number of days)",
           ggtheme = theme_bw())
dev.off()
fitleg

?ggsurvplot

# log rank test (sts test varname in stata)
survdiff(agencysurv ~ leg, data = agency)

# leg hazard
hazardleg <- -log(fitleg$surv)
png("hazardleg.png", width = 800, height = 600)
plot(fitleg$time, 
     hazardleg, 
     pch = 20,
     main = "Hazard function for US government agencies",
     xlab = "Time in days",
     ylab = "Hazard function")
text(x = 12000, y = 0.6, "leg = 1")
text(x = 9000, y = 1.5, "leg = 0")
dev.off()