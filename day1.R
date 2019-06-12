# don't forget to install the packages before loading them
library(survival)
library(KMsurv)
library(rio)

# import stata dataset
agency <- import("./agency.dta")

# create censorship indicator
agency <- cbind(agency, agency$enddate != "1997-12-31")
colnames(agency)[100] <- "terminated"

# convert logical to numeric
agency$terminated <- as.numeric(agency$terminated)

# number of agency terminations
sum(agency$terminated[!is.na(agency$terminated)])

# number of right-censored
length(agency$terminated) - sum(agency$terminated)

# creating a Surv object
agencysurv <- Surv(agency$enddate - agency$startdat, 
                   event = agency$terminated)

# fitting it 
fit <- survfit(agencysurv~1)
fit
summary(fit)

# plots (there is an automatic way to plot survival data, but I haven't managed to load the package yet)
png("KM.png", width = 800, height = 600)
plot(fit, 
     main = "Survival of US government agencies",
     xlab = "Time in days",
     ylab = "Probability of surviving (Kaplan-Meier)")
dev.off()

# computing the hazard function
hazard <- -log(fit$surv)

png("H.png", width = 800, height = 600)
plot(fit$time, 
     hazard, 
     type = "s",
     main = "Hazard function for US government agencies",
     xlab = "Time in days",
     ylab = "Hazard function")
dev.off()

fitleg <- survfit(agencysurv~leg)
png("leg.png", width = 800, height = 600)
plot(fitleg, 
     main = "Survival of US government agencies",
     xlab = "Time in days",
     ylab = "Probability of surviving (Kaplan-Meier)")
text(x = 15000, y = 0.5, "leg = 1")
text(x = 9000, y = 0.2, "leg = 0")
dev.off()
fitleg

# log rank test (sts test varname in stata)
survdiff(Surv(agency$enddate - agency$startdat, event = agency$terminated) ~ leg)

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
