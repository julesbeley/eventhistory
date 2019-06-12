library(survival)
library(KMsurv)
library(rio)
agency <- import("./agency.dta")
agency <- cbind(agency, agency$enddate != "1997-12-31")
colnames(agency)[100] <- "terminated"
agency$terminated <- as.numeric(agency$terminated)
agencysurv <- Surv(agency$enddate - agency$startdat, event = agency$terminated)
fit <- survfit(agencysurv~1)
plot(fit)
sum(agency$terminated[!is.na(agency$terminated)])
table(agency$term)
length(agency$terminated) - sum(agency$terminated)
agencysurv
