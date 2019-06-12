library(survival)
library(KMsurv)
library(rio)
agency <- import("./agency.dta")
agency$enddate-agency$startdat # dur



# to get rid of attributes
for (i in (1:ncol(agency))) {
    agency[,i] <- as.vector(agency[,i])
}