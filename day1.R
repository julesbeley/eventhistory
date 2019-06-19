# install the packages before loading them
rm(list = ls())
library(ggplot2)
library(reshape2)
library(survival)
library(rio)
library(survminer)
library(bshazard)
library(tidyverse)
library(eha)
library(flexsurv)
library(penalized)
library(glmnet)

# import stata dataset
agency <- import("./agency.dta")

# missing data map

png("missing.png", height = 500, width = 750)
agency %>% 
    is.na %>%
    melt %>% 
    ggplot(data = .,
           aes(x = Var2,
               y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_grey(name = "",
                    labels = c("Present","Missing")) +
    theme_minimal() + 
    theme(axis.text.x  = element_text(angle = 90, vjust = 0.2),
          text = element_text(size = 18)) + 
    labs(x = "Variables in Dataset",
         y = "Rows / observations")
dev.off()

agency <- agency[, -which(names(agency) %in% c("budget", "public", "adahmed", "adasenme"))]

sum(is.na(agency$budget))
sum(is.na(agency$publicl))
sum(is.na(agency$adahmed))
sum(is.na(agency$adasenme))

# create censorship indicator
agency <- cbind(agency, agency$enddate != "1997-12-31")
colnames(agency)[97] <- "terminated"
agency$terminated <- as.numeric(agency$terminated)

sum(agency$terminated[!is.na(agency$terminated)])
length(agency$terminated) - sum(agency$terminated)

# creating a Surv object and fitting it
agencysurv <- Surv(agency$enddate - agency$startdat, event = agency$terminated)
fit <- survfit(agencysurv ~ 1)
fit; summary(fit)

# plots (KM, hazard function, smoothed hazard, legislative survival and hazard function)
ggsurvplot(fit, data = agency, legend = "none",
           conf.int = TRUE, censor = FALSE,
           palette = "darkturquoise",
           title = "Survival function of US government agencies (KM)",
           xlab = "Time (number of days)",
           ggtheme = theme_bw()) -> KM

ggsurvplot(fit, data = agency, legend = "none",
           fun = "cumhaz", conf.int = TRUE,
           censor = FALSE, palette = "darkturquoise",
           title = "Cumulative hazard function of US government agencies",
           xlab = "Time (number of days)",
           ggtheme = theme_bw()) -> H

fitleg <- survfit(agencysurv ~ leg, data = agency)
fitleg

ggsurvplot(fitleg, data = agency, linetype = "solid",
           conf.int = TRUE, censor = FALSE,
           palette = c("gold", "darkturquoise"),
           title = "Survival function of US government agencies (KM)",
           xlab = "Time (number of days)",
           ggtheme = theme_bw()) -> leg

ggsurvplot(fitleg, data = agency, fun = "cumhaz",
           linetype = "solid", conf.int = TRUE,
           censor = FALSE, palette = c("gold", "darkturquoise"),
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

plot(bshazard(agencysurv ~ strata(leg), data = agency, lambda = 30), 
     conf.int = FALSE)

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

# strata, comparing survival for a dummy variable in a Cox PH model
coxph(agencysurv ~ strata(locat2), 
      data = agency) -> strata

# term has very few cases?
plot(survfit(strata))
plot(survfit(strata), 
     col = c("blue", "red", "black", "green", "orange"))

# recode location variable into three categories (low, medium, high)
agency$locat23[agency$locat2 == 1] <- 1
agency$locat23[agency$locat2 == 2 | agency$locat2 == 3] <- 2
agency$locat23[agency$locat2 == 4 | agency$locat2 == 5] <- 3
agency$locat23

coxph(agencysurv ~ strata(locat23), 
      data = agency) -> strata
plot(survfit(strata), 
     col = c("blue", "red", "black"))

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
# loglik doc: vector of length 2
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

colnames(agency)

# Cox PH model
coxph(agencysurv ~ leg + num + com + mem + locat2 + line + term + corp + exec, 
      data = agency) -> cox
cox
plot(survfit(cox))
coxph(agencysurv ~ leg + com + mem + locat23 + term + exec, 
      data = agency) -> bet
bet
plot(survfit(bet))
cox.zph(bet) -> zph
zph

par(mfrow = c(1,1))
plot(cox.zph(cox))
plot(zph[2])

png("Schoenfeld.png", width = 700, height = 700)
par(mfrow = c(2,2))
ggcoxzph(cox.zph(cox)) # equivalent with survminer
dev.off()

ggcoxdiagnostics(cox) # martingale residuals

# compare Cox PH model with Weibull and exponential
phreg(agencysurv ~ leg + num + com + mem, 
      data = agency, 
      shape = 1)
phreg(agencysurv ~ leg + num + com + mem, 
      data = agency)

# time-varying covariates
as.numeric(unique(agency$enddate[agency$terminated == 1])) # natural cut points

# how do you decide on cut points when creating a tvc?
# when is a tvc a good idea?
# why are p-values in Cox PH model with tvc all so high?
# when is a pwe model a good idea?

# wide to long (modify enddate so that it's equivalent to agencysurv[2]?)
# terminated variable is not correct (it is the same for all cut points)
survSplit(formula = agencysurv ~ leg + num + com + mem + exec, 
          data = agency,
          id = "agencyid",
          cut = seq(min(as.numeric(agency$enddate)), 
                    max(as.numeric(agency$enddate)), 
                    by = 365),
          end = "enddate",
          start = "startdat",
          event = "terminated") -> long

long$terminated <- long$agencysurv[,3]
long$time1 <- long$agencysurv[, 2]
long$time0 <- long$agencysurv[, 1]

long <- long[ , c("agencyid", "com", "leg", 
                  "num", "mem", "exec", "time0", 
                  "time1", "agencysurv", "terminated")]
long
# mem tvc
long$lmem <- long$mem * log(as.numeric(long$time1))
long

# fitting the model with the tvc / cluster 
coxph(agencysurv ~ leg + num + com + mem + exec, 
      data = agency) -> notvc
plot(survfit(notvc)) # estimated survival function
coxph(agencysurv ~ leg + num + com + mem + lmem + exec, 
      data = long, cluster(agencyid)) -> coxtvcc # clustering
coxph(agencysurv ~ leg + num + com + mem + lmem + exec, 
      data = long) -> coxtvcnc # no clustering
cox.zph(notvc)
cox.zph(coxtvcc)
cox.zph(coxtvcnc)

# piecewise exponential baseline model (computed as a poisson regression)
glm(terminated ~ leg + num + com + mem + exec, 
    data = long,
    family = "poisson") -> pwe
summary(pwe)

# comparing with stata results
coxph(agencysurv ~ leg + exec + bdivided, data = agency) 

coxph(agencysurv ~ leg + reorg + bdivided, 
      data = agency) -> test

coxph(agencysurv ~ leg + num + com + mem, 
      data = agency)

# final model (show that exec and corp (and mem?) uninteresting)
coxph(agencysurv ~ leg + com + mem + locat23 + term + corp + exec, 
      data = agency) -> bet
bet
plot(survfit(bet))
cox.zph(bet) -> zph
zph

coxph(agencysurv ~ leg + com + locat23 + term, 
      data = agency) -> bet2
bet2
plot(survfit(bet2))
cox.zph(bet2) -> zph2
zph2

#... tvc on PH assumption violator? (locat23)
survSplit(formula = agencysurv ~ ., 
          data = agency,
          cut = seq(min(as.numeric(agency$enddate)), 
                    max(as.numeric(agency$enddate)), 
                    by = 365),
          end = "enddate",
          start = "startdat",
          event = "terminated") -> long

long$terminated <- long$agencysurv[,3]
long$time1 <- long$agencysurv[, 2]
long$time0 <- long$agencysurv[, 1]

long

# locat23 tvc
long$llocat23 <- long$locat23 * as.numeric(long$time1)
long

coxph(agencysurv ~ exec + leg + com + locat23 + term + llocat23, 
      data = long, cluster(agencyid)) -> bet2tvc

bet2 # compare with no tvc
options(scipen = 999)
bet2tvc
cox.zph(bet2tvc)
plot(survfit(bet2tvc))

par(mfrow = c(2,2))
ggcoxzph(cox.zph(bet2tvc))

# apatables to export results to Word
# run frailty model as a test (as a diagnostic -> footnote) but no interpretation of coefficients
# use coxme to run frailty model
coxph(agencysurv ~ leg + com + locat23 + term + llocat23 + frailty(agencyid), 
      data = long, cluster(agencyid)) -> bet2tvc # doesn't work yet

# penalized Cox regression using "penalized" library
par(mfrow = c(1,1))
penalized(response = agencysurv, 
          penalized = ~ .,
          standardize = TRUE,
          data = na.omit(long),
          steps = "Park") -> pen
plotpath(pen, lwd = 2, labelsize = 1)

# the same with glmnet on agency
agencynet <- cbind(agency, agencysurv)
x <- na.omit(agencynet[,-which(names(agency) %in% c("agencyid",
                                                    "bureau",
                                                    "startdat",
                                                    "enddate",
                                                    "year",
                                                    "index",
                                                    "terminated",
                                                    "DO"))])
y <- cbind(x$agencysurv[,1], x$agencysurv[,2])
colnames(y) <- c("time", "status")
x <- x[,-which(names(x) %in% c("agencysurv"))]
x <- as.matrix(x)
net <- glmnet(x, y, family = "cox", alpha = 1, maxit = 500000)
cv <- cv.glmnet(x, y, family = "cox", alpha = 1)
plot(net, xvar = "lambda")
abline(v = log(cv$lambda.min), lty = 2)

# get coefficients for cross-validated estimate
coeff <- as.matrix(coef(cv, s = cv$lambda.min))
coeff <- as.data.frame(cbind(rownames(coeff), coeff))
row.names(coeff) <- NULL
colnames(coeff) <- c("name", "value")
coeff$value <- as.numeric(as.character(coeff$value))
coeff$value[coeff$value < 0.1 & coeff$value > - 0.1] <- NA # setting bounds
coeff <- na.omit(coeff)
coeff
