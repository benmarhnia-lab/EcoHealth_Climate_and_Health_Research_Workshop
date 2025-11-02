#-------------------------------------------------------------------------------#
#---------------Introduction to Climate Epidemiology Methods--------------------# 
#-------------------------Eco Health Workshop-----------------------------------#
#-----------------------------Date:10/28/25, 11/04/25---------------------------#
#--Tarik Benmarhnia (tbenmarhnia@health.ucsd.edu) & Yiqun Ma (yim022@ucsd.edu)--#
#-------------------------------------------------------------------------------#
#--------------code adapted from Gasparrini A et al., 2010----------------------#
#----https://github.com/gasparrini/2010_gasparrini_StatMed_Rcode/tree/master----#

library(dlnm)
library(splines)
library(lubridate)

##############################
# LOAD AND PREPARE THE DATASET
##############################
data <- read.csv(here("data", "df-train-test-sf.csv"))
data <- data %>%
  mutate(date = as.Date(date),
         dow = wday(date),
         temp = (Tmax + Tmin)/2) %>%
  select(date, dow, respiratory, temp)

n.year <- length(unique(year(data$date)))
cen <- median(data$temp)

##############################
# CROSSBASIS SPECIFICATION
##############################
# FIXING THE KNOTS AT EQUALLY SPACED VALUES OF TEMPERATURE AND
# AT EQUALLY-SPACED LOG-VALUES OF LAG
ktemp <- equalknots(data$temp,nk=4)
klag <- logknots(30,nk=3)
# CROSSBASIS MATRIX
ns.basis <- crossbasis(data$temp,argvar=list(knots=ktemp),
                       arglag=list(knots=klag),lag=30)
summary(ns.basis)

##############################
# MODEL FIT AND PREDICTION
##############################
ns <- glm(respiratory ~  ns.basis + dow + 
            ns(date,df=n.year*7),	family=quasipoisson(), data)
ns.pred <- crosspred(ns.basis,ns,at=-26:33,cen=cen)

##############################
# RESULTS AND PLOTS
##############################
# 3-D PLOT
plot(ns.pred,zlab="RR",xlab="Temperature")

# SLICES
percentiles <- round(quantile(data$temp,c(0.001,0.05,0.95,0.999)),1)
ns.pred <- crosspred(ns.basis,ns,at=-260:330/10,cen=cen)
plot(ns.pred,var=percentiles,lag=c(0,5,15,28))

# OVERALL EFFECT
plot(ns.pred,"overall",xlab="Temperature",
     main="Overall effect of temperature on respiratory hospitalizations\n(San Francisco 2009-2018)")

# RR AT CHOSEN PERCENTILES VERSUS median (AND 95%CI)
ns.pred$allRRfit[as.character(percentiles)]
cbind(ns.pred$allRRlow,ns.pred$allRRhigh)[as.character(percentiles),]
