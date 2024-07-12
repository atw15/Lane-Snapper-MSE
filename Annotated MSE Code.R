#Load Libraries 
library(snow)
library(snowfall)
library(MSEtool)
library(mvtnorm)
library(DLMtool)
library(readxl)
library(openxlsx)
library(ggplot2)

#Update Packages
update.packages("DLMtool")
update.packages("MSEtool")

# Example Function Usage
LstepCC3(plot=TRUE, Data = Atlantic_mackerel)

# Check Documentation
?matlenlim

# Create new Data object and set parameters
newdata <- new('Data')
newdata@vbLinf <- 50
newdata@L50 <- c(23.5, 24.5)
Atlantic_mackerel@vbLinf <- 155.23

# Run MSE to ensure DLMtool is working correctly
myMSE <- runMSE()

# Visualize MSE results
Pplot(myMSE)

# Import SEDAR Operating Model
NOAAOM <- XL2OM("Lane_Snapper_GOM_NOAA")

# Import Custom Operating Model
OMA <- XL2OM("MyLSOM")

# Add Truncated Parameters to Custom OM
OMA@Linf <- c(35.7, 36.4)
OMA@K <- c(0.39, 0.45)
OMA@t0 <- c(-0.85, -0.61)

# Define Management Procedures (MPs)
matlenlim8 <- function (x, Data, reps, plot = FALSE) {
  rec <- new("Rec")
  rec@LFR <- 20.32
  rec@LR5 <- rec@LFR * 0.95
  if (plot) 
    size_lim_plot(x, Data, rec)
  rec
}
# Assign MP class
class(matlenlim8) <- "MP"

matlenlim9 <- function (x, Data, reps, plot = FALSE) {
  rec <- new("Rec")
  rec@LFR <- 22.86
  rec@LR5 <- rec@LFR * 0.95
  if (plot) 
    size_lim_plot(x, Data, rec)
  rec
}
# Assign MP class
class(matlenlim9) <- "MP"

matlenlim10 <- function (x, Data, reps, plot = FALSE) {
  rec <- new("Rec")
  rec@LFR <- 25.4
  rec@LR5 <- rec@LFR * 0.95
  if (plot) 
    size_lim_plot(x, Data, rec)
  rec
}
# Assign MP class
class(matlenlim10) <- "MP"

# Define LstepCC0 Function
LstepCC0 <- function (x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.05, llim = c(0.96, 0.98, 1.05)) { 
  ylast <- (Data@LHYear - Data@Year[1]) + 1 
  ind <- c((ylast - 15): (ylast - 6)) # Reference period for CATCH (1999-2008)
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 5 years 
  C_dat <- Data@Cat[x, ind] 
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat / (yrsmth^0.5)) 
  } else { 
    TACstar <- rep(Data@MPrec[x], reps) 
  } 
  step <- stepsz * TACstar 
  Lrecent <- mean(Data@ML[ind2], na.rm = TRUE) # Mean of recent ML 
  Lave <- mean(Data@ML[ind]) # Mean of average ML - 1999-2008 
  rat <- Lrecent / Lave 
  if (rat < llim[1]) {
    TAC <- TACstar - 2 * step 
  } else if (rat < llim[2]) { 
    TAC <- TACstar - step 
  } else if (rat > llim[3]) { 
    TAC <- TACstar + step 
  } else { 
    TAC <- TACstar 
  } 
  TACfilter(TAC) 
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
# Assign MP class
class(LstepCC0) <- "MP"

# Define Itarget0.5_0.7_1.0 Function (Near MSY)
Itarget0.5_0.7_1.0 <- function (x, Data, reps = 100, yrsmth = 5, I0 = 0.7, xx = 0, Imulti = 1.0, w = 0.5) {
  dependencies = "Data@Cat, Data@CV_Cat"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- c((ylast - 15): (ylast - 6)) # Reference period for CATCH: 1999-2008
  ind2 <- ((ylast - (yrsmth - 1)):ylast) # Period for index: last 5 years
  C_dat <- Data@Cat[x, ind]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat, na.rm = TRUE), Data@CV_Cat / (yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind2], na.rm = TRUE) # Mean of recent Index
  Iave <- mean(Data@Ind[x, ind], na.rm = TRUE) # Mean of average Index - 1999-2008
  Itarget <- Iave * Imulti
  I0 <- I0 * Iave
  if (Irecent >= I0) {
    TAC <- TACstar * (w + (1 - w) * ((Irecent - I0) / (Itarget - I0)))
  } else {
    TAC <- w * TACstar * (Irecent / I0)^2
  }
  TACfilter(TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}
# Assign MP class
class(Itarget0.5_0.7_1.0) <- "MP"

# New Islope to try
Islope_0.4_10yr <- function(x, Data, reps = 100, yrsmth = 10, lambda = 0.4, xx = 0) {
  dependencies = "Data@Year, Data@Cat, Data@CV_Cat, Data@Ind"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 10 years
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat / (yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 10 years
  yind <- 1:yrsmth
  slppar <- summary(lm(I_hist ~ yind))$coefficients[2, 1:2]
  Islp <- rnorm(reps, slppar[1], slppar[2])
  TAC <- TACstar * (1 + lambda * Islp)
  TACfilter(TAC)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_0.4_10yr) <- "MP"

# MORE MPS
# Near MSY (Itarget = Iref), I0 = 0.7
Itarget0.5_0.7_1.0_8 <- function(x, Data, reps = 100, yrsmth = 5, I0 = 0.7, xx = 0, Imulti = 1.0, w = 0.5) {
  dependencies = "Data@Cat, Data@CV_Cat"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- c((ylast - 15):(ylast - 6)) # Reference period for CATCH: 1999-2008
  ind2 <- ((ylast - (yrsmth - 1)):ylast) # Period for index: last 5 years
  C_dat <- Data@Cat[x, ind]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat, na.rm = T), Data@CV_Cat / (yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind2], na.rm = T) # Mean of recent Index
  Iave <- mean(Data@Ind[x, ind], na.rm = T) # Mean of average Index - 1999-2008
  Itarget <- Iave * Imulti
  I0 <- I0 * Iave
  if (Irecent >= I0) {
    TAC <- TACstar * (w + (1 - w) * ((Irecent - I0) / (Itarget - I0)))
  } else {
    TAC <- w * TACstar * (Irecent / I0)^2
  }
  TACfilter(TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 20.3 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Itarget0.5_0.7_1.0_8) <- "MP"

# Near MSY (Itarget = Iref), I0 = 0.7
Itarget0.5_0.7_1.0_9 <- function(x, Data, reps = 100, yrsmth = 5, I0 = 0.7, xx = 0, Imulti = 1.0, w = 0.5) {
  dependencies = "Data@Cat, Data@CV_Cat"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- c((ylast - 15):(ylast - 6)) # Reference period for CATCH: 1999-2008
  ind2 <- ((ylast - (yrsmth - 1)):ylast) # Period for index: last 5 years
  C_dat <- Data@Cat[x, ind]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat, na.rm = T), Data@CV_Cat / (yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind2], na.rm = T) # Mean of recent Index
  Iave <- mean(Data@Ind[x, ind], na.rm = T) # Mean of average Index - 1999-2008
  Itarget <- Iave * Imulti
  I0 <- I0 * Iave
  if (Irecent >= I0) {
    TAC <- TACstar * (w + (1 - w) * ((Irecent - I0) / (Itarget - I0)))
  } else {
    TAC <- w * TACstar * (Irecent / I0)^2
  }
  TACfilter(TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 22.8 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Itarget0.5_0.7_1.0_9) <- "MP"

# Near MSY (Itarget = Iref), I0 = 0.7
Itarget0.5_0.7_1.0_10 <- function(x, Data, reps = 100, yrsmth = 5, I0 = 0.7, xx = 0, Imulti = 1.0, w = 0.5) {
  dependencies = "Data@Cat, Data@CV_Cat"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- c((ylast - 15):(ylast - 6)) # Reference period for CATCH: 1999-2008
  ind2 <- ((ylast - (yrsmth - 1)):ylast) # Period for index: last 5 years
  C_dat <- Data@Cat[x, ind]
  TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat, na.rm = T), Data@CV_Cat / (yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind2], na.rm = T) # Mean of recent Index
  Iave <- mean(Data@Ind[x, ind], na.rm = T) # Mean of average Index - 1999-2008
  Itarget <- Iave * Imulti
  I0 <- I0 * Iave
  if (Irecent >= I0) {
    TAC <- TACstar * (w + (1 - w) * ((Irecent - I0) / (Itarget - I0)))
  } else {
    TAC <- w * TACstar * (Irecent / I0)^2
  }
  TACfilter(TAC)
  
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 25.4 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Itarget0.5_0.7_1.0_10) <- "MP"

# Near MSY (Islope) with added retention size limit of 20.3
Islope_0.4_10yr_8 <- function(x, Data, reps = 100, yrsmth = 10, lambda = 0.4, xx = 0) {
  dependencies = "Data@Year, Data@Cat, Data@CV_Cat, Data@Ind"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 10 years
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat / (yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 10 years
  yind <- 1:yrsmth
  slppar <- summary(lm(I_hist ~ yind))$coefficients[2, 1:2]
  Islp <- rnorm(reps, slppar[1], slppar[2])
  TAC <- TACstar * (1 + lambda * Islp)
  TACfilter(TAC)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 20.3 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_0.4_10yr_8) <- "MP"

# This function defines the management procedure 'Islope_0.4_10yr_9'
Islope_0.4_10yr_9 <- function(x, Data, reps = 100, yrsmth = 10, lambda = 0.4, xx = 0) {
  dependencies = "Data@Year, Data@Cat, Data@CV_Cat, Data@Ind"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 10 years
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat / (yrsmth ^ 0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 10 years
  yind <- 1:yrsmth
  slppar <- summary(lm(I_hist ~ yind))$coefficients[2, 1:2]
  Islp <- rnorm(reps, slppar[1], slppar[2])
  TAC <- TACstar * (1 + lambda * Islp)
  TACfilter(TAC)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 22.8 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_0.4_10yr_9) <- "MP"

# This function defines the management procedure 'Islope_0.4_10yr_10'
Islope_0.4_10yr_10 <- function(x, Data, reps = 100, yrsmth = 10, lambda = 0.4, xx = 0) {
  dependencies = "Data@Year, Data@Cat, Data@CV_Cat, Data@Ind"
  ylast <- (Data@LHYear - Data@Year[1]) + 1
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 10 years
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * trlnorm(reps, mean(C_dat), Data@CV_Cat / (yrsmth ^ 0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 10 years
  yind <- 1:yrsmth
  slppar <- summary(lm(I_hist ~ yind))$coefficients[2, 1:2]
  Islp <- rnorm(reps, slppar[1], slppar[2])
  TAC <- TACstar * (1 + lambda * Islp)
  TACfilter(TAC)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 25.4 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_0.4_10yr_10) <- "MP"

# Running NOAA OM with MPS
FinalMPs <- c("matlenlim", "matlenlim8", "matlenlim9", "matlenlim10",
              "Islope_0.4_10yr", "Itarget0.5_0.7_1.0", "LstepCC0",
              "Islope_0.4_10yr_8", "Islope_0.4_10yr_9", "Islope_0.4_10yr_10",
              "Itarget0.5_0.7_1.0_8", "Itarget0.5_0.7_1.0_9", "Itarget0.5_0.7_1.0_10")

# MSE-2
NOAAOM@nsim <- 500 # Increase/Decrease simulations and try
MSE2 <- runMSE(OM = NOAAOM, MPs = FinalMPs)
Converge(MSE2)
summary(MSE2)
Tplot(MSE2, lab.size = 1 / 15)

# Get results and plot tradeoffs
Results2 <- summary(MSE2)
Results2
NOAA_plot2(MSE2)
plot(MSE2)

slotNames(MSE2)
boxplot(MSE2@Catch[, 2, ], xlab = "Year", ylab = "Catch")

# Get customized results (e.g., different year intervals)
Results2.1 <- cbind(Results2, NOAA_plot2(MSE2))
Results2.1
write.csv(Results2, "MSE2_Results.csv")
save("NOAAOM", "MSE2", "Results2", file = "MSE2.RData")

# MSE-3
OMA@nsim <- 1400
MSE3 <- runMSE(OM = OMA, MPs = FinalMPs)

Converge(MSE3)
summary(MSE3)
Tplot(MSE3)
dev.off()

write.csv(Results3, "MSE3_Results6_23.csv")
save("OMA", "MSE3", "Results3", file = "MSE3_6_23.RData")

P100(MSE3)
Pplot(MSE3)
# Get results and plot tradeoffs
Results3 <- summary(MSE3)
Results3
NOAA_plot(MSE3)

# Get customized results (e.g., different year intervals)
Results3.1 <- cbind(Results3, NOAA_plot2(MSE3))
Results3.1

NOAA_plot2(MSE3)

NOAAMSE2 <- runMSE(OM = NOAAMSE, MPs = MyMPs2)
LSOMMSE2 <- runMSE(OM = LSOM, MPs = NewMPs)
NewNOAAMSE <- runMSE(OM = NOAAMSE, MPs = NewMPs)

MyMPs3 <- c('matlenlim', "matlenlim8", "matlenlim9", "matlenlim10", "Itarget_SEDAR49",
            "Itarget_SEDAR49_8", "Itarget_SEDAR49_9", "Itarget_SEDAR49_10", "Islope_SEDAR49",
            "Islope_SEDAR49_8", "Islope_SEDAR49_9", "Islope_SEDAR49_10")
NOAAMSE@nsim <- 700 # Increase/Decrease simulations and try
NOAAMSE3 <- runMSE(OM = NOAAMSE, MPs = MyMPs3)
Converge(LSOMMSE2)
summary(LSOMMSE2)
Pplot(MSE2)
Pplot2(MSE2)

LSOM@nsim <- 700
LSOMMSE3 <- runMSE(OM = LSOM, MPs = MyMPs3)
Converge(LSOMMSE3)

# Another Run with Lstep and without Itarget
MyMPs4 <- c('matlenlim', "matlenlim8", "matlenlim9", "matlenlim10", "LstepCC0", "Islope_SEDAR49",
            "Islope_SEDAR49_8", "Islope_SEDAR49_9", "Islope_SEDAR49_10")
NOAAMSE@nsim <- 600
NOAAMSE4 <- runMSE(OM = NOAAMSE, MPs = MyMPs4)
Converge(NOAAMSE4)

LSOM@nsim <- 800
LSOMMSE4 <- runMSE(OM = LSOM, MPs = MyMPs4)
Converge(LSOMMSE4)

# Calculate the convergence statistics
Converge(NOAAMSE2)
Converge(LSOMMSE2)
Converge(NewNOAAMSE)
?Converge

# The summary function can be used to generate a table of MP performance with respect to a set of performance metrics:
summary(NOAAMSE2)
summary(LSOMMSE2)
summary(NOAAMSE3)
summary(LSOMMSE3)
summary(NewNOAAMSE)
Pplot(NewNOAAMSE)

# Summary statistics for MSE without Itarget and with LstepCC0
summary(NOAAMSE4)
STY(NOAAMSE4)
AAVY(NOAAMSE4, Ref = 0.15)
NOAA_plot(NOAAMSE4)
Pplot2(NOAAMSE4, traj = "quant", quants = c(0.2, 0.8))
PWhisker(NOAAMSE4)

summary(LSOMMSE4)
STY(MSE3)
P10(MSE3)
P50(MSE3)
P100(MSE3)
AAVY(LSOMMSE4, Ref = 0.15)
NOAA_plot(LSOMMSE4)
Pplot2(LSOMMSE4, traj = "quant", quants = c(0.2, 0.8))

# The Tplot function includes minimum acceptable risk thresholds indicated by the horizontal and vertical gray shading.
# MPs that fail to meet one or both of the risk thresholds for each axis are shown in
# italics text. The Tplot function returns a data frame showing the performance
# of each MP with respect to the 5 performance metrics, and whether the MP is
# Satisficed, i.e., if it meets the minimum performance criteria for all performance metrics.
Tplot(NOAAMSE2)
Tplot(LSOMMSE3, lab.size = NULL)
Tplot(NOAAMSE3, lab.size = NULL)

Tplot(NOAAMSE4)
Tplot(LSOMMSE4)

# The Pplot function plots the trajectories of biomass, fishing mortality, and
# relative yield for the Management Procedures.
Pplot(NOAAMSE2)
Pplot(LSOMMSE3)
Pplot(NOAAMSE3)

# Kobe plots are often used in stock assessment and MSE to examine the proportion
# of time the stock spends in different states.
Kplot(NOAAMSE2)
Kplot(LSOMMSE3)
Kplot(NOAAMSE3)

# Subset by Performance
statsNOAA = summary(NOAAMSE3)
statsLS = summary(LSOMMSE3)

acceptNOAA = which(statsNOAA$P50 > .70) # Greater than 70% probability that biomass is above 0.5Bmsy
MPsNOAA = statsNOAA[acceptNOAA, "MP"] # Acceptable MPs
subMSENOAA = Sub(NOAAMSE3, MPs = MPsNOAA)
Tplot(subMSENOAA)

acceptLS = which(statsLS$P50 > .70) # Greater than 70% probability that biomass is above 0.5Bmsy
MPsLS = statsLS[acceptLS, "MP"] # Acceptable MPs
subMSELS = Sub(LSOMMSE3, MPs = MPsLS)
Tplot(subMSELS)

## LH2OM Test

# Import SEDAR OM
?str
NOAAOM <- XL2OM("Lane_Snapper_GOM_NOAA")
?LH2OM
# MPs
NewMPs <- c("Islope_0.4_10yr", "Itarget0.5_0.7_1.0", "LstepCC0")
str(NOAALH2OM@cpars)
str(LS_OM@cpars)
NOAAOM@Linf

NOAAOM@nsim
NOAAOM@M <- c(0.33, 0.366)
NOAAOM@K <- c(0.116, 0.219)
NOAAOM@Linf <- c(42.2, 49.3)
NOAAOM@L50 <- c(23.5, 24.5)
NOAAOM <- LH2OM(NOAAOM)

NewNOAAMSE <- runMSE(OM = NOAAOM, MPs = NewMPs)
Converge(NewNOAAMSE)
summary(NewNOAAMSE)

# Adds 8 in size limit to I target
Itarget_SEDAR49_8 <- function (x, Data, reps = 100, plot = FALSE, yrsmth = 5, xx = 0, 
                               Imulti = 1.5) 
{
  ylast <- (Data@LHYear[1] - Data@Year[1]) + 1   # BB: moved up one line
  ind <- c((ylast - 14): (ylast - 6)) # Reference period for CATCH (2000-2008)
  ind2 <- ((ylast - (yrsmth - 1)): ylast) # Period for index: last 5 years
  ind3 <- ((ylast - (yrsmth * 2 - 1)): ylast)
  C_dat <- Data@Cat[x, ind]
  TACstar <- (1 - xx) * MSEtool::trlnorm(reps, mean(C_dat, 
                                                    na.rm = TRUE), Data@CV_Cat[x, 1]/(yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind2], na.rm = TRUE) # Mean of recent Index 
  Iave <- mean(Data@Ind[x, ind], na.rm = TRUE) # Mean of average Index - 2000-2008
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    TAC <- TACstar * (1 + ((Irecent - I0)/(Itarget - I0)))
  } else {
    TAC <- TACstar * (Irecent/I0)^2
  }
  TAC <- MSEtool::TACfilter(TAC)
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(op)
    par(mfrow = c(1, 2))
    ylim <- range(c(Data@Ind[x, ], Itarget, I0))
    plot(Data@Year, Data@Ind[x, ], type = "l", lwd = 2, 
         bty = "l", xlab = "Year", ylab = "Index", 
         ylim = ylim)
    points(max(Data@Year), mean(Data@Ind[x, ind], na.rm = TRUE), 
           cex = 2, pch = 16, col = "blue")
    text(max(Data@Year), mean(Data@Ind[x, ind], na.rm = TRUE), 
         cex = 1, "Irecent", pos = 3, col = "blue", 
         xpd = NA)
    lines(Data@Year[ind3], rep(mean(Data@Ind[x, ind3], na.rm = TRUE), 
                               length(ind3)), lty = 2, col = "orange")
    text(mean(Data@Year[ind3]), mean(Data@Ind[x, ind3], na.rm = TRUE), 
         "Iave", col = "orange", pos = 1)
    points(max(Data@Year), Itarget, cex = 2, pch = 16, col = "green")
    text(max(Data@Year), Itarget, cex = 1, "Itarget", 
         pos = 3, col = "green", xpd = NA)
    points(max(Data@Year), I0, cex = 2, pch = 16, col = "red")
    text(max(Data@Year), I0, cex = 1, "I0", pos = 3, 
         col = "red", xpd = NA)
    ylim <- range(c(C_dat, TACstar, TAC))
    Years <- Data@Year[ind2]
    if (max(Years) != max(Data@Year)) {
      Years <- c(Years, (max(Data@Year[ind2]) + 1): max(Data@Year))
    }
    yrs <- length(Years) - length(C_dat)
    Cdat <- c(C_dat, rep(NA, yrs))
    plot(c(Years, max(Years) + 1), c(Cdat, NA), type = "l", 
         lwd = 2, bty = "l", xlab = "Year", ylab = paste0("Catch (", 
                                                          Data@Units, ")"), ylim = ylim)
    abline(v = max(Data@Year[ind2]), col = "gray", 
           lty = 3)
    points(max(Data@Year[ind2]), mean(TACstar, na.rm = TRUE), 
           cex = 2, col = "orange", pch = 16)
    text(max(Data@Year[ind2]), mean(TACstar, na.rm = TRUE), 
         "mean TAC*", pos = 2, xpd = NA, col = "orange")
    boxplot(TAC, at = max(Years) + 1, add = TRUE, col = "gray", 
            axes = FALSE)
  }
  list(TAC = TAC)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 20.32 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Itarget_SEDAR49_8) <- "MP"

# Adds 9 in size limit to I target
Itarget_SEDAR49_9 <- function (x, Data, reps = 100, plot = FALSE, yrsmth = 5, xx = 0, 
                               Imulti = 1.5) 
{
  ylast <- (Data@LHYear[1] - Data@Year[1]) + 1   # BB: moved up one line
  ind <- c((ylast - 14): (ylast - 6)) # Reference period for CATCH (2000-2008)
  ind2 <- ((ylast - (yrsmth - 1)): ylast) # Period for index: last 5 years
  ind3 <- ((ylast - (yrsmth * 2 - 1)): ylast)
  C_dat <- Data@Cat[x, ind]
  TACstar <- (1 - xx) * MSEtool::trlnorm(reps, mean(C_dat, 
                                                    na.rm = TRUE), Data@CV_Cat[x, 1]/(yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind2], na.rm = TRUE) # Mean of recent Index 
  Iave <- mean(Data@Ind[x, ind], na.rm = TRUE) # Mean of average Index - 2000-2008
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    TAC <- TACstar * (1 + ((Irecent - I0)/(Itarget - I0)))
  }
  else {
    TAC <- TACstar * (Irecent/I0)^2
  }
  TAC <- MSEtool::TACfilter(TAC)
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(op)
    par(mfrow = c(1, 2))
    ylim <- range(c(Data@Ind[x, ], Itarget, I0))
    plot(Data@Year, Data@Ind[x, ], type = "l", lwd = 2, 
         bty = "l", xlab = "Year", ylab = "Index", 
         ylim = ylim)
    points(max(Data@Year), mean(Data@Ind[x, ind], na.rm = TRUE), 
           cex = 2, pch = 16, col = "blue")
    text(max(Data@Year), mean(Data@Ind[x, ind], na.rm = TRUE), 
         cex = 1, "Irecent", pos = 3, col = "blue", 
         xpd = NA)
    lines(Data@Year[ind3], rep(mean(Data@Ind[x, ind3], na.rm = TRUE), 
                               length(ind3)), lty = 2, col = "orange")
    text(mean(Data@Year[ind3]), mean(Data@Ind[x, ind3], na.rm = TRUE), 
         "Iave", col = "orange", pos = 1)
    points(max(Data@Year), Itarget, cex = 2, pch = 16, col = "green")
    text(max(Data@Year), Itarget, cex = 1, "Itarget", 
         pos = 3, col = "green", xpd = NA)
    points(max(Data@Year), I0, cex = 2, pch = 16, col = "red")
    text(max(Data@Year), I0, cex = 1, "I0", pos = 3, 
         col = "red", xpd = NA)
    ylim <- range(c(C_dat, TACstar, TAC))
    Years <- Data@Year[ind2]
    if (max(Years) != max(Data@Year)) {
      Years <- c(Years, (max(Data@Year[ind2]) + 1): max(Data@Year))
    }
    yrs <- length(Years) - length(C_dat)
    Cdat <- c(C_dat, rep(NA, yrs))
    plot(c(Years, max(Years) + 1), c(Cdat, NA), type = "l", 
         lwd = 2, bty = "l", xlab = "Year", ylab = paste0("Catch (", 
                                                          Data@Units, ")"), ylim = ylim)
    abline(v = max(Data@Year[ind2]), col = "gray", 
           lty = 3)
    points(max(Data@Year[ind2]), mean(TACstar, na.rm = TRUE), 
           cex = 2, col = "orange", pch = 16)
    text(max(Data@Year[ind2]), mean(TACstar, na.rm = TRUE), 
         "mean TAC*", pos = 2, xpd = NA, col = "orange")
    boxplot(TAC, at = max(Years) + 1, add = TRUE, col = "gray", 
            axes = FALSE)
  }
  list(TAC = TAC)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 22.86 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Itarget_SEDAR49_9) <- "MP"

# Adds 10 in size limit to I target
Itarget_SEDAR49_10 <- function (x, Data, reps = 100, plot = FALSE, yrsmth = 5, xx = 0, 
                                Imulti = 1.5) 
{
  ylast <- (Data@LHYear[1] - Data@Year[1]) + 1   # BB: moved up one line
  ind <- c((ylast - 14): (ylast - 6)) # Reference period for CATCH (2000-2008)
  ind2 <- ((ylast - (yrsmth - 1)): ylast) # Period for index: last 5 years
  ind3 <- ((ylast - (yrsmth * 2 - 1)): ylast)
  C_dat <- Data@Cat[x, ind]
  TACstar <- (1 - xx) * MSEtool::trlnorm(reps, mean(C_dat, 
                                                    na.rm = TRUE), Data@CV_Cat[x, 1]/(yrsmth^0.5))
  Irecent <- mean(Data@Ind[x, ind2], na.rm = TRUE) # Mean of recent Index 
  Iave <- mean(Data@Ind[x, ind], na.rm = TRUE) # Mean of average Index - 2000-2008
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    TAC <- TACstar * (1 + ((Irecent - I0)/(Itarget - I0)))
  }
  else {
    TAC <- TACstar * (Irecent/I0)^2
  }
  TAC <- MSEtool::TACfilter(TAC)
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(op)
    par(mfrow = c(1, 2))
    ylim <- range(c(Data@Ind[x, ], Itarget, I0))
    plot(Data@Year, Data@Ind[x, ], type = "l", lwd = 2, 
         bty = "l", xlab = "Year", ylab = "Index", 
         ylim = ylim)
    points(max(Data@Year), mean(Data@Ind[x, ind], na.rm = TRUE), 
           cex = 2, pch = 16, col = "blue")
    text(max(Data@Year), mean(Data@Ind[x, ind], na.rm = TRUE), 
         cex = 1, "Irecent", pos = 3, col = "blue", 
         xpd = NA)
    lines(Data@Year[ind3], rep(mean(Data@Ind[x, ind3], na.rm = TRUE), 
                               length(ind3)), lty = 2, col = "orange")
    text(mean(Data@Year[ind3]), mean(Data@Ind[x, ind3], na.rm = TRUE), 
         "Iave", col = "orange", pos = 1)
    points(max(Data@Year), Itarget, cex = 2, pch = 16, col = "green")
    text(max(Data@Year), Itarget, cex = 1, "Itarget", 
         pos = 3, col = "green", xpd = NA)
    points(max(Data@Year), I0, cex = 2, pch = 16, col = "red")
    text(max(Data@Year), I0, cex = 1, "I0", pos = 3, 
         col = "red", xpd = NA)
    ylim <- range(c(C_dat, TACstar, TAC))
    Years <- Data@Year[ind2]
    if (max(Years) != max(Data@Year)) {
      Years <- c(Years, (max(Data@Year[ind2]) + 1): max(Data@Year))
    }
    yrs <- length(Years) - length(C_dat)
    Cdat <- c(C_dat, rep(NA, yrs))
    plot(c(Years, max(Years) + 1), c(Cdat, NA), type = "l", 
         lwd = 2, bty = "l", xlab = "Year", ylab = paste0("Catch (", 
                                                          Data@Units, ")"), ylim = ylim)
    abline(v = max(Data@Year[ind2]), col = "gray", 
           lty = 3)
    points(max(Data@Year[ind2]), mean(TACstar, na.rm = TRUE), 
           cex = 2, col = "orange", pch = 16)
    text(max(Data@Year[ind2]), mean(TACstar, na.rm = TRUE), 
         "mean TAC*", pos = 2, xpd = NA, col = "orange")
    boxplot(TAC, at = max(Years) + 1, add = TRUE, col = "gray", 
            axes = FALSE)
  }
  list(TAC = TAC)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 25.4 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Itarget_SEDAR49_10) <- "MP"

# Slope based on SEDAR49 data
Islope_SEDAR49 <-function (x, Data, reps = 100, yrsmth = 5, lambda = 0.4, xx = 0.0) 
{
  ylast <- (Data@LHYear - Data@Year[1]) + 1   # BB
  ind <- c((ylast - 14): (ylast - 6)) # Reference period for CATCH (2000-2008) 
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 5 years
  Years <- Data@Year[ind]
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * MSEtool::trlnorm(reps, mean(C_dat, 
                                                      na.rm = TRUE), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 5 years
  yind <- 1:yrsmth
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  if (reps > 1) {
    Islp <- rnorm(reps, slppar[1], slppar[2])
  } else {
    Islp <- slppar[1]
  }
  TAC <- TACstar * (1 + lambda * Islp)
  list(TAC = TAC, TACstar = TACstar, I_hist = I_hist, Islp = Islp, 
       C_dat = C_dat, Data = Data, Years = Years)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_SEDAR49) <- "MP"

# Slope with 8in size limit
Islope_SEDAR49_8 <-function (x, Data, reps = 100, yrsmth = 5, lambda = 0.4, xx = 0.0) 
{
  ylast <- (Data@LHYear - Data@Year[1]) + 1   # BB
  ind <- c((ylast - 14): (ylast - 6)) # Reference period for CATCH (2000-2008) 
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 5 years
  Years <- Data@Year[ind]
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * MSEtool::trlnorm(reps, mean(C_dat, 
                                                      na.rm = TRUE), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 5 years
  yind <- 1:yrsmth
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  if (reps > 1) {
    Islp <- rnorm(reps, slppar[1], slppar[2])
  } else {
    Islp <- slppar[1]
  }
  TAC <- TACstar * (1 + lambda * Islp)
  list(TAC = TAC, TACstar = TACstar, I_hist = I_hist, Islp = Islp, 
       C_dat = C_dat, Data = Data, Years = Years)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 20.32 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_SEDAR49_8) <- "MP"

# Slope with 9in size limit
Islope_SEDAR49_9 <-function (x, Data, reps = 100, yrsmth = 5, lambda = 0.4, xx = 0.0) 
{
  ylast <- (Data@LHYear - Data@Year[1]) + 1   # BB
  ind <- c((ylast - 14): (ylast - 6)) # Reference period for CATCH (2000-2008) 
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 5 years
  Years <- Data@Year[ind]
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * MSEtool::trlnorm(reps, mean(C_dat, 
                                                      na.rm = TRUE), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 5 years
  yind <- 1:yrsmth
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  if (reps > 1) {
    Islp <- rnorm(reps, slppar[1], slppar[2])
  } else {
    Islp <- slppar[1]
  }
  TAC <- TACstar * (1 + lambda * Islp)
  list(TAC = TAC, TACstar = TACstar, I_hist = I_hist, Islp = Islp, 
       C_dat = C_dat, Data = Data, Years = Years)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 22.86 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_SEDAR49_9) <- "MP"

# Slope with 10in size limit
Islope_SEDAR49_10 <-function (x, Data, reps = 100, yrsmth = 5, lambda = 0.4, xx = 0.0) 
{
  ylast <- (Data@LHYear - Data@Year[1]) + 1   # BB
  ind <- c((ylast - 14): (ylast - 6)) # Reference period for CATCH (2000-2008) 
  ind2 <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year) # Period for index: last 5 years
  Years <- Data@Year[ind]
  C_dat <- Data@Cat[x, ind]
  if (is.na(Data@MPrec[x]) || length(Data@Year) == ylast + 1) {
    TACstar <- (1 - xx) * MSEtool::trlnorm(reps, mean(C_dat, 
                                                      na.rm = TRUE), Data@CV_Cat/(yrsmth^0.5))
  } else {
    TACstar <- rep(Data@MPrec[x], reps)
  }
  I_hist <- Data@Ind[x, ind2] # Period for index: last 5 years
  yind <- 1:yrsmth
  slppar <- summary(lm(log(I_hist) ~ yind))$coefficients[2, 1:2]
  if (reps > 1) {
    Islp <- rnorm(reps, slppar[1], slppar[2])
  } else {
    Islp <- slppar[1]
  }
  TAC <- TACstar * (1 + lambda * Islp)
  list(TAC = TAC, TACstar = TACstar, I_hist = I_hist, Islp = Islp, 
       C_dat = C_dat, Data = Data, Years = Years)
  Rec <- new("Rec")
  Rec@TAC <- TAC
  Rec@LFR <- 25.4 # Adds length at full retention
  Rec@LR5 <- Rec@LFR * 0.95
  Rec
}

# Assign it the class 'MP' so that DLMtool recognizes the function as a management procedure
class(Islope_SEDAR49_10) <- "MP"
