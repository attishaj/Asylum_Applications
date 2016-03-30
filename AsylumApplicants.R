########## Asylum Applicants #########

##### Libraries ####

#load libraries
library(tseries)
library(fBasics)
library(forecast)
library(lmtest)
library(fUnitRoots)

##### Data Load and Exploratory Analysis for Syria, Afghanistan, Iraq to EU ####
data=read.table("SyrAfgIrq2EU.csv",header=T, sep=',') 
head(data)
rate = data[,1]
head(rate)
ratets = ts(rate, start=c(1999,1), freq=12)
ratets

#time plot of original data
plot(ratets, xlab="Time", ylab="# of Applications", lwd=2,
     main="Asylum Applications Jan 1999 to Oct 2015\nOrigin: Syria/Afg/Iraq    Destination: EU")

#take log
lratets=log(ratets)
plot(lratets, xlab="Time", ylab="Log (# of Applications)", lwd=2,
     main="Asylum Applications Jan 1999 to Oct 2015\nOrigin: Syria/Afg/Iraq    Destination: EU")


##### Data Load and Exploratory Analysis for Syria to Germany ####

myd <- read.table("Syr2Ger.csv",header=T, sep=',') 
head(myd)
rate <- myd$Value
head(rate)
ratets <- ts(rate, start=c(1999,1), freq=12)
ratets

#time plot of original data
plot(ratets, type="l", xlab="Time", ylab="# of Applications", lwd=2,
     main="Asylum Applications Jan 1999 to Oct 2015\nOrigin: Syria     Destination: Germany")

#take the log and plot time series and ACF
lratets = log(ratets)

plot(lratets, xlab="Time", ylab="Log (# of Applications)", lwd=2,
     main="Asylum Applications Jan 1999 to Oct 2015\nOrigin: Syria     Destination: Germany")



### Proceed with Syria to Germany data for remainder of code

a <- acf(coredata(lratets), plot = F)
plot(a, main="ACF of Log Data: Syria to Germany", lwd=2)

#dickey fuller test
adfTest(lratets, lags=3, type=c("ct"))
adfTest(lratets, lags=5, type=c("ct"))
adfTest(lratets, lags=7, type=c("ct"))
#all p-values were close to 1. time process is unit-root non-stationary


#let's apply the first difference and then re-run the plotting and
#DF tests to see if we can now reject the H0.
dlratets=diff(lratets)

plot(dlratets, type="l", xlab="Time", ylab="Log (# of Applications)", 
     main="First Differenced Log Data\nOrigin: Syria     Destination: Germany")
abline(h=mean(dlratets), col="green4", lty=5, lwd=2)

b <- acf(as.vector(dlratets),lag.max=30, plot=F)
plot(b, main="ACF of First Differenced Log Data: Syria to Germany", lwd=2)

c <- pacf(as.vector(dlratets),lag.max=30, plot=F)
plot(c, main="PACF of First Differenced Log Data: Syria to Germany", lwd=2)

adfTest(dlratets, lags=3, type=c("c"))
adfTest(dlratets, lags=5, type=c("c"))
adfTest(dlratets, lags=7, type=c("c"))
#can reject H0 now.


#basic statistics and normality testing for differenced log series
basicStats(dlratets)

hist(dlratets, xlab="Difference of Log Rate of Applications", prob=TRUE, main="Histogram")
xfit<-seq(min(dlratets),max(dlratets),length=40)
yfit<-dnorm(xfit,mean=mean(dlratets),sd=sd(dlratets))
lines(xfit, yfit, col="blue", lwd=2)

qqnorm(dlratets)
qqline(dlratets, col = 2) 

normalTest(dlratets,method=c("jb"))  


##### Model Creation ####

#model m1
m1=auto.arima(lratets, trace=T, ic="bic", allowdrift = T, allowmean = T)
coeftest(m1)
m1
#we get an ARIMA(2,2,2)(0,0,1)[12] model
#both AR coefficients are not significant.


#try removing both non-significant AR coefficients
#try fitting ARIMA(0,2,2)(0,0,1)[12] model
#model m2
m2=Arima(x=lratets, order=c(0,2,2), seasonal=c(0,0,1), method="ML")
coeftest(m2)
m2
#all coefficients are now significant
#bic value (81.03) was lower than that of model m1 (91.63)

### SIDEBAR ###
#since the function is suggesting a 2nd difference, let's go back and 
#check the ACF and plot for the second difference

#take the second difference of the logged rate and plot the ACF values
d2lratets <- diff(dlratets)

plot(d2lratets, type="l", xlab="Time", ylab="Log(# Applications)", 
     main="Second Differenced Log Data\nOrigin: Syria     Destination: Germany")
abline(h=mean(d2lratets), col="green4", lty=5, lwd=2)

b2 <- acf(as.vector(d2lratets),lag.max=30, plot=F)
plot(b2, main="ACF of Second Differenced Log Data: Syria to Germany", lwd=2)
#looks stationary

#basic statistics and normality testing for second differenced log series
basicStats(d2lratets)

hist(d2lratets, xlab="Second Difference of Log Rate of Applications", prob=TRUE, main="Histogram")
xfit<-seq(min(d2lratets),max(d2lratets),length=40)
yfit<-dnorm(xfit,mean=mean(dlratets),sd=sd(d2lratets))
lines(xfit, yfit, col="blue", lwd=2)

qqnorm(d2lratets)
qqline(d2lratets, col = 2) 

normalTest(d2lratets,method=c("jb"))  

### END SIDEBAR ###


#it was interesting that auto.arima decided to take second difference,
#and not just the first. what if we restrict it to just the first difference?
#(dickey fuller test on first differenced data rejected H0)

#model m3
m3=auto.arima(lratets, d=1, trace=T, ic="bic", allowdrift = T, allowmean = T)
coeftest(m3)
m3
#we get an ARIMA(0,1,1)(0,0,1)[12] model
#both parameters are significant


##### Residual Analysis and Forecasting ####

### Model m1 ###
#diagnostics

#analyze the residuals
acf(as.vector(m1$resid))
#acf plot shows residual correlations are zero with borderline case at lag 10

Box.test(m1$residuals, 7, "Ljung-Box", fitdf=5)
Box.test(m1$residuals, 9, "Ljung-Box", fitdf=5)
Box.test(m1$residuals, 12, "Ljung-Box", fitdf=5)
Box.test(m1$residuals, 15, "Ljung-Box", fitdf=5)
#p-values for LJB test do not confirm white noise series for all tested lags
#perhaps due to large df and small sample size

hist(m1$residuals)
qqnorm(m1$residuals)
qqline(m1$residuals, col=2)
#residuals appear to be fairly normal

#Forecast for Model m1
f1=forecast.Arima(m1, h=5)
f1exp=exp(f1$mean)
ratets
f1exp
plot(f1)
plot(f1, include=100)
plot(f1, include=50)
lines(ts(tail(f1$fitted,50), frequency=12,start=c(2011,9)), lty=5, lwd=2, col="blue")
#forecasts appear to be consistent with the overall dynamic behavior of the process

#backtesting
source("backtest.R")
#there's 202 values in the dataset, so we'll use 182 for training (~90%)
backtest(m1, lratets, h=1, orig=182)
#MAPE is 1.8%. Pretty good




### Model m2 ###
#diagnostics

#analyze the residuals
pResid <- acf(as.vector(m2$resid), plot =F)
plot(pResid, main="ACF of Model 2 Residuals", lwd=2)
#acf plot shows residual correlations are zero with borderline case at lag 10

Box.test(m2$residuals, 5, "Ljung-Box",fitdf=3) 
Box.test(m2$residuals, 7, "Ljung-Box", fitdf=3)
Box.test(m2$residuals, 9, "Ljung-Box", fitdf=3)
Box.test(m2$residuals, 12, "Ljung-Box", fitdf=3)
Box.test(m2$residuals, 15, "Ljung-Box", fitdf=3)
#p-values indicate H0 can't be rejected at 5% significance level

hist(m2$residuals, xlab = "Residuals", prob=T, main="Histogram of Model 2 Residuals")
xfit <-seq(min(m2$residuals), max(m2$residuals), length=40)
yfit <- dnorm(xfit, mean = mean(m2$residuals), sd=sd(m2$residuals))
lines(xfit, yfit, col="blue", lwd=2)

qqnorm(m2$residuals)
qqline(m2$residuals, col=2)
#residuals appear to be fairly normal


#Forecast for Model m2
f2=forecast.Arima(m2, h=5)
f2exp=exp(f2$mean)
ratets
f2exp

plot(f2, ylab="Log(# of Applications)", xlab="Time", main="Model M2 Forecasts")
lines(ts(f2$fitted, frequency=12,start=c(1999,1)), lty=5, col="blue")

plot(tail(ratets,50), type="l")
plot(f2, include=50, ylab="Log(# of Applications)", xlab="Time", main="Model M2 Forecasts")
lines(ts(tail(f2$fitted,50), frequency=12,start=c(2011,9)), lty=5, lwd=2, col="blue")
#forecasts appear to be consistent with the overall dynamic behavior of the process


#backtesting
source("backtest.R")
#there's 202 values in the dataset, so we'll use 182 for training (~90%)
backtest(m2, lratets, h=1, orig=182)
#MAPE is 1.9%. Pretty good



### Model m3 ###
#diagnostics
#analyze the residuals
pResid3 <- acf(as.vector(m3$resid), plot =F)
plot(pResid3, main="ACF of Model 3 Residuals", lwd=2)
#acf plot shows residuals are white noise

Box.test(m3$residuals, 5, "Ljung-Box",fitdf=2) 
Box.test(m3$residuals, 7, "Ljung-Box", fitdf=2)
Box.test(m3$residuals, 9, "Ljung-Box", fitdf=2)
Box.test(m3$residuals, 12, "Ljung-Box",fitdf=2) 
Box.test(m3$residuals, 15, "Ljung-Box",fitdf=2) 
#acf of the residuals exhibit white noise

hist(m3$residuals, xlab = "Residuals", prob=T, main="Histogram of Model 3 Residuals")
xfit <-seq(min(m3$residuals), max(m3$residuals), length=40)
yfit <- dnorm(xfit, mean = mean(m3$residuals), sd=sd(m3$residuals))
lines(xfit, yfit, col="blue", lwd=2)

qqnorm(m3$residuals)
qqline(m3$residuals, col=2)
#residuals appear to be fairly normal

#Forecast for Model m3
f3=forecast.Arima(m3, h=5)
f3exp=exp(f3$mean)
ratets
f3exp

plot(f3, ylab="Log(# of Applications)", xlab="Time", main="Model M3 Forecasts")
lines(ts(f3$fitted, frequency=12,start=c(1999,1)), lty=5, col="blue")

plot(f3, include=50, ylab="Log(# of Applications)", xlab="Time", main="Model M3 Forecasts")
lines(ts(tail(f3$fitted,50), frequency=12,start=c(2011,9)), lty=5, lwd=2, col="blue")
#forecasts do NOT appear to be consistent with the overall dynamic behavior of the process



#backtesting
source("backtest.R")
#there's 202 values in the dataset, so we'll use 182 for training (~90%)
backtest(m3, lratets, h=1, orig=182)
#MAPE is 2.2%, still pretty good, but forecasts are not consistent with time process
