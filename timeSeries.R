# SETTINGS ####
rm(list = ls())
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

install_packages_if_not_present <- function(x) {
    if (sum(!x %in% installed.packages())) {
        install.packages(x[!x %in% installed.packages()])
    }
}
packages <- c("fBasics", "forecast")
install_packages_if_not_present(packages)
sapply(packages, require, character.only = TRUE)

source("./FTS_functions.R")

# 1. Find at least two linear time series models, using the Box-Jenkins methodology, for the quarterly
#    earnings per share of Coca-Cola Company from the first quarter of 1983 to the third quarter of 2009.
#    Identify your models using the entire available sample (coca_cola_earnings.csv) ####

# Quarterly data -> s=4 if there is seasonal behaviour
data <- read.csv("./coca_cola_earnings.csv", header = TRUE, sep = ";", dec = ",") # 107 rows, 2 columns
y <- data[, 2]
tmp <- as.Date(paste(substr(data[,1], 1, 4), substr(data[,1], 5, 6), substr(data[,1], 7, 8), sep = "-"))

nlags <- 80     # i.e. 10-20 years seems reasonable

ts_acf_pacf(y, tmp, seriesName = "y", nlags = nlags)
# It does NOT look stationary
# Take logs to stabilise the variance (no precautions needed as y is +ve throughout)
log.y <- log(y)
ts_acf_pacf(log.y, tmp, seriesName = "log(y)", nlags = nlags)
# Now stable variance

# Seasonal parameter 
s <- 4   

ndiffs(log.y, alpha=0.05, test=c("adf")) # 1, regular differences d=1
nsdiffs(log.y, m=s, test=c("ocsb"))      # 1, seasonal differences D=1

# Apply the differences (regular and seasonal)
fit <- arima(log.y, order = c(0,1,0), seasonal = list(order = c(0,1,0), period = s)) 
ts_acf_pacf(fit$residuals, tmp, seriesName = "(0,1,0)x(0,1,0)", nlags = nlags)
# It looks stationary in the mean

# ACF has: lag 1 out of limits
#          lag 3 out of limits
#          lag 4 out of limits <-- seasonal (but not lag 8/12/16...)
#          lag 5 out of limits

# MODEL 1 ####
# Let's start with (0,1,0)x(0,1,1)_4 because of lag 4 out of limits
fit <- arima(log.y, order = c(0,1,0), seasonal = list(order = c(0,1,1), period = s)) 
check_significance(fit)
# Model is significant, check the residuals
ts_acf_pacf(fit$residuals, tmp, seriesName = "(0,1,0)x(0,1,1) residuals", nlags = nlags)
# NOT white noise!

# Let's include the regular MA
fit <- arima(log.y, order = c(0,1,1), seasonal = list(order = c(0,1,1), period = s)) 
check_significance(fit)
# Model is significant, check the residuals
ts_acf_pacf(fit$residuals, tmp, seriesName = "(0,1,1)x(0,1,1) residuals", nlags = nlags)
# Formal test (NOTE: `tsdiag(fit, gof.lag = nlags)` plots a diagnostic)
Box.test(fit$residuals, lag=18) # p-value = 0.4205
Box.test(fit$residuals, lag=9)  # p-value = 0.2367
# Residuals are white noise
# Hence, (0,1,1)x(0,1,1)_4 is a valid model
model1 <- fit


# MODEL 2 ####
# Try with regular and seasonal AR instead of the differences
# So, instead of taking differences, try (1,0,0)x(1,0,0)_4 and check
fit <- arima(log.y, order = c(1,0,0), seasonal = list(order = c(1,0,0), period = s)) 
fit # Seasonal coeff is close-ish to 1 but the regular coeff isn't
# Let's try then with AR and MA parts keeping d=1 but D=0
fit <- arima(log.y, order = c(1,1,1), seasonal = list(order = c(1,0,1), period = s)) 
check_significance(fit)
# ar1 coefficient is not significant, force it to zero
fit <- arima(log.y, order = c(1,1,1), fixed = c(0,NA,NA,NA), seasonal = list(order = c(1,0,1), period = s)) 
# Model is significant, check the residuals
ts_acf_pacf(fit$residuals, tmp, seriesName = "(0,1,1)x(1,0,1) residuals", nlags = nlags)
# Formal test (NOTE: `tsdiag(fit, gof.lag = nlags)` plots a diagnostic)
Box.test(fit$residuals, lag=18) # p-value = 0.5171
Box.test(fit$residuals, lag=9)  # p-value = 0.2542
# Residuals are white noise
# Hence, (0,1,1)x(1,0,1)_4 is a valid model
model2 <- fit


# MODEL 3 ####
# Inspecting PACF of model (0,1,0)x(0,1,0)_4, seasonal lags 4 and 8 are out of limits, hence check model (0,1,0)x(2,0,0)_4
fit <- arima(log.y, order = c(0,1,0), seasonal = list(order = c(2,0,0), period = s))
check_significance(fit)
# Model is significant, check residuals
ts_acf_pacf(fit$residuals, tmp, seriesName = "(0,1,0)x(2,0,0) residuals", nlags = nlags)
# They are NOT white noise, but ACF lags 1 (regular) and 5 (polynomial) are out of limits
fit <- arima(log.y, order = c(0,1,1), seasonal = list(order = c(2,0,0), period = s))
check_significance(fit)
# Model is significant, check residuals
ts_acf_pacf(fit$residuals, tmp, seriesName = "(0,1,1)x(2,0,0) residuals", nlags = nlags)
# Not white noise, lag 5 still not captured
fit <- arima(log.y, order = c(0,1,5), seasonal = list(order = c(2,0,0), period = s))
check_significance(fit)
# Coefficients ma2, ma4 are not significant, let's force them to be zero
fit <- arima(log.y, order = c(0,1,5), fixed = c(NA, 0, NA, 0, NA, NA, NA), seasonal = list(order = c(2,0,0), period = s))
# Check the residuals
ts_acf_pacf(fit$residuals, tmp, seriesName = "(0,1,5)x(2,0,0) residuals", nlags = nlags) # They seem WN
# Formal test (NOTE: `tsdiag(fit, gof.lag = nlags)` plots a diagnostic)
Box.test(fit$residuals, lag=10) # p-value = 0.7433
# Hence, (0,1,5)x(2,0,0)_4 is a valid model
model3 <- fit


# MODEL 1 PREDICTIONS ####
# Point predictions and standard errors
log.y.pred <- predict(model1, n.ahead=12)

# Plotting real data with point predictions 
plot_forecast(y, exp(log.y.pred$pred), main = "Predictions (0,1,1)x(0,1,1)_4",
              ylab = "earning/share Coca-Cola")


# MODEL 2 PREDICTIONS ####
# Point predictions and standard errors
log.y.pred <- predict(model2, n.ahead=12)

# Plotting real data with point predictions
plot_forecast(y, exp(log.y.pred$pred), main = "Predictions (1,1,1)x(1,0,1)_4",
              ylab = "earning/share Coca-Cola")


# MODEL 3 PREDICTIONS ####
# Point predictions and standard errors
log.y.pred <- predict(model3, n.ahead=12)

# Plotting real data with point predictions
plot_forecast(y, exp(log.y.pred$pred), main = "Predictions (0,1,5)x(2,0,0)_4",
              ylab = "earning/share Coca-Cola")


# 2. For the models identified in the previous step, leave for example the last 24 real values to compare ####
#    all the models in terms of forecasting (out of sample forecasting exercise).
#    What is the best model and why is this your choice? ####

# MODEL 1 ####
# (0,1,1)x(0,1,1)_4
model1_error <- out_sample_backtesting(y, n.est = 83, n.hor = 1, log.transf = T,
                                       model1, scheme = "rolling")
model1_error$MAPE


# MODEL 2 ####
# (1,1,1)x(1,0,1)_4
model2_error <- out_sample_backtesting(y, n.est = 83, n.hor = 1, log.transf = T,
                                       model2, scheme = "rolling")
model2_error$MAPE


# MODEL 3 ####
# (0,1,5)x(2,0,0)_4
model3_error <- out_sample_backtesting(y, n.est = 83, n.hor = 1, log.transf = T,
                                       model3, fit.fixed = c(NA, 0, NA, 0, NA, NA, NA), scheme = "rolling")
model3_error$MAPE


