# Supporting functions for all Forecasting Time Series scripts

ts_acf_pacf <- function(y, x.time = 1:length(y), seriesName = "s", nlags = length(y)-1, marg = 5) {
    par(mar = c(marg, marg, marg, marg)) 
    par(mfrow = c(3,1))
    plot(x.time, y, xlab = "Time", ylab = seriesName)
    lines(x.time, y)
    plot(acf(y, plot = F, lag.max = nlags)[1:nlags], main = paste0("Series ", seriesName))
    plot(pacf(y, plot = F, lag.max = nlags), main = paste0("Series ", seriesName))
}

norm_graph_check <- function(y, y_lim) {
    mu <- mean(y)
    sigma <- sd(y)
    
    xx <- seq(mu - 3 * sigma, mu + 3 * sigma, length = 100)
    yy <- dnorm(xx, mu, sigma)
    
    par(mar = c(3,3,3,3)) 
    par(mfrow = c(1,1))
    hist(y, prob = T, ylim = c(0,y_lim), xlim = c(mu - 3 * sigma, mu + 3 * sigma), col = "grey")
    lines(density(y), lwd = 2)
    lines(xx, yy, lwd = 2, col = "blue")
}

plot_forecast <- function (y, series.pred, main = "Predictions", ylab = "feature", legend.loc = "topleft", marg = 4) {
    new <- c(y, series.pred)
    par(mar = c(marg, marg, marg, marg))
    par(mfrow = c(1,1))
    plot.ts(new, main = main, ylab = ylab, col = 3, lwd = 2)
    lines(y, col = 4, lwd = 2)
    legend(legend.loc, legend = c("Predictions","Historical"),
           col = c(3,4), bty = "n", lwd = 2)
}

# plot_forecast_data <- function (prediction.y, prediction.x, real.y, real.x, leg.loc = "topleft") {
#     par(mfrow = c(1,1))
#     plot(real.x, real.y, type = "b", col = "red")
#     lines(real.x, real.y, type = "b", col = "blue")
#     points(prediction.x, prediction.y)
#     legend(leg.loc, c("real","forecast"),
#            col = c("red","blue"), pch = c(1,1), bty = "n" )
# }

check_significance <- function (fit) {
    sig_mat <- matrix(0, nrow = 3, ncol = length(names(fit$coef)),
                      dimnames = list(c("CI upper bound", "CI lower bound", "Significant"),
                                      names(fit$coef)))
    for (cn in colnames(fit$var.coef)) {
        sig_mat["CI upper bound", cn] <- fit$coef[cn] + 1.96 * sqrt(fit$var.coef[cn, cn])
        sig_mat["CI lower bound", cn] <- fit$coef[cn] - 1.96 * sqrt(fit$var.coef[cn, cn])
        sig_mat["Significant", cn] <- (sig_mat["CI lower bound", cn] > 0) || (sig_mat["CI upper bound", cn] < 0)
    }
    return(sig_mat)
}

out_sample_backtesting <- function (y, n.est, n.hor, log.transf = FALSE,
                                    model, fit.fixed = NULL,
                                    scheme = "recursive") {
    n <- length(y)
    n.estimation <- n.est               
    n.forecasting <- n - n.estimation
    horizons <- n.hor
    
    predicc <- matrix(0, nrow = n.forecasting, ncol = horizons)
    real <- matrix(0, nrow = n.forecasting, ncol = 1)
    real[] <- y[(n.estimation + 1):n]
    MSFE <- matrix(0, nrow = horizons, ncol = 1)
    MAPE <- matrix(0, nrow = horizons, ncol = 1)
    
    for (Periods_ahead in 1:horizons) {
        
        for (i in 1:n.forecasting) {
            
            if (scheme == "recursive") {
                aux.y <- y[1:(n.estimation - Periods_ahead + i)]
            } else {
                aux.y <- y[i:(n.estimation - Periods_ahead + i)]
            }
            
            if (log.transf) {
                fit <- arima(log(aux.y), order = model$arma[c(1,6,2)], fixed = fit.fixed,
                             seasonal = list(order = model$arma[c(3,7,4)], period = model$arma[5]))
            } else {
                fit <- arima(aux.y, order = model$arma[c(1,6,2)], fixed = fit.fixed,
                             seasonal = list(order = model$arma[c(3,7,4)], period = model$arma[5]))
            }
            
            y.pred <- predict(fit, n.ahead = Periods_ahead)
            if (log.transf) {
                predicc[i, Periods_ahead] <- exp(y.pred$pred[Periods_ahead])
            } else {
                predicc[i, Periods_ahead] <- y.pred$pred[Periods_ahead]
            }
        }
        error <- real - predicc[,Periods_ahead]
        MSFE[Periods_ahead,] <- mean(error^2)
        MAPE[Periods_ahead,] <- mean(abs(error/real)) * 100
    }
    return(list("MSFE" = MSFE, "MAPE" = MAPE))
}