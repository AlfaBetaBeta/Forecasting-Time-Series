# Forecasting time series - Box-Jenkins methodology

**Tl;dr**: SARIMA forecasting exercise based on the dataset `coca_cola_earnings.csv`. For the main script `timeSeries.R` to execute properly, it must be saved in a directory alongside the source data file `coca_cola_earnings.csv` and the additional R script `FTS_functions.R` with all the sourced-in user defined functions. This supporting R script may be opened for function inspection but should not be edited.

## Summary

The present tutorial focuses on the analysis of the time series comprising the evolution of the quarterly earnings per share of the Coca-Cola Company, in the timeframe spanning between 1983 and 2009. The aim of the tutorial is the selection of a set of linear models that may capture appropriately the structure of the time series, following the Box-Jenkins methodology. Ultimately, once a suitable set of models is proposed, and after comparing their out-of-sample forecasting performance, a final model recommendation is made.

For the purpose of model selection, the entire sample **y<sub>t</sub>** is used, comprising 107 values. For reference, the evolution of **y<sub>t</sub>** is shown below:

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/yt.png" width=100% height=100%>

By inspection, a salient feature of **y<sub>t</sub>** is the lack of stationarity. Considering that its values are positive throughout, the trend of increasing variance and the overall shape resembling an exponential, it is convenient to transform the original data via logging. By doing so, the evolution of log(**y<sub>t</sub>**) changes as follows:

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/log_yt.png" width=100% height=100%>

Once the variance is stabilised, and setting the seasonality parameter to `s <- 4` given the quarterly character of the data, inspection of the ACF below suggests not only seasonal but also cyclic effects, though these are initially ignored in the (augmented) Dickey-Fuller test.

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/log_yt_ACF_PACF.png" width=100% height=100%>

Following the Dickey-Fuller and the Osborn-Chui-Smith-Birchenhall tests for regular and seasonal differencing, respectively, the difference parameters are set to `d <- 1` and `D <- 1`. Fitting a seasonal ARIMA (0,1,0) x (0,1,0)<sub>4</sub> on the logged data leads to mean stationarity and provides the baseline ACF and PACF shown below on which to ground the decisions regarding all subsequent models.

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_010_010_4.png" width=100% height=100%>

Inspection of the ACF above reveals regular lags 1 and 3 being out of limits, alongside seasonal lag 4 and *polynomial* lag 5. Accommodating the seasonal effect at first, fitting a seasonal ARIMA (0,1,0) x (0,1,1)<sub>4</sub> leads to a significant `sma1` coefficient but also to residuals that are not white noise. Hence, the regular moving average part is incorporated in the model as a seasonal ARIMA (0,1,1) x (0,1,1)<sub>4</sub>, whereby the resulting coefficients are significant and the model residuals are indeed white noise. This can be assessed qualitatively in the figure below, either from the residual ACF/PACF or more formally from the p-values arising from the Box test for a suitable range of lags. Significance is assessed with the usual Z-value for 95% confidence after graphically checking that it constitutes a reasonable proxy for normality.

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_011_011_4.png" width=100% height=100%>
<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_011_011_4_Ljung_Box.png" width=100% height=100%>

If instead of taking differences, a regular and seasonal autoregressive component of order 1 is introduced in the model, the resulting fit leads to a seasonal coefficient close to 1 but not so for the regular coefficient, which requires setting `d <- 1` whilst leaving `D <- 0`. The resulting model has a non- significant `ar1` coefficient, effectively becoming a seasonal ARIMA (0,1,1) x (1,0,1)<sub>4</sub> with white noise residuals, displaying a very similar set of residual ACF/PACF and p-value distribution as in the previous case (qualitatively similar to the figure above).

Following a different approach, if the baseline PACF (from SARIMA (0,1,0) x (0,1,0)<sub>4</sub>) is revisited, it is noteworthy that seasonal lags 4 and 8 are out of limits. If these are accounted for with a seasonal moving average component, and considering that solely taking `d <- 1` whilst keeping `D <- 0` may be sufficient to attain mean stationarity, a seasonal ARIMA (0,1,0) x (2,0,0)<sub>4</sub> is fitted. This model has significant coefficients but its residuals are not white noise, as can be checked by inspection of the residual ACF below.

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_010_200_4.png" width=100% height=100%>

Since regular lag 1 and *polynomial* lag 5 are out of limits, it is considered convenient to incorporate a regular autoregressive component into the model, initially of order 1. Although this leads to a seasonal ARIMA (0,1,1) x (2,0,0)<sub>4</sub> with significant coefficients, the residual ACF does not fully represent white noise, as can be seen below.

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_011_200_4.png" width=100% height=100%>

Since lag 5 above does not seem to be captured by the regular and seasonal interaction, it is decided to increase the order of the regular moving average component to 5. The resulting seasonal ARIMA (0,1,5) x (2,0,0)<sub>4</sub> has significant coefficients except for `ma2` and `ma4`, which are fixed to zero. The fitted model residuals are white noise, as can be assessed by inspection of the residual ACF/PACF or the p-value distribution arising from the Box test below.

<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_015_200_4.png" width=100% height=100%>
<img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_015_200_4_Ljung_Box.png" width=100% height=100%>

Hence, the models shortlisted for out-of-sample forecasting are:

- SARIMA (0,1,1) x (0,1,1)<sub>4</sub>
- SARIMA (0,1,1) x (1,0,1)<sub>4</sub>
- SARIMA (0,1,5) x (2,0,0)<sub>4</sub>

In order to proceed with the **out-of-sample forecasting**, the last 24 values of the time series are reserved to compare with the predictions (1 horizon). A rolling scheme (of 84 shifting values) is applied here, as it is considered that it may accommodate better than the recursive scheme the instability occurring between the last quarter of 1999 and the first quarter of 2000. The logging and the inverse transformation of data are handled appropriately, and the selected metric is the mean average percentage error (MAPE), which attains the following values for the selected models:

- SARIMA (0,1,1) x (0,1,1)<sub>4</sub> &#8594; MAPE = 5.23%
- SARIMA (0,1,1) x (1,0,1)<sub>4</sub> &#8594; MAPE = 5.29%
- SARIMA (0,1,5) x (2,0,0)<sub>4</sub> &#8594; MAPE = 5.57%

Based on this metric, the recommended model is the seasonal ARIMA (0,1,1) x (0,1,1)<sub>4</sub>, although all the linear models shortlisted here showcase a similar forecasting performance. Additionally, the recommended model is simple and easier to interpret, requiring only the determination of two coefficients. The latter factor may be particularly relevant when historical data is scarce.

There is one last consideration, however, initially outside of the scope of this tutorial. When comparing the point predictions of the shortlisted models for 12 horizons beyond the last observed value of the original series, there is a strong contrast between the predicted trends, as shown below.

<p align="middle">
  <img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_011_011_4_predictions.png" width=30% height=30%>
  <img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_111_101_4_predictions.png" width=30% height=30%>
  <img src="https://github.com/AlfaBetaBeta/Forecasting-Time-Series/blob/master/img/model_015_200_4_predictions.png" width=30% height=30%>
</p>

Considering how similar the MAPE values are for all models, such outcome was not anticipated and suggests that SARIMA (0,1,5) x (2,0,0)<sub>4</sub> may be incorporating undesirable memory effects from the instability between years 1999 and 2000. This could be caused by the high order of the regular moving average component, beyond the value of the seasonality parameter, although this remains the subject of further analysis.
