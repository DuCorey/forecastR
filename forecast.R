#' Forecasting utility functions and models

#' Compact forecastings reduces the memory of the resulting model.
#' Simply call the compact forecasting model with the original data to use it for forecasting
#' Avoids keeping in memory the data twice

#' Result forecast returns the final numerical prediction which varies across forecasting
#' frameworks

#' Update allows you to speed add new data to a forecasting model while keeping the parameters
#' the same

#' Update forecast methods so that different forecasting frameworks work better together

#' imports
library(forecast)
library(tsintermittent)
library(MAPA)

#' functions
compact_forecast <- function(object, ...) {
    UseMethod("compact_forecast")
}


result_forecast <- function(object, ...) {
    UseMethod("result_forecast")
}


result_forecast.forecast <- function(fcast) {
    #' Return the numerical forecast result only from the forecast framework
    return(as.numeric(fcast$mean))
}


ets <- function(x, ...) {
    return(forecast::ets(as.ts(x), ...))
}


forecast.ets <- function(obj, h, ...) {
    fcast <- forecast::forecast.ets(object = obj, h = h, ...)
    fcast$h <- h
    return(fcast)
}


update.ets <- function(model, newdata, ...) {
    #' How to update an ets forcast update with new data
    #' Reuse the component type from the initial model but reestimate the coefficients
    modelcomponents <- paste(model$components[1], model$components[2],
                             model$components[3], sep = "")
    damped <- (model$components[4] == "TRUE")

    return(ets(newdata, model = modelcomponents, damped = damped, ...))
}


compact_forecast.ets <- function(x) {
    #' A compact form for the model which contains it's parameters so that one
    #' use it for other forecasts in the future. This avoids having to save
    #' the heavy time series data in the object call used to make the object.
    #' We save the parameters which are optimised for later use.

    #' Unfortunately these are the only parameters that can be kept.
    #' Any other parameters when passed back to the model will change the output.

    structure(
        list(
            components = x$components
        ),
        class = "etscompact"
    )
}


forecast.etscompact <- function(obj, x, h) {
    #' Forecast a compact ets this requires the original data to be sent again
    #' obj : compactets object
    #' x : original data
    #' h : forecast horizon
    modeltype <- paste(obj$components[1], obj$components[2],
                       obj$components[3], sep = "")
    damped <- as.logical(obj$components[4])

    model <- ets(x, model = modeltype, damped = damped)
    fcast <- forecast(model, h)

    return(fcast)
}


ses <- function(x, ...) {
    fcast <- forecast::ets(as.ts(x), "ANN", opt.crit = "mse", ...)
    fcast$method <- fcast$model$method <- "Simple exponential smoothing"
    fcast$model$call <- match.call()
    fcast$series <- deparse(substitute(x))

    structure(
        fcast,
        class = "ses"
    )
}


forecast.ses <- function(obj, h, ...) {
    out <- forecast::forecast.ets(object = obj, h = h, ...)
    out$h <- h
    return(out)
}


update.ses <- function(obj, newdata) {
    #' How we want to update a ses object with new data
    #' We can't since the ses already has a bunch of parameters already set in the ets call
    return(ses(newdata))
}


compact_forecast.ses <- function(obj) {
    structure(
        NA,
        class = "sescompact"
    )
}


forecast.sescompact <- function(obj, x, h) {
    #' Forecast a compact ses this requires the original data to be sent again
    #' obj : compactets object
    #' x : original data
    #' h : forecast horizon
    model <- ses(x)
    fcast <- forecast(model, h)

    return(fcast)
}


theta <- function(x, h = 1, ...) {
    #' The object returned by a thetaf is a forecast object. This means that
    #' there isn't a forecasting function that can be called on it. This object
    #' allows us to then call a forecast on it.
    #' Default h = 1 because we don't care about the forecast yet. This is fast.
    model <- forecast::thetaf(as.ts(x), h = h, ...)

    structure(
        model,
        class = "Theta"
    )
}


forecast.Theta <- function(obj, h, ...) {
    if (length(obj$mean) == h) {
        out <- x
    } else {
        ## Create and return the theta object with the correct forecast horizon
        return(forecast(theta(obj$x, h = h), h = h))
    }
    out$h <- h

    structure(
        out,
        class = "forecast"  # Manually set the class back to forecast
    )
}


update.Theta <- function(obj, newdata) {
    #' How we want to update a Theta object with newdata
    #' We can't since Theta needs to recalculate all the parameteres
    return(theta(newdata))
}


ARIMA <- function(x, ...) {
    #' ARIMA forecast
    return(forecast::auto.arima(as.ts(x), ...))
}


forecast.ARIMA <- function(obj, h, ...) {
    fcast <- forecast::forecast(object = obj, h = h, ...)
    fcast$h <- h
    return(fcast)
}


update.ARIMA <- function(model, newdata, ...) {
    #' Needs to reupdate everything
    #' Unfortunately there's no easy way to get the values of
    #' p,d,q from an existing model
    return(ARIMA(newdata))
}


compact_forecast.ARIMA <- function(obj) {
    structure(
        list(
            NA
        ),
        class = "ARIMAcompact"
    )
}


snaive <- function(x, ...) {
    #' snaive forecast
    fcast <- forecast::snaive(as.ts(x), ...)
    fcast$series <- as.ts(x)

    structure(
        fcast,
        class = "snaive"
    )
}


forecast.snaive <- function(obj, h, ...) {
    out <- forecast::snaive(y = obj$series, h = h)
    out$h <- h
    return(out)
}


update.snaive <- function(model, newdata, ...) {
    return(snaive(newdata))
}


compact_forecast.snaive <- function(obj) {
    structure(
        NA,
        class = "snaivecompact"
    )
}


croston <- function(data, f.type = c("SBA.base", "SBA.opt"), ...) {
    f.type <- match.arg(f.type, c("SBA.base", "SBA.opt"))
    x <- as.numeric(data)

    model <- switch(f.type,
                    SBA.base = tsintermittent::crost(x, h = 0, w = 0.05,
                                                     type = "sba",
                                                     init.opt = FALSE, ...),
                    SBA.opt = tsintermittent::crost(x, h = 0, type = "sba", ...))
    model$x <- x
    model$type <- f.type
    model$fitted <- data

    structure(
        model,
        class = "croston"
    )
}


forecast.croston <- function(obj, h, ...) {
    #' Since we already have the croston object, we don't need to redo any of
    #' the optimisations in the crosotn method
    out <- tsintermittent::crost(obj$x, h = h, w = obj$weights, init = obj$initial,
                                 type = obj$model, init.opt = FALSE, ...)
    out$h <- h

    structure(
        out,
        class = "crostonforecast"
    )
}


update.croston <- function(obj, newdata, ...) {
    #' How we want to update a croston object with new data
    #' We have to recalculate all the parameters
    return(croston(newdata, obj$type, ...))
}


compact_forecast.croston <- function(obj) {
    #' Compact representation of the croston object.
    #' This is all we need to recreate the original object given the same data
    structure(
        list(
            weights = obj$weights,
            initial = obj$initial,
            model = obj$model
        ),
        class = "crostoncompact"
    )
}


compact_forecast.crostoncompact <- function(obj) {
    return(obj)
}


forecast.crostoncompact <- function(obj, x, h) {
    obj$x <- x
    return(forecast.croston(obj, h))
}


result_forecast.crostonforecast <- function(fcast) {
    #' Return the numerical forecast result only
    return(as.numeric(fcast$frc.out))
}


MAPA <- function(data, agg, ...) {
    mapafit <- MAPA::mapaest(as.ts(data), ppy = agg, type = "ets", model = "ZZZ", ...)

    structure(
        list(
            x = as.ts(data),
            mapafit = mapafit,
            agg = agg
        ),
        class = "MAPA"
    )
}


fitted.MAPA <- function(obj) {
    return(obj$mapafit[, "fitted"])
}


forecast.MAPA <- function(obj, h, ifh = 0, ...) {
    fcast <- mapafor(obj$x, obj$mapafit, fh = h, ifh = ifh, comb = "w.mean")

    structure(
        fcast,
        class = "MAPAforecast"
    )
}


update.MAPA <- function(obj, newdata, ...) {
    return(MAPA(newdata, obj$agg))
}


compact_forecast.MAPA <- function(obj) {
    structure(
        list(
            mapafit = obj$mapafit
        ),
        class = "MAPAcompact"
    )
}


forecast.MAPAcompact <- function(obj, x, h, ...) {
    #' How we forecast a compacted MAPA object
    fcast <- mapafor(x, obj$mapafit, fh = h, ifh = 0, comb = "w.mean")

    structure(
        fcast,
        class = "MAPAforecast"
    )
}


result_forecast.MAPAforecast <- function(obj) {
    return(as.numeric(obj$outfor))
}
