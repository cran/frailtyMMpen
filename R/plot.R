#' Plot the regularization path
#' 
#' @description Plot the whole regularization path run by frailtyMMpen
#' @usage
#' ##S3 method for class "fpen"
#' @param x Object with class "fpen"
#' @param ... Further arguments pass to or from other methods
#' @method plot fpen
#' @export
#' 
plot.fpen <- function(x, ...) {
  xaxis = log(x$tune)
  yaxis = t(x$coef)

  n = ncol(yaxis)
  matplot(xaxis, yaxis, lty = 1, xlab = "log(tune)", ylab = "coefficients", type = "l")
  abline(v = log(x$tune.min), col = "red", lwd = 3, lty = 2)
}

#' Plot the baseline hazard or the predicted hazard based on the new data
#' 
#' @description Both the cumulative hazard and the survival curves can be plotted.
#' @usage
#' ##S3 method for class "fmm"
#' @param object Object with class "fmm"
#' @param newdata The new data for prediction of hazard
#' @param surv Plot survival curve instead of cumulative hazard, the default is \code{FALSE}
#' @param ... Further arguments pass to or from other methods
#' @method plot fmm
#' @export
#' 
#' @details If parameter \code{newdata} is given, the plot is based on the predicted hazard while if it is not given,
#' the plot is based on the baseline hazard. To construct the new data, please refer to the detailed description from 
#' function \code{predict.fmm} and the following example.
#' 
#' @seealso \code{\link{predict.fmm}}
#' 
#' @examples 
#' 
#' \donttest{
#' gam_re = frailtyMM(Surv(tstart, tstop, status) ~  sex + treat + cluster(id), cgd, frailty = "gamma")
#' 
#' # Plot the survival curve based on baseline hazard
#' plot(gam_re, surv = TRUE)
#' 
#' # Construct new data and plot the cumulative hazard based on new data
#' newre = c(1, 1, 2)
#' names(newre) = c(gam_re$coefname, "id")
#' plot(gam_re, newdata = newre)
#' }
plot.fmm <- function(x, newdata = NULL, surv = FALSE, ...) {
  
  df = predict(object = x, newdata = newdata, surv = surv)
  
  if (is.null(newdata)) {
    if (surv == TRUE) {
      yax = "Baseline Survival curve"
    } else {
      yax = "Baseline Cumulative hazard"
    }
  } else {
    if (surv == TRUE) {
      yax = "Survival curve"
    } else {
      yax = "Cumulative hazard"
    }
  }
  
  plot(df$time, df$estimate, type = "l", col = "red",
       xlab = "time", ylab = yax, ylim = c(0, max(df$estimateub) + 0.1), lwd = 2)
  lines(df$time, df$estimateub, col = "red", lty = 3, lwd = 1.5)
  lines(df$time, df$estimatelb, col = "red", lty = 3, lwd = 1.5)
  
}