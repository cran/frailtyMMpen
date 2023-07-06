#' Provide the summary for the model fitting 
#' 
#' @importFrom numDeriv hessian
#' @usage
#' ##S3 method for class "fmm"
#' @param object Object with class "fmm", generated from \code{frailtyMM}.
#' @param ... Ignored
#' 
#' @details the summary for the model of frailtyMM.
#' The standard error and p-value of estimated parameters are based on Fisher Information matrix.
#'
#' @method summary fmm
#' @export

summary.fmm <- function(object, ...) {
  
  model = object
  frailtyc = switch(model$frailty, "Gamma" = 0, "Log-Normal" = 1, "Inverse Gaussian" = 2, "PVF" = 3)
  datatype = switch(model$datatype, "Cluster" = 1, "Multi-event" = 2, "Recurrent" = 3)
  p = length(model$coef)
  
  if (is.null(model$power)) {
    power = 0.0
  } else {
    power = model$power
  }
  
  Fisher = -numDeriv::hessian(logLik,
                              x = c(model$est.tht, model$coef),
                              method="Richardson",
                              data = model$input,
                              lambda = model$lambda,
                              frailtyc = frailtyc,
                              id = model$id - 1,
                              N = model$N,
                              a = model$a,
                              p = p,
                              power = power,
                              type = datatype)
  
  sd = diag(solve(Fisher))
  
  th_sd = sqrt(sd[1])
  coef_sd = sqrt(sd[-1])
  
  zth = model$est.tht / th_sd
  zcoef = model$coef / coef_sd
  
  pval_th = 2 - 2*pnorm(abs(zth))
  pval_coef = 2 - 2*pnorm(abs(zcoef))
  
  ret = list(est.tht = model$est.tht,
             coef = model$coef,
             th_sd = th_sd,
             coef_sd = coef_sd,
             zth = zth,
             zcoef = zcoef,
             pval_th = pval_th,
             pval_coef = pval_coef,
             frailty = model$frailty,
             power = model$power,
             datatype = model$datatype,
             iter = model$iter,
             coefname = model$coefname,
             convergence = model$convergence,
             likelihood = model$likelihood)
  
  class(ret) = "fmm_summary"
  attr(ret, "call") = attr(model, "call")
  return(ret)
}