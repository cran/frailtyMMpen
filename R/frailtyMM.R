#' Fitting frailty models with clustered, multi-event and recurrent data using MM algorithm
#' 
#' @description {
#' This formula is used to fit the non-penalized regression. 3 types of the models can be fitted, shared frailty model where
#' hazard rate of \eqn{j^{th}} object in \eqn{i^{th}} cluster is
#' \deqn{\lambda_{ij}(t|\omega_i) = \lambda_0(t) \omega_i \exp(\boldsymbol{\beta}' \mathbf{X_{ij}}).}
#' The multi-event frailty model with different baseline hazard of different event and the hazard rate of \eqn{j^{th}} event for individual \eqn{i^{th}} is 
#' \deqn{\lambda_{ij}(t|\omega_i) = \lambda_{0j}(t) \omega_i \exp(\boldsymbol{\beta}' \mathbf{X_{ij}}).}
#' }
#' The recurrent event model where the \eqn{j^{th}} event of individual \eqn{i} has observed feature \eqn{\mathbf{X_{ij}}},
#' \deqn{\lambda_{ij}(t|\omega_i) = \lambda_0(t) \omega_i \exp(\boldsymbol{\beta}' \mathbf{X_{ij}}).}
#' 
#' @param formula Formula where the left hand side is an object of the type \code{Surv}
#' and the right hand side contains the variables and additional specifications. 
#' \code{+cluster()} function specify the group id for clustered data or individual id for recurrent data.
#' \code{+event()} function specify the event id for multi-event data (only two events are allowed).
#' @param data The \code{data.frame} where the formula argument can be evaluated.
#' @param frailty The frailty used for model fitting. The default is "LogN", other choices are
#' "InvGauss", "Gamma" and "PVF". (Note that the computation time for PVF family will be slow 
#' due to the non-explicit expression of likelihood function)
#' @param power The power used if PVF frailty is applied.
#' @param tol The tolerance level for convergence.
#' @param maxit Maximum iterations for MM algorithm.
#' @param ... additional arguments pass to the function.
#' @import mgcv survival
#' @importFrom graphics abline lines matplot
#' @importFrom stats model.frame model.matrix pnorm predict printCoefmat
#' @export
#' 
#' @details To run the shared frailty model, \code{Surv(tstop, status)} formula should be applied along with \code{+cluster()} to specify the
#' corresponding clusters, if you want to run the simple frailty model without shared frailty, you do not need to use \code{+cluster()} and the
#' formula only contains the name of the covariates. To run the multi-event model, 
#' \code{Surv(tstop, status)} formula should be applied along with \code{+event()} to specify the corresponding events. If multi-event data
#' is fitted, please use {1,2...,K} to denote the event id from the input data. To run the recurrent event model, 
#' \code{Surv(tstart, tstop, status)} formula should be applied along with \code{+cluster()} where the cluster here denotes the individual id and
#' each individual may have many observed events at different time points.
#' 
#' The default frailty will be log-normal frailty, in order to fit other frailty models, simply set parameter \code{frailty} as "InvGauss", "Gamma" or "PVF",
#' the parameter \code{power} is only used when \code{frailty}=PVF and since the likelihood of PVF (tweedie) distribution is approximated using 
#' \code{Tweedie} function from package mgcv, 1<\code{power}<2.
#' 
#' @return An object of class \code{fmm} that contains the following fields:
#' \item{coef}{coefficient estimated from a specific model.}
#' \item{est.tht}{frailty parameter estimated from a specific model.}
#' \item{lambda}{frailty for each observation estimated from a specific model.}
#' \item{likelihood}{The observed log-likelihood given estimated parameters.}
#' \item{input}{The input data re-ordered by cluster id. \code{y} is the event time, \code{X} is covariate matrix and \code{d} is the status while 0 indicates censoring.}
#' \item{frailty}{frailty used for model fitting.}
#' \item{power}{power used for model fitting is PVF frailty is applied.}
#' \item{iter}{total number of iterations.}
#' \item{convergence}{convergence threshold.}
#' \item{formula}{formula applied as input.}
#' \item{coefname}{name of each coefficient from input.}
#' \item{id}{id for individuals or clusters, {1,2...,a}. Note that, since the original id may not be the sequence starting from 1, this output
#' id may not be identical to the original id. Also, the order of id is corresponding to the returned \code{input}.}
#' \item{N}{total number of observations.}
#' \item{a}{total number of individuals or clusters.}
#' \item{datatype}{model used for fitting.}
#' 
#' @references 
#' \itemize{
#' \item Huang, X., Xu, J. and Zhou, Y. (2022). Profile and Non-Profile MM Modeling of Cluster Failure Time and Analysis of ADNI Data. \emph{Mathematics}, 10(4), 538.
#' \item Huang, X., Xu, J. and Zhou, Y. (2023). Efficient algorithms for survival data with multiple outcomes using the frailty model. \emph{Statistical Methods in Medical Research}, 32(1), 118-132.
#' }
#' 
#' @examples 
#' 
#' # Kidney data fitted by Clustered Inverse Gaussian Frailty Model
#' 
#' \donttest{
#' InvG_real_cl = frailtyMM(Surv(time, status) ~ age + sex + cluster(id),
#'                          kidney, frailty = "InvGauss")
#' InvG_real_cl
#' 
#' # Cgd data fitted by Recurrent Log-Normal Frailty Model
#' 
#' logN_real_re = frailtyMM(Surv(tstart, tstop, status) ~ sex + treat + cluster(id),
#'                          cgd, frailty = "Gamma")
#' logN_real_re
#' }
#' 
#' # Simulated data example
#' 
#' data(simdataCL)
#'
#' # Parameter estimation under different model structure and frailties
#' 
#' # Clustered Gamma Frailty Model
#' gam_cl = frailtyMM(Surv(time, status) ~ . + cluster(id), 
#'                    simdataCL, frailty = "Gamma")
#' 
#' \donttest{
#' # Clustered Log-Normal Frailty Model
#' logn_cl = frailtyMM(Surv(time, status) ~ . + cluster(id), 
#'                     simdataCL, frailty = "LogN")
#' 
#' # Clustered Inverse Gaussian Frailty Model
#' invg_cl = frailtyMM(Surv(time, status) ~ . + cluster(id), 
#'                     simdataCL, frailty = "InvGauss")
#'                    
#' data(simdataME)
#' 
#' # Multi-event Gamma Frailty Model
#' gam_me = frailtyMM(Surv(time, status) ~ . + cluster(id), 
#'                    simdataCL, frailty = "Gamma")
#' 
#' 
#' # Multi-event Log-Normal Frailty Model
#' logn_me = frailtyMM(Surv(time, status) ~ . + event(id), 
#'                     simdataME, frailty = "LogN")
#' 
#' # Multi-event Inverse Gaussian Frailty Model
#' invg_me = frailtyMM(Surv(time, status) ~ . + event(id),
#'                     simdataME, frailty = "InvGauss")
#' 
#' data(simdataRE)
#' 
#' # Recurrent event Gamma Frailty Model
#' gam_re = frailtyMM(Surv(start, end, status) ~ . + cluster(id),
#'                    simdataRE, frailty = "Gamma")
#' 
#' # Recurrent event Log-Normal Frailty Model
#' logn_re = frailtyMM(Surv(start, end, status) ~ . + cluster(id),
#'                    simdataRE, frailty = "LogN")
#' 
#' # Recurrent event Inverse Gaussian Frailty Model
#' invg_re = frailtyMM(Surv(start, end, status) ~ . + cluster(id), 
#'                     simdataRE, frailty = "InvGauss")
#' }
#' 
#' # Obtain the summary statistics under fitted model
#' 
#' coef(gam_cl)
#' summary(gam_cl)
#' 
frailtyMM <- function(formula, data, frailty = "LogN", power = NULL, tol = 1e-5, maxit = 200, ...) {
  
  Call <- match.call()
  
  if(!inherits(formula, "formula")) {
    stop("please provide formula object for formula")
  }
  
  if(!inherits(data, "data.frame")) {
    stop("please provide data.frame type for data")
  }
  
  m <- model.frame(formula, data)
  mx <- model.matrix(formula, data)
  
  if (ncol(m[[1]]) == 2) {
    
    cluster_id <- grep("cluster", names(m))
    event_id <- grep("event", names(m))
    
    if (length(cluster_id) == 0 && length(event_id) == 0) {
      
      type = "Cluster"
      mx1 = mx[, -c(1)]
      coef_name = colnames(mx1)
      
      N = nrow(mx1)
      newid = seq(0, N-1, 1)
      
      if (N <= 2) {
        stop("Please check the sample size of data")
      }
      
      y = m[[1]][, 1]
      X = mx1
      d = m[[1]][, 2]
      a = N
      
      initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 1)
      
      output = frailtyMMcal(y, X, d, N, a, newid,
                            coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                            frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 1)
      
      ret = list(coef = output$coef,
                 est.tht = output$est.tht,
                 lambda = output$lambda,
                 likelihood = output$likelihood,
                 Ar = output$Ar,
                 input = output$input,
                 frailty = frailty,
                 power = power,
                 iter = output$iter,
                 convergence = output$convergence,
                 formula = formula,
                 coefname = coef_name,
                 id = newid + 1,
                 N = N,
                 a = a,
                 datatype = "Cluster")
    }
    
    if (length(cluster_id) == 1) {
      
      type = "Cluster"
      pb = unlist(gregexpr('\\(', names(m)[cluster_id])) + 1
      pe = unlist(gregexpr('\\)', names(m)[cluster_id])) - 1
      clsname = substr(names(m)[cluster_id], pb, pe)
      remove_cluster_id <- grep(clsname, names(m))
      mx1 = mx[, -c(1, remove_cluster_id)]
      mxid = mx[, cluster_id]
      coef_name = colnames(mx1)
      
      nord = order(mxid)
      mxid = mxid[nord]
      
      N = length(mxid)
      newid = rep(0, N)
      
      if (N <= 2) {
        stop("Please check the sample size of data")
      }
      
      for (i in 2:N) {
        if (mxid[i] > mxid[i-1]) {
          newid[i:N] = newid[i:N] + 1
        }
      }
      
      y = m[[1]][nord, 1]
      X = mx1[nord, ]
      d = m[[1]][nord, 2]
      a = max(newid) + 1
      
      initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 1)
      
      output = frailtyMMcal(y, X, d, N, a, newid,
                            coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                            frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 1)
      
      ret = list(coef = output$coef,
                 est.tht = output$est.tht,
                 lambda = output$lambda,
                 likelihood = output$likelihood,
                 Ar = output$Ar,
                 input = output$input,
                 frailty = frailty,
                 power = power,
                 iter = output$iter,
                 convergence = output$convergence,
                 formula = formula,
                 coefname = coef_name,
                 id = newid + 1,
                 N = N,
                 a = a,
                 datatype = "Cluster")
      
    }
    
    if (length(event_id) == 1) {
      
      type = "Multiple"
      pb = unlist(gregexpr('\\(', names(m)[event_id])) + 1
      pe = unlist(gregexpr('\\)', names(m)[event_id])) - 1
      evsname = substr(names(m)[event_id], pb, pe)
      remove_event_id <- grep(evsname, names(m))
      mx1 = mx[, -c(1, remove_event_id)]
      mxid = mx[, event_id]
      coef_name = colnames(mx1)
      
      mxid_info = table(mxid)
      n = length(mxid_info)
      b = min(mxid_info)
      p = ncol(mx1)
      if (b != max(mxid_info)) {
        stop("every subject should have same number of events")
      }
      
      nord = order(mxid)
      N = length(nord)
      mx1 = mx1[nord, ]
      X = mx1[nord, ]
      y = m[[1]][nord, 1]
      d = m[[1]][nord, 2]
      
      initGam = frailtyMMcal(y, X, d, N, b, NULL, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 2)
      
      output = frailtyMMcal(y, X, d, N, b, NULL,
                            coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                            frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 2)
      
      ret = list(coef = output$coef,
                 est.tht = output$est.tht,
                 lambda = output$lambda,
                 likelihood = output$likelihood,
                 Ar = output$Ar,
                 input = output$input,
                 frailty = frailty,
                 power = power,
                 iter = output$iter,
                 convergence = output$convergence,
                 formula = formula,
                 coefname = coef_name,
                 id = NULL,
                 N = N,
                 a = b,
                 datatype = "Multi-event")
      
    }
  }
  
  if (ncol(m[[1]]) == 3) {
    
    type = "Recurrent"
    cluster_id <- grep("cluster", names(m))
    pb = unlist(gregexpr('\\(', names(m)[cluster_id])) + 1
    pe = unlist(gregexpr('\\)', names(m)[cluster_id])) - 1
    clsname = substr(names(m)[cluster_id], pb, pe)
    remove_cluster_id <- grep(clsname, names(m))
    mx1 = mx[, -c(1, remove_cluster_id)]
    mxid = mx[, cluster_id]
    coef_name = colnames(mx1)
    
    mxid_info = table(mxid)
    p = ncol(mx1)

    nord = order(mxid)
    mxid = mxid[nord]
    
    N = length(mxid)
    newid = rep(0, N)
    
    if (N <= 2) {
      stop("Please check the sample size of data")
    }
    
    for (i in 2:N) {
      if (mxid[i] > mxid[i-1]) {
        newid[i:N] = newid[i:N] + 1
      }
    }
    
    y = m[[1]][nord, 2]
    X = mx1[nord, ]
    d = m[[1]][nord, 3]
    a = max(newid) + 1
    
    initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 3)
    
    output = frailtyMMcal(y, X, d, N, a, newid,
                          coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                          frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 3)
    
    ret = list(coef = output$coef,
               est.tht = output$est.tht,
               lambda = output$lambda,
               likelihood = output$likelihood,
               Ar = output$Ar,
               input = output$input,
               frailty = frailty,
               power = power,
               iter = output$iter,
               convergence = output$convergence,
               formula = formula,
               coefname = coef_name,
               id = newid + 1,
               N = N,
               a = a,
               datatype = "Recurrent")
    
  }
  
  attr(ret, "call") <-  Call
  class(ret) <- "fmm"
  return(ret)
}