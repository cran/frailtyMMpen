#' Fitting penalized frailty models with clustered, multi-event and recurrent data using MM algorithm
#' 
#' @description {This formula is used to fit the penalized regression. 3 types of the models can be fitted similar to the function
#' \code{frailtyMM}. In addition, variable selection can be done by three types of penalty, LASSO, MCP and SCAD with the following
#' objective function where \eqn{\lambda} is the tuning parameter and \eqn{q} is the dimension of \eqn{\boldsymbol{\beta}},
#'  \deqn{l(\boldsymbol{\beta},\Lambda_0|Y_{obs}) + n\sum_{p=1}^{q} p(|\beta_p|, \lambda).}
#'  The BIC is computed using the following equation,
#'  \deqn{-2l(\hat{\boldsymbol{\beta}}, \hat{\Lambda}_0) + G_n(\hat{S}+1)\log(n),}
#'  where \eqn{G_n=\max{1, \log(\log(q+1))}} and \eqn{\hat{S}} is the degree of freedom.
#' }
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
#' @param penalty The penalty used for regularization, the default is "LASSO", other choices are "MCP" and "SCAD".
#' @param gam The tuning parameter for MCP and SCAD which controls the concavity of the penalty. For MCP, 
#' \deqn{p^{\prime}(\beta, \lambda)=sign(\beta)(\lambda - \frac{|\beta|}{\gamma})} and for "SCAD",
#' \deqn{p^{\prime}(\beta, \lambda)=\lambda\{I(|\beta| \leq \lambda)+\frac{(\gamma \lambda-|\beta|)_{+}}{(\gamma-1) \lambda} I(|\beta|>\lambda)\}.}
#' The default value of \eqn{\gamma} for MCP is 3 and SCAD is 3.7.
#' @param tune The sequence of tuning parameters provided by user. If not provided, the default grid will be applied.
#' @param tol The tolerance level for convergence.
#' @param maxit Maximum iterations for MM algorithm.
#' @param ... additional arguments pass to the function.
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib frailtyMMpen, .registration = TRUE
#' 
#' @details Without a given \code{tune}, the default sequence of tuning parameters are used to provide the regularization path.
#' The formula is same as the input for function \code{frailtyMM}.
#' 
#' @return An object of class \code{fmm} that contains the following fields:
#' \item{coef}{matrix of coefficient estimated from a specific model.}
#' \item{est.tht}{vector of frailty parameter estimated from a specific model.}
#' \item{lambda}{list of frailty for each observation estimated from a specific model.}
#' \item{likelihood}{vector of the observed log-likelihood given estimated parameters.}
#' \item{BIC}{vector of the BIC given estimated parameters.}
#' \item{tune}{vector of tuning parameters used for penalized regression.}
#' \item{tune.min}{tuning parameter where minimal of BIC is obtained.}
#' \item{convergence}{convergence threshold.}
#' \item{input}{The input data re-ordered by cluster id. \code{y} is the event time, \code{X} is covariate matrix and \code{d} is the status while 0 indicates censoring.}
#' \item{y}{input stopping time.}
#' \item{X}{input covariate matrix.}
#' \item{d}{input censoring indicator.}
#' \item{formula}{formula applied as input.}
#' \item{coefname}{name of each coefficient from input.}
#' \item{id}{id for individuals or clusters, {1,2...,a}. Note that, since the original id may not be the sequence starting from 1, this output
#' id may not be identical to the original id. Also, the order of id is corresponding to the returned \code{input}.}
#' \item{N}{total number of observations.}
#' \item{a}{total number of individuals or clusters.}
#' \item{datatype}{model used for fitting.}
#' 
#' @seealso \code{\link{frailtyMM}}
#' 
#' @references 
#' \itemize{
#' \item Huang, X., Xu, J. and Zhou, Y. (2022). Profile and Non-Profile MM Modeling of Cluster Failure Time and Analysis of ADNI Data. \emph{Mathematics}, 10(4), 538.
#' \item Huang, X., Xu, J. and Zhou, Y. (2023). Efficient algorithms for survival data with multiple outcomes using the frailty model. \emph{Statistical Methods in Medical Research}, 32(1), 118-132.
#' }
#' 
#' @examples 
#' 
#' data(simdataCL)
#' 
#' # Penalized regression under clustered frailty model
#' 
#' # Clustered Gamma Frailty Model
#' 
#' # Using default tuning parameter sequence
#' gam_cl1 = frailtyMMpen(Surv(time, status) ~ . + cluster(id),
#'                        simdataCL, frailty = "Gamma")
#' 
#' \donttest{
#' # Using given tuning parameter sequence
#' gam_cl2 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), 
#'                        simdataCL, frailty = "Gamma", tune = 0.1)
#' 
#' # Obtain the coefficient where minimum BIC is obtained
#' coef(gam_cl1)
#' 
#' # Obtain the coefficient with tune = 0.2.
#' coef(gam_cl1, tune = 0.2)
#' 
#' # Plot the regularization path
#' plot(gam_cl1)
#' 
#' # Get the degree of freedom and BIC for the sequence of tuning parameters provided
#' print(gam_cl1)
#' 
#' }
#' 
frailtyMMpen <- function(formula, data, frailty = "LogN", power = NULL, penalty = "LASSO", gam = NULL, tune = NULL, tol = 1e-5, maxit = 200, ...) {
  
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
      p = ncol(mx1)
      newid = seq(0, N-1, 1)
      
      if (N <= 2) {
        stop("Please check the sample size of data")
      }
      
      y = m[[1]][, 1]
      X = mx1
      d = m[[1]][, 2]
      a = N
      
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
      p = ncol(mx1)
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
    nord = order(mxid)
    mxid = mxid[nord]
    p = ncol(mx1)
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
    
  }
  
  threshold = tol
  
  if (is.null(tune)) {
    tuneseq = exp(seq(-5.5, 1, 0.25))
  } else {
    tuneseq = tune
  }
  
  if (type == "Cluster") {
    
    initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = threshold, type = 1)
    
    ini = frailtyMMcal(y, X, d, N, a, newid,
                       coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                       frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 1)
    
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMMcal(y, X, d, N, a, newid,
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0,
                         frailty = frailty, power = power, penalty = penalty, gam.val = gam, tune = tuneseq[z], maxit = maxit, threshold = threshold, type = 1)
      
      coef0 = cur$coef
      est.tht0 = cur$est.tht
      lambda0 = cur$lambda
      likelihood0 = cur$likelihood
      
      coef_all[[z]] = coef0
      est.tht_all[[z]] = est.tht0
      lambda_all[[z]] = lambda0
      likelihood_all[[z]] = likelihood0
      BIC_all[[z]] = -2*likelihood0 + max(1, log(log(p + 1)))*(sum(abs(coef0) > threshold) + 1)*log(N)
      
      if (sum(abs(coef0)) < threshold) {
        # cat(sum(abs(coef0)), "????\n")
        break
      }
      
      # cat(z, "---------\n")
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(coef_all)
    likelihood_all = unlist(likelihood_all)
    BIC_all = unlist(BIC_all)
    
    output = list(coef = coef_all,
                  est.tht = est.tht_all,
                  lambda = lambda_all,
                  likelihood = likelihood_all,
                  BIC = BIC_all,
                  tune = tuneseq[seq_len(z)],
                  tune.min = tuneseq[which.min(BIC_all)],
                  Ar = ini$Ar,
                  input = initGam$input,
                  y = y,
                  X = X,
                  d = d,
                  formula = formula,
                  coefname = coef_name,
                  id = newid + 1,
                  N = N,
                  a = a,
                  datatype = "Cluster")
  }
  
  if (type == "Multiple") {
    
    initGam = frailtyMMcal(y, X, d, N, b, NULL, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = tol, type = 2)
    
    ini =  frailtyMMcal(y, X, d, N, b, NULL,
                        coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                        frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 2)
    
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMMcal(y, X, d, N, b, NULL,
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0,
                         frailty = frailty, power = power, penalty = penalty, gam.val = gam, tune = tuneseq[z], maxit = maxit, threshold = tol, type = 2)
      
      coef0 = cur$coef
      est.tht0 = cur$est.tht
      lambda0 = cur$lambda
      likelihood0 = cur$likelihood
      
      coef_all[[z]] = coef0
      est.tht_all[[z]] = est.tht0
      lambda_all[[z]] = lambda0
      likelihood_all[[z]] = likelihood0
      BIC_all[[z]] = -2*likelihood0 + max(1, log(log(p + 1)))*(sum(abs(coef0) > 1e-6) + 1)*log(b)
      
      if (sum(abs(coef0)) < 1e-6) {
        # cat(sum(abs(coef0)), "????\n")
        break
      }
      
      # cat(z, "---------\n")
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(coef_all)
    likelihood_all = unlist(likelihood_all)
    BIC_all = unlist(BIC_all)
    
    output = list(coef = coef_all,
                  est.tht = est.tht_all,
                  lambda = lambda_all,
                  likelihood = likelihood_all,
                  BIC = BIC_all,
                  tune = tuneseq[seq_len(z)],
                  tune.min = tuneseq[which.min(BIC_all)],
                  Ar = ini$Ar,
                  input = initGam$input,
                  y = y,
                  X = X,
                  d = d,
                  formula = formula,
                  coefname = coef_name,
                  id = NULL,
                  N = N,
                  a = b,
                  datatype = "Multi-event")
  } 
  
  if (type == "Recurrent") {
    
    
    initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = maxit, threshold = threshold, type = 3)
    
    ini = frailtyMMcal(y, X, d, N, a, newid,
                       coef.ini = initGam$coef, est.tht.ini = initGam$est.tht, lambda.ini = initGam$lambda,
                       frailty = frailty, power = power, penalty = NULL, maxit = maxit, threshold = tol, type = 3)
    
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMMcal(y, X, d, N, a, newid,
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0,
                         frailty = frailty, power = power, penalty = penalty, gam.val = gam, tune = tuneseq[z], maxit = maxit, threshold = threshold, type = 3)
      
      coef0 = cur$coef
      est.tht0 = cur$est.tht
      lambda0 = cur$lambda
      likelihood0 = cur$likelihood
      
      coef_all[[z]] = coef0
      est.tht_all[[z]] = est.tht0
      lambda_all[[z]] = lambda0
      likelihood_all[[z]] = likelihood0
      BIC_all[[z]] = -2*likelihood0 + max(1, log(log(p + 1)))*(sum(abs(coef0) > 1e-6) + 1)*log(a)
      
      if (sum(abs(coef0)) < 1e-6) {
        # cat(sum(abs(coef0)), "????\n")
        break
      }
      
      # cat(z, "---------\n")
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(coef_all)
    likelihood_all = unlist(likelihood_all)
    BIC_all = unlist(BIC_all)
    
    output = list(coef = coef_all,
                  est.tht = est.tht_all,
                  lambda = lambda_all,
                  likelihood = likelihood_all,
                  BIC = BIC_all,
                  tune = tuneseq[seq_len(z)],
                  tune.min = tuneseq[which.min(BIC_all)],
                  Ar = ini$Ar,
                  input = initGam$input,
                  y = y,
                  X = X,
                  d = d,
                  formula = formula,
                  coefname = coef_name,
                  id = newid + 1,
                  N = N,
                  a = a,
                  datatype = "Recurrent")
  } 
  
 
  
  attr(output, "call") <-  Call
  class(output) = "fpen"
  output
}