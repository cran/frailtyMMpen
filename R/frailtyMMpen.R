#' Fitting penalized frailty models with clustered, multi-event and recurrent data using MM algorithm
#' 
#' @description {This formula is used to fit the penalized regression. 3 types of the models can be fitted similar to the function
#' \code{frailtyMM}. In addition, variable selection can be done by three types of penalty, LASSO, MCP and SCAD with the following
#' objective function where \eqn{\lambda} is the tuning parameter and \eqn{q} is the dimension of \eqn{\boldsymbol{\beta}},
#'  \deqn{l(\boldsymbol{\beta},\Lambda_0|Y_{obs}) - n\sum_{p=1}^{q} p(|\beta_p|, \lambda).}
#'  The BIC is computed using the following equation,
#'  \deqn{-2l(\hat{\boldsymbol{\beta}}, \hat{\Lambda}_0) + G_n(\hat{S}+1)\log(n),}
#'  where \eqn{G_n=\max\{1, \log(\log(q+1))\}} and \eqn{\hat{S}} is the degree of freedom.
#'  
#'  Surrogate function is also derived for penalty part for efficient estimation of penalized regression, similar to the notation used
#'  in \code{\link{frailtyMM}}, we let \eqn{\boldsymbol{\alpha}} be the collection of all parameters and baseline hazard function. Given that,
#'  
#'  \if{html}{\figure{fig15.png}{options: style="width:750px;max-width:75\%;"}}
#'  \if{latex}{\out{\begin{center}}{\figure{fig15.png}{options: width=5in}}\out{\end{center}}}
#'  
#'  by local quadratic approximation,
#'  
#'  \if{html}{\figure{fig16.png}{options: style="width:750px;max-width:75\%;"}}
#'  \if{latex}{\out{\begin{center}}{\figure{fig16.png}{options: width=5in}}\out{\end{center}}}
#'  
#'  And thus, the surrogate function given \eqn{k^{th}} iteration result is as follows,
#'  
#'  \if{html}{\figure{fig17.png}{options: style="width:750px;max-width:75\%;"}}
#'  \if{latex}{\out{\begin{center}}{\figure{fig17.png}{options: width=5in}}\out{\end{center}}}
#' }
#' 
#' @param formula Formula where the left hand side is an object of the type \code{Surv}
#' and the right hand side contains the variables and additional specifications. 
#' \code{+cluster()} function specify the group id for clustered data or individual id for recurrent data.
#' \code{+event()} function specify the event id for multi-event data (only two events are allowed).
#' @param data The \code{data.frame} where the formula argument can be evaluated.
#' @param frailty The frailty used for model fitting. The default is "lognormal", other choices are
#' "invgauss", "gamma" and "pvf". (Note that the computation time for PVF family will be slow 
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
#' \item{coef}{matrix of coefficient estimated from a specific model where each column correponds to an input tuning parameter.}
#' \item{est.tht}{vector of frailty parameters estimated from a specific model with respect to each tuning parameter.}
#' \item{lambda}{list of frailty for each observation estimated from a specific model with respect to each tuning parameter.}
#' \item{likelihood}{vector of the observed log-likelihood given estimated parameters with respect to each tuning parameter.}
#' \item{BIC}{vector of the BIC given estimated parameters with respect to each tuning parameter.}
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
#'                        simdataCL, frailty = "gamma")
#' 
#' \donttest{
#' # Using given tuning parameter sequence
#' gam_cl2 = frailtyMMpen(Surv(time, status) ~ . + cluster(id), 
#'                        simdataCL, frailty = "gamma", tune = 0.1)
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
frailtyMMpen <- function(formula, data, frailty = "gamma", power = NULL, penalty = "LASSO", gam = NULL, tune = NULL, tol = 1e-5, maxit = 200, ...) {
  
  Call <- match.call()
  
  if(!inherits(formula, "formula")) {
    stop("please provide formula object for formula")
  }
  
  if(!inherits(data, "data.frame")) {
    stop("please provide data.frame type for data")
  }
  
  m <- model.frame(formula, data)
  mx <- model.matrix(formula, data)
  
  lower_frailty = tolower(frailty)
  
  frailty = switch(lower_frailty, "gamma" = "Gamma", "lognormal" = "LogN", "invgauss" = "InvGauss", "pvf" = "PVF",
                   stop("Invalid frailty specified, please check the frailty input"))
  
  out_frailty = switch(frailty, "Gamma" = "Gamma", "LogN" = "Log-Normal", "InvGauss" = "Inverse Gaussian", "PVF" = "PVF")
  
  if (ncol(m[[1]]) == 2) {
    
    cluster_id <- grep("^cluster\\(", colnames(mx))
    event_id <- grep("^event\\(", colnames(mx))
    
    
    if (length(cluster_id) == 0 && length(event_id) == 0) {
      
      type = "Cluster"
      mx1 = mx[, -c(1), drop = FALSE]
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
      
      neworder = order(y, decreasing = TRUE)
      newrank = seq(1, N, 1)[order(neworder)]
      
      y = y[neworder]
      X = X[neworder, ]
      d = d[neworder]
      newid = newid[neworder]
      
    }
    
    if (length(cluster_id) == 1) {
      
      type = "Cluster"
      pb = unlist(gregexpr('\\(', colnames(mx)[cluster_id])) + 1
      pe = unlist(gregexpr('\\)', colnames(mx)[cluster_id])) - 1
      clsname = substr(colnames(mx)[cluster_id], pb, pe)
      remove_cluster_id = c(which(colnames(mx) == clsname), cluster_id)
      mx1 = mx[, -c(1, remove_cluster_id), drop = FALSE]
      mxid = mx[, cluster_id]
      
      coef_name = colnames(mx1)
      nord = order(mxid)
      mxid = mxid[nord]
      N = length(mxid)
      p = ncol(mx1)
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
      X = mx1[nord, , drop = FALSE]
      d = m[[1]][nord, 2]
      a = max(newid) + 1
      
      neworder = order(y, decreasing = TRUE)
      newrank = seq(1, N, 1)[order(neworder)]
      
      y = y[neworder]
      X = X[neworder, ]
      d = d[neworder]
      newid = newid[neworder]
      
    }
    
    
    if (length(event_id) == 1) {
      
      type = "Multiple"
      pb = unlist(gregexpr('\\(', colnames(mx)[event_id])) + 1
      pe = unlist(gregexpr('\\)', colnames(mx)[event_id])) - 1
      evsname = substr(colnames(mx)[event_id], pb, pe)
      remove_event_id = c(which(colnames(mx) == evsname), event_id)
      mx1 = mx[, -c(1, remove_event_id), drop = FALSE]
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
      X = mx1[nord, , drop = FALSE]
      y = m[[1]][nord, 1]
      d = m[[1]][nord, 2]
      
    }
  }
  
  if (ncol(m[[1]]) == 3) {
    
    type = "Recurrent"
    cluster_id <- grep("^cluster\\(", colnames(mx))
    pb = unlist(gregexpr('\\(', colnames(mx)[cluster_id])) + 1
    pe = unlist(gregexpr('\\)', colnames(mx)[cluster_id])) - 1
    clsname = substr(colnames(mx)[cluster_id], pb, pe)
    remove_cluster_id = c(which(colnames(mx) == clsname), cluster_id)
    mx1 = mx[, -c(1, remove_cluster_id), drop = FALSE]
    mxid = mx[, cluster_id]
    
    coef_name = colnames(mx1)
    nord = order(mxid)
    mxid = mxid[nord]
    N = length(mxid)
    p = ncol(mx1)
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
    X = mx1[nord, , drop = FALSE]
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
    
    initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", maxit = 10, threshold = threshold, type = 1)
    
    if (sum(abs(initGam$coef) > 1e-6) == 0) {
      initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = frailty, maxit = 10, threshold = threshold, type = 1)
    }

    ini = initGam
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    width <- options()$width/2
    len_tune = length(tuneseq)
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMMcal(y, X, d, N, a, newid,
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0, safe.ini = list(coef = ini$coef, est.tht = ini$est.tht, lambda = ini$lambda),
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
      
      progress_width = floor(z/len_tune*width)
      cat('[', paste0(c(rep('=', progress_width), '>', rep('-', width - progress_width)),  collapse=''), ']', round(z/len_tune*100),'%\r')
      
      if (sum(abs(coef0)) < threshold) {
        progress_width = width
        cat('[', paste0(c(rep('=', progress_width), '>', rep('-', width - progress_width)),  collapse=''), ']', round(100),'%\r')
        break
      }
      
      if (p > N) {
        coef0 = ini$coef
        est.tht0 = ini$est.tht
        lambda0 = ini$lambda
      } 
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(est.tht_all)
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
    
    initGam = frailtyMMcal(y, X, d, N, b, NULL, frailty = "Gamma", power = NULL, penalty = NULL, maxit = 10, threshold = tol, type = 2)
    
    if (sum(abs(initGam$coef) > 1e-6) == 0) {
      initGam = frailtyMMcal(y, X, d, N, b, NULL, frailty = frailty, power = NULL, penalty = NULL, maxit = 10, threshold = tol, type = 2)
    }
    
    ini = initGam
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    width <- options()$width/2
    len_tune = length(tuneseq)
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMMcal(y, X, d, N, b, NULL,
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0, safe.ini = list(coef = ini$coef, est.tht = ini$est.tht, lambda = ini$lambda),
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
      
      progress_width = floor(z/len_tune*width)
      cat('[', paste0(c(rep('=', progress_width), '>', rep('-', width - progress_width)),  collapse=''), ']', round(z/len_tune*100),'%\r')
      
      if (p > N) {
        coef0 = ini$coef
        est.tht0 = ini$est.tht
        lambda0 = ini$lambda
      } 
      
      if (sum(abs(coef0)) < 1e-6) {
        progress_width = width
        cat('[', paste0(c(rep('=', progress_width), '>', rep('-', width - progress_width)),  collapse=''), ']', round(100),'%\r')
        break
      }
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(est.tht_all)
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
    
    
    initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = "Gamma", power = NULL, penalty = NULL, maxit = 10, threshold = threshold, type = 3)
    
    if (sum(abs(initGam$coef) > 1e-6) == 0) {
      initGam = frailtyMMcal(y, X, d, N, a, newid, frailty = frailty, power = NULL, penalty = NULL, maxit = 10, threshold = threshold, type = 3)
    }
    
    ini = initGam
    coef0 = ini$coef
    est.tht0 = ini$est.tht
    lambda0 = ini$lambda
    likelihood0 = ini$likelihood
    
    coef_all = list()
    est.tht_all = list()
    lambda_all = list()
    likelihood_all = list()
    BIC_all = list()
    
    width <- options()$width/2
    len_tune = length(tuneseq)
    
    for (z in seq_len(length(tuneseq))) {
      cur = frailtyMMcal(y, X, d, N, a, newid,
                         coef.ini = coef0, est.tht.ini = est.tht0, lambda.ini = lambda0, safe.ini = list(coef = ini$coef, est.tht = ini$est.tht, lambda = ini$lambda),
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
      
      progress_width = floor(z/len_tune*width)
      cat('[', paste0(c(rep('=', progress_width), '>', rep('-', width - progress_width)),  collapse=''), ']', round(z/len_tune*100),'%\r')
      
      if (p > N) {
        coef0 = ini$coef
        est.tht0 = ini$est.tht
        lambda0 = ini$lambda
      } 
      
      if (sum(abs(coef0)) < 1e-6) {
        progress_width = width
        cat('[', paste0(c(rep('=', progress_width), '>', rep('-', width - progress_width)),  collapse=''), ']', round(100),'%\r')
        break
      }
    }
    
    
    coef_all = data.frame(matrix(unlist(coef_all), nrow = length(coef0)))
    est.tht_all = unlist(est.tht_all)
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