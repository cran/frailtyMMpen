#' @method print fmm_summary
#' @export

print.fmm_summary <- function(x, ...) {
  
  cat("Call:\n")
  dput(attr(x, "call"))
  cat("\n")
  
  cat(paste("  ", x$frailty, "frailty model with", x$datatype, "data, estimated using MM algorithm\n"))
 
  cat("\n")
  
  coef_summary = data.frame(coef = round(x$coef, 3),
                            "exp(coef)" = round(exp(x$coef), 3),
                            "se(coef)" = round(x$coef_sd, 3),
                            z = round(x$zcoef, 3),
                            p = x$pval_coef)
  row.names(coef_summary) = x$coefname
  colnames(coef_summary) = c("coef", "exp(coef)", "se(coef)", "z", "p")
  
  printCoefmat(coef_summary, digits = 3, signif.stars = TRUE)
  
  cat("\n")
  cat(paste("   Frailty parameter, theta:", round(x$est.tht, 3), "(",
            round(x$th_sd, 3), ")", "with p-value = ", sprintf("%.3e", x$pval_th), "\n"))
  cat("\n")
  cat(paste("   Observed likelihood:", round(x$likelihood, 3), "\n"))
  cat("\n")
  cat("   Convergence criteria:")
  cat(paste(" thershold:", sprintf("%.3e", x$convergence), ", number of iterations:", x$iter))

}

#' print a non-penalized regression object
#' 
#' Print the summary of a non-penalized regression fitted by any model with function \code{frailtyMM}
#' 
#' @param x Object with class "fmm" fitted by function \code{frailtyMM}.
#' @param ... Ignored
#' @method print fmm
#' @return No return value, called to print the summary for non-penalized regression.
#' @export
#' @seealso \code{\link{frailtyMM}}
print.fmm <- function(x, ...) {
  
  out = summary(x)
  print(out)
  
}

#' print a penalized regression object
#' 
#' Print the summary of a non-penalized regression fitted by any model with function \code{frailtyMMpen}.
#' The first column is the tuning parameter sequence, the second column is the degree of freedom and the third column is the BIC.
#'
#' @param x Object with class "fpen" fitted by function \code{frailtyMMpen}.
#' @param ... Ignored
#' @method print fpen
#' @return No return value, called to print the summary for penalized regression.
#' @export
#'
#' @seealso \code{\link{frailtyMMpen}}
print.fpen <- function(x, ...) {
  
  cat("Call:\n")
  dput(attr(x, "call"))
  cat("\n")
  
  tune = x$tune
  df = colSums(abs(x$coef) > 1e-5)
  
  tune_summary = data.frame(tune = round(tune, 5),
                            Df = df,
                            BIC = round(x$BIC, 3))
  
  print(tune_summary)
  
}