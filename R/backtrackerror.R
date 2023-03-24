backtrackerror <- function(model, coef, est.tht, lambda, lambda2 = NULL) {
  coefn = model$coef
  est.thtn = model$est.tht
  
  c1 = abs((sum(coefn) - sum(coef)) / sum(coef))
  c2 = abs((est.thtn - est.tht) / est.tht)
  
  if (is.null(lambda2) && !is.null(dim(lambda))) {
    lambdan = model$lambda
    c3 = abs((sum(lambdan) - sum(lambda)) / sum(lambda))
  } 
  
  if (is.null(dim(lambda))) {
    lambdan = unlist(model$lambda)
    c3 = abs((sum(lambdan) - sum(unlist(lambda))) / sum(unlist(lambda)))
  }
  
  if (!is.null(lambda2)){
    lambdan1 = model$lambda1
    lambdan2 = model$lambda2
    c3 = max(abs((sum(lambdan1) - sum(lambda)) / sum(lambda)), abs((sum(lambdan2) - sum(lambda2)) / sum(lambda2)))
  }

  if (c1 > length(coef)*2 || c2 > 3 || c3 > 10) {
    return(1)
  }
  return(0)
}