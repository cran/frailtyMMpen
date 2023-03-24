frailtyMMcal <- function(y, X, d, N, a, id, coef.ini = NULL, est.tht.ini = NULL, lambda.ini = NULL, frailty = "Gamma", power = NULL, penalty = NULL, gam.val = NULL, tune = NULL, maxit = 200, threshold = 1e-5, type = 0) {
  
  p = dim(X)[2]
  
  frailtyc = switch(frailty, "Gamma" = 0, "LogN" = 1, "InvGauss" = 2, "PVF" = 3)
  
  if (is.null(penalty)) {
    penaltyc = 0
    gam = 0
    tune = 0
  } else {
    penaltyc = switch(penalty, "LASSO" = 1, "MCP" = 2, "SCAD" = 3)
    gam = switch(penalty, "LASSO" = 0, "MCP" = 3, "SCAD" = 3.7)
  }
  
  if (!is.null(gam.val)) {
    gam = gam.val
  }
  
  if (is.null(power)) {
    power = 0.0
  }
  
  coef = rep(0, p)
  est.tht = 1
  lambda = rep(1/N, N)
  
  if (!is.null(coef.ini)) {
    coef = coef.ini
  }
  
  if (!is.null(est.tht.ini)) {
    est.tht = est.tht.ini
  }
  
  if (!is.null(lambda.ini)) {
    lambda = lambda.ini
  } 
  
  coef.ini = coef
  lambda.ini = lambda
  est.tht.ini = est.tht
  
  l0 = logLikcal(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type)
  error = 3

  ## MM iteration
  num = 0
  SQS1 = 1
  p0 = 0
  threshold = threshold * (p+1)
  
  while(error > threshold && num < maxit) {
    
    coef0 = coef
    est.tht0 = est.tht
    lambda0 = lambda

    p1 = proc.time()[1]
    rs1 = MMroutine(y, X, d, coef, lambda, est.tht, frailtyc, penaltyc, gam, tune, id, N, a, p, power, type)
    p2 = proc.time()[1]
    
    coef1 = rs1$coef
    est.tht = rs1$est.tht
    lambda = rs1$lambda
    Ar = rs1$Ar
    
    if (sum(abs(coef1)) < threshold) {
      break
    }

    u_be = coef1 - coef

    p3 = proc.time()[1]
    rs2 = MMroutine(y, X, d, coef1, lambda, est.tht, frailtyc, penaltyc, gam, tune, id, N, a, p, power, type)
    p4 = proc.time()[1]
    
    coef2 = rs2$coef
    est.tht = rs2$est.tht
    lambda = rs2$lambda

    if (!SQS1) {
      coef = coef2
    }
  
    if (SQS1) {
      
      v_be = coef2 - 2*coef1 + coef
      al_be = sum(u_be*v_be)/sum(v_be^2)
      if (al_be > -1) {al_be = -1}
      
      for (k1 in 0:7) {
        dc = (2*al_be*u_be - al_be^2*v_be) * 2^(-k1)
  
        if (sum(abs(dc)) > p) {
          next
        }
  
        coef.temp = coef - dc
  
        p5 = proc.time()[1]
        rs = MMroutine(y, X, d, coef.temp, lambda, est.tht, frailtyc, penaltyc, gam, tune, id, N, a, p, power, type)
        p6 = proc.time()[1]
  
        if (rs$error) {
          coef = coef.ini
          est.tht = est.tht.ini
          lambda = lambda.ini
          SQS1 = 0
          num = 0
          break
        }
        
        if (!backtrackerror(model = rs, coef = coef.temp, est.tht = est.tht, lambda = lambda)) {
          break
        }
      }
      
      if (!rs$error) {
        coef = rs$coef
        est.tht = rs$est.tht
        lambda = rs$lambda
        Ar = rs$Ar
      }
    }
 
    error = sum(abs(coef - coef0)) + sum(abs(est.tht - est.tht0))
    
    num = num + 1
    # cat(error, " ", est.tht, " ", " ", al_be, " ", num, '\n')
    p0 = p0 + (p2 - p1) + (p4 - p3) + (p6 - p5)
  }
  
  l1 = logLikcal(y, X, d, coef, lambda, est.tht, frailtyc, id, N, a, p, power, type)
  
  output = list(coef = coef,
                est.tht = est.tht,
                lambda = lambda,
                likelihood = l1,
                Ar = Ar,
                input = list(y = y, X = X, d = d),
                iter = num,
                convergence = error)
  
  # cat("->", p0)
  return(output)
}

