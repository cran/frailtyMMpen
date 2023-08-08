frailtyMMcal <- function(y, X, d, N, a, id, coef.ini = NULL, est.tht.ini = NULL, lambda.ini = NULL, safe.ini = NULL, frailty = "Gamma", power = NULL, penalty = NULL, gam.val = NULL, tune = NULL, maxit = 200, threshold = 1e-5, type = 0, SQS1 = 1) {
  
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
  p0 = 0
  threshold = threshold * (p+1)
  error_count = 0
  
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
    
    if (rs1$error) {
      
      if (error_count >= 1 && !is.null(safe.ini)) {
        
        coef = safe.ini$coef
        coef1 = coef
        est.tht = safe.ini$est.tht
        lambda = safe.ini$lambda
        
      } else {
        
        coef = coef.ini
        coef1 = coef
        est.tht = est.tht.ini
        lambda = lambda.ini
        SQS1 = 0
        num = 0
      }
      error_count = error_count + 1
      
    }
    
    if (sum(abs(coef1)) < threshold) {
      coef = rep(0, p)
      break
    }

    u_be = coef1 - coef

    p3 = proc.time()[1]
    rs2 = MMroutine(y, X, d, coef1, lambda, est.tht, frailtyc, penaltyc, gam, tune, id, N, a, p, power, type)
    p4 = proc.time()[1]
    
    coef2 = rs2$coef
    est.tht = rs2$est.tht
    lambda = rs2$lambda
    
    al_be = -1
    
    if (rs2$error) {

      if (error_count >= 1 && !is.null(safe.ini)) {
        
        coef = safe.ini$coef
        coef2 = coef
        est.tht = safe.ini$est.tht
        lambda = safe.ini$lambda
        
      } else {
        
        coef = coef.ini
        coef2 = coef
        est.tht = est.tht.ini
        lambda = lambda.ini
        SQS1 = 0
        num = 0
        error_count = error_count + 1
      }
      error_count = error_count + 1
      
      
    }

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
          error_count = error_count + 1
          break
          
        }
        
        temp = backtrackerror(model = rs, coef = coef.temp, est.tht = est.tht, lambda = lambda)
        
        if (temp == 2) {
          
          coef = coef.ini
          est.tht = est.tht.ini
          lambda = lambda.ini
          SQS1 = 0
          num = 0
          error_count = error_count + 1
          break
          
        }
        
        if (!temp) {
          break
        }
      }
      
      if (!rs$error && temp != 2) {
        coef = rs$coef
        est.tht = rs$est.tht
        lambda = rs$lambda
        Ar = rs$Ar
      }
    }
 
    error = sum(abs(coef - coef0)) + sum(abs(est.tht - est.tht0))
    
    if (error > 1e+5 && SQS1 == 0) {
      coef = coef.ini
      est.tht = est.tht.ini
      lambda = lambda.ini
      num = 0
      error_count = error_count + 2
    }
    
    if (error_count > 2) {
      break
    }
    
    num = num + 1
    # cat(error, " ", est.tht, " ", " ", al_be, " ", num, " ", sum(abs(coef)), " ", error_count, '\n')
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

