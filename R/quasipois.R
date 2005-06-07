quasipois <- function(formula, data, phi = NULL, tol = 0.001){  
# check call validity
  CALL <- match.call(expand.dots = FALSE)
  mf <- model.frame(formula = formula, data = data)
  y <- model.response(mf)
# computations
  if(is.null(phi)){
    phi <- 1e-04
    X2 <- 0
    delta <- X2 + 2 * tol
    data$w <- rep(1, length(y))
    fm <- glm(formula = formula, family = poisson, data = data, weights = w)                  
    ok <- TRUE                                                                      
    while(ok){
      X2 <- sum(residuals(fm, type = "pearson")^2)
      delta <- X2 - df.residual(fm)
      if(delta <= tol)                                                              
        ok <- FALSE                                                                 
        else{                                                                         
          phi <- phi * sum(residuals(fm, type = "pearson")^2) / df.residual(fm)
          data$w <- 1 / (1 + phi * fitted(fm))
          fm <- update(fm, weights = w)                                               
          }                                                                           
        }
      }        
  else{
    fm <- glm(formula = formula, family = poisson, data = data)
    delta.logL <- 1
    while(delta.logL > 1e-06){                 
      data$w <- 1 / (1 + phi * fitted(fm))    
      fm.new <- update(fm, weights = w)
      delta.logL <- logLik(fm.new) - logLik(fm)
      fm <- fm.new
      }
    }

# results
  data$w <- 1 / (1 + phi * fitted(fm))
  fm <- glm(formula = formula, family = poisson, data = data, weights = w)
# outputs
  new(Class = "glimQL", CALL = CALL, fm = fm, phi = phi)
  }

