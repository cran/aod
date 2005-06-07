if(!isGeneric("summary"))
  setGeneric(name = "summary", def = function(object, ...) standardGeneric("summary"))

### summary method for glimML objects
setMethod("summary", signature = "glimML", 
  function(object){
# function to insert NAs in se when there are NAs in b
    insNA <- function(b, se){
      if(any(is.na(b))){
        nb <- length(b)
        SE <- rep(NA, nb)
        j <- 1
        for(i in seq(nb))
          if(!is.na(b[i])){
            SE[i] <- se[j]
            j <- j + 1
            }
        se <- SE
        }
      se
      }
# check whether any fixed-effect coef was set to a fixed value
    nb <- length(coef(object))
    param <- object@param
    vpar <- object@varparam
# position of fixed-effect coefficients
    pos1 <- seq(nb)
# position of parameters set to a fixed value, if any
    fp <- match("fixpar", table = names(object@CALL))
    pos2 <- NA
    if(!is.na(fp))
      pos2 <- eval(object@CALL$fixpar[[2]])
# remove fixed parameters, if any
    pos3 <- setdiff(pos1, pos2)
    Coef <- data.frame()
# compute new var-cov mat, coef vector and position of term(s) to be tested
    if(length(pos3) > 0){
      vp3 <- as.matrix(vpar[pos3, pos3])
      b3 <- param[pos3]
# coef, se, z and t test
      se3 <- sqrt(diag(vp3))
      se3 <- insNA(b3, se3)
      Coef <- data.frame(b = b3, se = se3, z = b3 / se3, P = 2 * (1 - pnorm(abs(b3) / se3)))
      nam <- names(b3)
      rownames(Coef) <- nam
      colnames(Coef) <- c("Estimate", "Std. Error", "z value", "Pr(> |z|)")
      }
# fixed-effect coefficients which were set to a fixed value, if any
    pos4 <- setdiff(pos1, pos3)
    FixedCoef <- data.frame()
    if(length(pos4) > 0){
      FixedCoef <- data.frame(Value = param[pos4])
      }

# overdispersion coefficients phi
# need to check whether any coef was set to a fixed value
    nb <- length(coef(object))
    param <- object@param
    vpar <- object@varparam
# position of overdispersion parameters
    pos1 <- (nb + 1):length(param)
# position of parameters set to a fixed value, if any
    fp <- match("fixpar", table = names(object@CALL))
    pos2 <- NA
    if(!is.na(fp))
      pos2 <- eval(object@CALL$fixpar[[2]])
# remove fixed parameters, if any
    pos3 <- setdiff(pos1, pos2)
    Phi <- data.frame()
# compute new var-cov mat, coef vector and position of term(s) to be tested
    if(length(pos3) > 0){
      vp3 <- as.matrix(vpar[pos3, pos3])
      b3 <- param[pos3]
      se3 <- sqrt(diag(vp3))
      se3 <- insNA(b3, se3)
# coef, se, z and t test for phi
# beware: unilateral test for phi because it cannot be negative
      if(any(b3 < 0))
        warning("Negative values for phi.")
      Phi <- data.frame(b = b3, se = se3, z = b3 / se3, P = 1 - pnorm(abs(b3) / se3))
      nam <- names(b3)
      rownames(Phi) <- nam
      colnames(Phi) <- c("Estimate", "Std. Error", "z value", "Pr(> z)")
      }
# print random coefficients which were set to a fixed value, if any
    pos4 <- setdiff(pos1, pos3)
    FixedPhi <- data.frame()
    if(length(pos4) > 0){
      cat("\nOverdispersion coefficients set to fixed values:\n")
      FixedPhi <- data.frame(Value = param[pos4])
      }
    res <- new(Class = "summary.glimML",
               object = object, Coef = Coef, FixedCoef = FixedCoef, Phi = Phi, FixedPhi = FixedPhi)
    res
    })

### show method for objects of class summary.glimML
setMethod("show", signature = "summary.glimML", 
  function(object){
  Object <- object@object
# Title and call
    switch(Object@method,
      BB = cat("Beta-binomial model\n", "-------------------\n", sep = ""),
      NB = cat("Negative-binomial model\n", "-----------------------\n", sep = ""))
    print(Object@CALL)

# Checks whether convergence problems occurred
    n <- Object@code
    iter <- Object@iterations
    if(n < 3)
      cat("\nConvergence was obtained after " , iter, " iterations.\n", sep = "")
    else
      cat("Possible convergence problem. Optimization process code:", n,"(see ?optim).\n")

# Print estimated fixed effects, if any
    Coef <- object@Coef
    if(nrow(Coef) > 0){
      nam <- rownames(Coef)
      List <- vector(mode = "list", length = 4)
      for(i in 1:4){
        x <- Coef[,i]
        List[[i]] <- if(i < 4)
                       format(round(x, 4), nsmall = 3)
                     else
                       ifelse(x < 1e-4, "< 1e-4", format(round(x, 4), nsmall = 3))
        }
      Coeftext <- as.data.frame(t(do.call("rbind", List)))
      rownames(Coeftext) <- nam
      colnames(Coeftext) <- c("Estimate", "Std. Error", "z value", "Pr(> |z|)")
      cat("\nFixed-effect coefficients:\n")
      print(Coeftext)
      }
# print fixed-effect coefficients which were set to a fixed value, if any
    FixedCoef <- object@FixedCoef
    if(nrow(FixedCoef) > 0){
      cat("\nFixed-effect coefficients set to fixed values:\n")
      print(FixedCoef)
      }

# overdispersion coefficients phi
    Phi <- object@Phi
# compute new var-cov mat, coef vector and position of term(s) to be tested
    if(nrow(Phi) > 0){
      nam <- rownames(Phi)
      List <- vector(mode = "list", length = 4)
      for(i in 1:4){
        x <- Phi[,i]
        List[[i]] <- if(i < 4)
                       format(round(x, 4), nsmall = 3)
                     else
                       ifelse(x < 1e-4, "< 1e-4", format(round(x, 4), nsmall = 3))
        }
      Phitext <- as.data.frame(t(do.call("rbind", List)))
      rownames(Phitext) <- nam
      colnames(Phitext) <- c("Estimate", "Std. Error", "z value", "Pr(> z)")
      cat("\nOverdispersion coefficients:\n")
      print(Phitext)
      }
# print overdispersion coefficients which were set to a fixed value, if any
    FixedPhi <- object@FixedPhi
    if(nrow(FixedPhi) > 0){
      cat("\nOverdispersion coefficients set to fixed values:\n")
      print(FixedPhi)
      }

    cat("\nLog-likelihood = ", format(round(Object@logL, 3), nsmall = 3), "; ", sep = "")
    cat("nbpar = ", Object@nbpar, "; ", sep = "")
    cat("df.residual = ", df.residual(Object), "; ", sep = "")
    cat("Deviance = ", format(round(deviance(Object), 3), nsmall = 3), "; ", sep = "")
    cat("AIC = ", format(round(-2 * Object@logL + 2 * Object@nbpar, 3), nsmall = 3), "\n\n", sep = "")
    invisible(object)
    })


### show method for glimML objects
setMethod("show", signature = "glimML",  function(object) show(summary(object)))
