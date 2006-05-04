negbin <- function(formula, random, data, phi.ini = NULL, warnings = FALSE, na.action = na.omit,
                   fixpar = list(), hessian = TRUE, ...){
  CALL <- mf <- match.call(expand.dots = FALSE)

# formula for the correlation parameters: random
  if(length(random) == 3){
    form <- deparse(random)
    warning("The formula for rho (", form, ") contains a response which is ignored.")
    random <- random[-2]
    }
  explain <- as.character(attr(terms(random), "variables"))[-1]
  if(length(explain) > 1){
    warning("The formula for rho contains several explanatory variables (", paste(explain, collapse = ", "), ").\n",
            "Only the first one (", explain[1], ") was considered.")
    explain <- explain[1]
    }

# if data is not given in the call, a global formula must be built from fixed and random formulas
# to build the data frame which will be returned by the function

# global formula
  gf3 <- if(length(explain) == 1) paste(as.character(formula[3]), explain, sep = " + ") else as.character(formula[3])
  gf <- formula(paste(formula[2], "~", gf3))

# get the data
  if(missing(data)) data <- environment(gf)
  
# model frame and model matrix for the fixed effects
  mb <- match(c("formula", "data", "na.action"), names(mf), 0)
  mfb <- mf[c(1, mb)]
  mfb$drop.unused.levels <- TRUE
  mfb[[1]] <- as.name("model.frame")
  names(mfb)[2] <- "formula"
  mfb <- eval(mfb, parent.frame())
  mt <- attr(mfb, "terms")
  y <- model.response(mfb, "numeric")
  modmatrix.b <- if(!is.empty.model(mt)) model.matrix(mt, mfb, contrasts) else matrix(, NROW(y), 0)
  offset <- model.offset(mfb)

### change on 12th July 2005 (check lines with weight = 0)

  if(any(is.infinite(offset)))
    warning("The data set contains at least one line with weight = 0.\n")

### end change

# model frame and model matrix for the correlation structure
  mr <- match(c("random", "data", "na.action"), names(mf), 0)
  mr <- mf[c(1, mr)]
  mr$drop.unused.levels <- TRUE
  mr[[1]] <- as.name("model.frame")
  names(mr)[2] <- "formula"
  mr <- eval(mr, parent.frame())
  if(length(explain) == 0)
    modmatrix.phi <- model.matrix(object = ~ 1, data = mr)
  else{
    express <- paste("model.matrix(object = ~ -1 + ", explain, ", data = mr", 
                     ", contrasts = list(", explain, " = 'contr.treatment'))", sep = "")
    if(is.ordered(data[ , match(explain, table = names(mr))]))
      warning(explain, " is an ordered factor.\n", "Treatment contrast was used to build model matrix for phi.")
    modmatrix.phi <- eval(parse(text = express))
    }

  fm <- glm(formula = formula, family = poisson, data = data)
  b <- coef(fm)
  if(any(is.na(b))){
    print(nab <- b[is.na(b)])
    stop("Initial values for the fixed effects contain at least one missing value.")
    }
  nb.b <- ncol(modmatrix.b)
  nb.phi <- ncol(modmatrix.phi)
  if(is.null(phi.ini)) phi.ini <- rep(0.1, nb.phi)
  param.ini <- c(b, phi.ini) 
  if(!is.null(unlist(fixpar)))
    param.ini[fixpar[[1]]] <- fixpar[[2]]  

  # minuslogL
  minuslogL <- function(param){
    if(!is.null(unlist(fixpar)))
      param[fixpar[[1]]] <- fixpar[[2]]  
    b <- param[1:nb.b]
    eta <- as.vector(modmatrix.b %*% b)
    mu <- if(is.null(offset)) exp(eta) else exp(offset + eta)
    phi <- as.vector(modmatrix.phi %*% param[(nb.b + 1):(nb.b + nb.phi)])
    k <- 1 / phi 

## old code that caused problems with R 2.3.0 (ifelse producing many warnings ==> slow execution)
#    fn <- ifelse(phi == 0, dpois(x = y, lambda = mu, log = TRUE), 
#                           lgamma(y + k) - lfactorial(y) - lgamma(k) + k * log(k / (k + mu)) + y * log(mu / (k + mu)))
#    fn <- sum(fn)

# new code (2006-05-04)
  cnd <- phi == 0
  f1 <- dpois(x = y[cnd], lambda = mu[cnd], log = TRUE) 
  y2 <- y[!cnd]; k2 <- k[!cnd]; mu2 <- mu[!cnd]
  f2 <- lgamma(y2 + k2) - lfactorial(y2) - lgamma(k2) + k2 * log(k2 / (k2 + mu2)) + y2 * log(mu2 / (k2 + mu2))
  fn <- sum(c(f1, f2))
# end new code      

    if(!is.finite(fn))
      fn <- -1e20
    -fn
    }

# Fit
  withWarnings <- function(expr){
    myWarnings <- NULL
    wHandler <- function(w){
      myWarnings <<- c(myWarnings, list(w))
      invokeRestart("muffleWarning")
      }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
    }
  reswarn <- withWarnings(optim(par = param.ini, fn = minuslogL, hessian = hessian, ...))
  res <- reswarn$value
  if(warnings){
    if(length(reswarn$warnings) > 0){
      v <- unlist(lapply(reswarn$warnings, as.character))
      tv <- data.frame(message = v, freq = rep(1, length(v)))
      cat("Warnings during likelihood maximisation:\n")
      print(aggregate(tv[, "freq", drop = FALSE], list(warning = tv$message), sum))
      }
    }

# Results
  param <- res$par
  namb <- colnames(modmatrix.b)
  namphi <- paste("phi", colnames(modmatrix.phi), sep = ".")
  nam <- c(namb, namphi)
  names(param) <- nam

  if(!is.null(unlist(fixpar)))
    param[fixpar[[1]]] <- fixpar[[2]]
  varparam <- NA
  if(hessian){
    H <- res$hessian
    if(is.null(unlist(fixpar)))
      varparam <- qr.solve(H)
    else{
      idparam <- 1:(nb.b + nb.phi)
      idestim <- idparam[-fixpar[[1]]]
      Hr <- H[-fixpar[[1]], -fixpar[[1]]]
      Vr <- solve(Hr)
#      varparam <- H
#      varparam[idestim, idestim] <- Vr
      varparam <- matrix(rep(NA, NCOL(H) * NROW(H)), ncol = NCOL(H))
      varparam[idestim, idestim] <- Vr
      }
    }
  if(!is.null(dim(varparam)))
    dimnames(varparam) <- list(nam, nam)
    
  nbpar <- if(is.null(unlist(fixpar)))
             sum(!is.na(param))
           else
             sum(!is.na(param[-fixpar[[1]]]))
  logL.max <- sum(dpois(x = y, lambda = y, log = TRUE))
  logL <- -res$value
  dev <- -2 * (logL - logL.max)
  df.residual <- length(y) - nbpar
  iterations <- res$counts[1]
  code <- res$convergence

# Output
  new(Class = "glimML", CALL = CALL, link = "log", method = "NB", data = data, formula = formula, random = random, 
      param = param, varparam = varparam, logL = logL, logL.max = logL.max, dev = dev, df.residual = df.residual, 
      nbpar = nbpar, iterations = iterations, code = code, param.ini = param.ini, na.action = na.action)
  }
