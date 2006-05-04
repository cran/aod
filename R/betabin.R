### fitting function
betabin <- function(formula, random, data = NULL, link = c("logit", "cloglog"), phi.ini = NULL, warnings = FALSE, 
                    na.action = na.omit, fixpar = list(), hessian = TRUE, ...){
# get the call
  CALL <- mf <- match.call(expand.dots = FALSE)

# utility function to remove white spaces at the beginning and at the end of character strings
  tr <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)

# check the link
  link <- match.arg(link)

# formula for the fixed effects: formula
  if(length(formula) != 3)                                                            
    stop(paste(tr(deparse(formula)), collapse = " "), "is not a valid formula.")                          
  else                                                                          
    if(substring(deparse(formula)[1], 1, 5) != "cbind")                          
      stop(paste(tr(deparse(formula)), collapse = ""), " is not a valid formula.\n",               
           "The response must be a matrix of the form cbind(success, failure)")

# formula for the correlation parameters: random
  if(length(random) == 3){
    form <- deparse(random)
    warning("The formula for phi (", form, ") contains a response which is ignored.")
    random <- random[-2]
    }
  explain <- as.character(attr(terms(random), "variables"))[-1]
  if(length(explain) > 1){
    warning("The formula for phi contains several explanatory variables (", paste(explain, collapse = ", "), ").\n",
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
  Y <- model.response(mfb, "numeric")
  modmatrix.b <- if(!is.empty.model(mt)) model.matrix(mt, mfb, contrasts) else matrix(, NROW(Y), 0)
  weights <- model.weights(mfb)
  if(!is.null(weights) && any(weights < 0))
    stop("Negative wts not allowed")

### change on 12th July 2005 (check lines with weight = 0)
  n <- rowSums(Y)
  y <- Y[, 1]
  
  if(any(n == 0))
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

# Data
  fam <- eval(parse(text = paste("binomial(link =", link,")")))
  fm <- glm(formula = formula, family = fam, data = data, na.action = na.action)

# initial b are the ML estimates of the pure binomial model
  b <- coef(fm)
  if(any(is.na(b))){
    print(nab <- b[is.na(b)])
    stop("Initial values for the fixed effects contain at least one missing value.")
    }
  
## Initial values
  nb.b <- ncol(modmatrix.b)
  nb.phi <- ncol(modmatrix.phi)
# check phi.ini
  if(!is.null(phi.ini) && !(phi.ini < 1 & phi.ini > 0))
    stop("phi.ini was set to ", phi.ini, ".\nphi.ini should verify 0 < phi.ini < 1")
  else
# intial values for phi.ini
    if(is.null(phi.ini))
      phi.ini <- rep(.1, nb.phi)

  param.ini <- c(b, phi.ini)

  if(!is.null(unlist(fixpar))) 
    param.ini[fixpar[[1]]] <- fixpar[[2]]  
  # minuslogL
  minuslogL <- function(param){
   if(!is.null(unlist(fixpar)))
     param[fixpar[[1]]] <- fixpar[[2]]  
   b <- param[1:nb.b]
   eta <- as.vector(modmatrix.b %*% b)
   p <- invlink(eta, type = link)
   phi <- as.vector(modmatrix.phi %*% param[(nb.b + 1):(nb.b + nb.phi)])

# old code changed on 4th May 2006
#   fn <- ifelse(phi == 0,
#           dbinom(x = y, size = n, prob = p, log = TRUE),
#           lchoose(n,y) + lbeta(p * (1-phi)/phi + y, (1-p) * (1-phi)/phi + n-y) - lbeta(p*(1-phi)/phi, (1-p)*(1-phi)/phi))
#    fn <- sum(fn)

# new code 4th May 2006
   cnd <- phi == 0
   f1 <- dbinom(x = y[cnd], size = n[cnd], prob = p[cnd], log = TRUE) 
   n2 <- n[!cnd] ; y2 <- y[!cnd] ; p2 <- p[!cnd] ; phi2 <- phi[!cnd]
   f2 <- lchoose(n2, y2) + lbeta(p2 * (1 - phi2)/phi2 + y2, (1 - p2) * (1 - phi2)/phi2 + n2 - y2) - lbeta(p2 * (1 - phi2)/phi2, (1 - p2) * (1 - phi2)/phi2)
   fn <- sum(c(f1, f2))       
# end new code 4th May 2006

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

## Results
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
      dimnames(Vr) <- list(idestim, idestim)
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

  logL.max <- sum(dbinom(x = y, size = n, prob = y / n, log = TRUE))
  logL <- -res$value
  dev <- -2 * (logL - logL.max)
  df.residual <- sum(n > 0) - nbpar
  iterations <- res$counts[1]
  code <- res$convergence

# Output
  new(Class = "glimML", CALL = CALL, link = link, method = "BB", data = data, formula = formula, random = random, 
      param = param, varparam = varparam,  logL = logL, logL.max = logL.max, dev = dev, df.residual = df.residual,
      nbpar = nbpar, iterations = iterations, code = code, param.ini = param.ini, na.action = na.action)
  }
  
