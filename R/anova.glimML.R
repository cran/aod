### likelihood ratio test for models of formal class glimML
if(!isGeneric("anova"))
  setGeneric(name = "anova", def = function(object, ...) standardGeneric("anova"))

setMethod(f = "anova", signature(object = "glimML"), definition = function(object, ..., digits = 4){
  List <- as.list(substitute(list(object, ...)))[-1]
  nam <- sapply(List, function(x) deparse(x)[1])
  method <- get(nam[1])@method
# management of methods "BB" and "NB"
  is.meth <- unlist(lapply(nam, function(x){
    fm <- get(x)
    meth <- fm@method
    meth == method
    }))
  nam <- nam[is.meth]
  n <- length(nam)
  if(n < 2)
    stop("At least 2 valid models are needed.")
  dfr <- data.frame("logL" = rep(NA, n), k = rep(NA, n), AIC = rep(NA, n), BIC = rep(NA, n), "Resid. dev." = rep(NA, n), 
                    "Resid. Df" = rep(NA, n), Test = rep(NA, n), Deviance = rep(NA, n), Df = rep(NA, n), 
                    "P(> Chi2)" = rep(NA, n), check.names = FALSE)
  mod <- vector(mode = "character", length = length(nam))
# utility function to remove white spaces at the beginning and at the end of character strings
  tr <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)
  for(i in seq(n)){
    fm <- get(nam[i])
# model name and formula
# paste + tr are necessary when there are interaction terms in the formula
    fb <- paste(tr(deparse(fm@formula)), collapse = " ")
    ft <- deparse(fm@random)
    mod[i] <- paste(nam[i], ": ", "fixed = ", fb, "; random = ", ft, sep = "")
# anova table
    k <- fm@nbpar
    dfres <- df.residual(fm)
    dfr$logL[i] <- fm@logL
    dfr$k[i]    <- k
    dfr$AIC[i]  <- -2 * fm@logL + 2 * k
    dfr$BIC[i]  <- -2 * fm@logL + log(k + dfres) * k
    dfr[i, 5]   <- deviance(fm)
    dfr[i, 6]   <- dfres
    }
  dfr$Test <- c(NA, paste(nam[-nrow(dfr)], nam[-1], sep = "-"))
  dfr$Deviance  <- c(NA, 2 * diff(dfr[,1]))
  dfr$Df <- c(NA, diff(dfr[,2]))
  for(i in 2:n)
    dfr[i , 10] <- 1 - pchisq(abs(dfr$Deviance[i]), df = abs(dfr$Df[i]))
  rownames(dfr) <- nam
  type <- switch(method, BB = "beta-binomial", NB = "negative-binomial")
  new(Class = "anova.glimML", models = mod, anova.table = dfr, type = type, digits = digits)
  })


#anova.glimML <-  function(object, ..., digits = 4){
#    List <- as.list(substitute(list(object, ...)))[-1]
#    nam <- sapply(List, function(x) deparse(x)[1])
#    method <- get(nam[1])@method
## management of methods "BB" and "NB"
#    is.meth <- unlist(lapply(nam, function(x){
#      fm <- get(x)
#      meth <- fm@method
#      meth == method
#      }))
#    nam <- nam[is.meth]
#    n <- length(nam)
#    if(n < 2)
#      stop("At least 2 valid models are needed.")
#    dfr <- data.frame("logL" = rep(NA, n), k = rep(NA, n), AIC = rep(NA, n), BIC = rep(NA, n), "Resid. dev." = rep(NA, n), 
#                      "Resid. Df" = rep(NA, n), Test = rep(NA, n), Deviance = rep(NA, n), Df = rep(NA, n), 
#                      "P(> Chi2)" = rep(NA, n), check.names = FALSE)
#    mod <- vector(mode = "character", length = length(nam))
## utility function to remove white spaces at the beginning and at the end of character strings
#    tr <- function(string) gsub("^[[:space:]]+|[[:space:]]+$", "", string)
#    for(i in seq(n)){
#      fm <- get(nam[i])
## model name and formula
## paste + tr are necessary when there are interaction terms in the formula
#      fb <- paste(tr(deparse(fm@formula)), collapse = " ")
#      ft <- deparse(fm@random)
#      mod[i] <- paste(nam[i], ": ", "fixed = ", fb, "; random = ", ft, sep = "")
## anova table
#      k <- fm@nbpar
#      dfres <- df.residual(fm)
#      dfr$logL[i] <- c(logLik(fm))
#      dfr$k[i]    <- k
#      dfr$AIC[i]  <- AIC(fm)
#      dfr$BIC[i]  <- -2 * c(logLik(fm)) + log(k + dfres) * k
#      dfr[i, 5]   <- deviance(fm)
#      dfr[i, 6]   <- dfres
#      }
#    dfr$Test <- c(NA, paste(nam[-nrow(dfr)], nam[-1], sep = "-"))
#    dfr$Deviance  <- c(NA, 2 * diff(dfr[,1]))
#    dfr$Df <- c(NA, diff(dfr[,2]))
#    for(i in 2:n)
#      dfr[i , 10] <- 1 - pchisq(abs(dfr$Deviance[i]), df = abs(dfr$Df[i]))
#    rownames(dfr) <- nam
#    type <- switch(method, BB = "beta-binomial", NB = "negative-binomial")
#    structure(list(models = mod, anova.table = dfr, type = type, digits = digits), class = "anova.glimML")
#    }

setMethod("show", signature = "anova.glimML", definition = function(object){
  mod <- object@models
  dfr <- object@anova.table
  type <- object@type
  digits <- object@digits
  nam <- rownames(dfr)
  cat("Analysis of Deviance Table (", type, " models)\n\n", sep = "")
  sapply(mod, function(x) cat(x, "\n"))
  cat("\n")
  List <- lapply(dfr, function(x) ifelse(is.na(x), "", format(x, digits = digits)))
  dfr <- as.data.frame(t(do.call("rbind", List)))
  rownames(dfr) <- nam
  print(dfr)
  invisible(dfr)
  })

#### method print for objects of class "anova.glimML"
#print.anova.glimML <-   function(x, ...){
#    mod <- x[["models"]]
#    dfr <- x[["anova.table"]]
#    type <- x[["type"]]
#    digits <- x[["digits"]]
#    nam <- rownames(dfr)
#    cat("Analysis of Deviance Table (", type, " models)\n\n", sep = "")
#    sapply(mod, function(x) cat(x, "\n"))
#    cat("\n")
#    List <- lapply(dfr, function(x) ifelse(is.na(x), "", format(x, digits = digits)))
#    dfr <- as.data.frame(t(do.call("rbind", List)))
#    rownames(dfr) <- nam
#    print(dfr, ...)
#    invisible(dfr)
#    }
