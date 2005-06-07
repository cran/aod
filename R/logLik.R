if(!isGeneric("logLik"))
  setGeneric(name = "logLik", def = function(object, ...) standardGeneric("logLik"))

setMethod(f = "logLik", signature = "glimML", definition = function(object, ...){
  val <- object@logL
  attr(val, "df") <- object@nbpar
  class(val) <- "logLik"
  val
  })

if(!isGeneric("AIC"))
  setGeneric("AIC", function(object, ..., k = 2) standardGeneric("AIC"))

setMethod(f = "AIC", signature = "logLik", definition = function(object, ..., k = 2)
  -2 * c(object) + k * attr(object, "df"))

setMethod(f = "AIC", signature = "glimML", definition = function(object, ..., k = 2){
  if(length(list(...))){
    object <- list(object, ...)
    val <- lapply(object, logLik)
    val <- as.data.frame(t(sapply(val, function(el) c(attr(el, "df"), AIC(el, k = k)))))
    names(val) <- c("df", "AIC")
    Call <- match.call()
    Call$k <- NULL
    row.names(val) <- as.character(Call[-1])
    val
    }
  else
    AIC(logLik(object), k = k)
  })
