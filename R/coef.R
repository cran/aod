if(!isGeneric("coef"))
  setGeneric(name = "coef", def = function(object, ...) standardGeneric("coef"))

## model coefficients for class glimQL (functions quasibin and quasipois)
setMethod(f = "coef", signature(object = "glimQL"), definition = function(object, ...) coef(object@fm))

## model coefficients for class glimML (functions betabin and betapois)
setMethod(f = "coef", signature = "glimML", function(object, ...) {
  f <- object@formula[-2]
  mm <- model.matrix(f, data = object@data)
  v <- object@param[seq(ncol(mm))]
  names(v) <- colnames(mm)
  v
  })
