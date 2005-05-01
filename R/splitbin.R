splitbin <- function(formula, data, id = "id"){
# convert any variable of mode character into a factor
  data <- as.data.frame(lapply(data, function(x) if(mode(x) == "character") factor(x) else x))

# check formula
  if(length(formula) != 3)                                                            
    stop(paste(deparse(formula), collapse = " "), "is not a valid formula.")                          
  else                                                                          
    if(substring(deparse(formula)[1], 1, 5) != "cbind")                          
      stop(paste(deparse(formula), collapse = ""), " is not a valid formula.\n",               
           "The response must be a matrix of the form cbind(success, failure)")

# check data w.r.t. formula
  f <- all.vars(formula)
  v <- is.element(f, names(data))
  if(any(!v)){
    err <- paste(f[!v], collapse = ", ")
    stop("Variable(s) ", dQuote(err), " in formula were not found in ", dQuote(deparse(substitute(data))))
    }

# build the data set
  dfr <- model.frame(formula = formula, data = data)
  rownam <- rownames(data)
  resp <- model.response(dfr)
  n <- rowSums(resp)
  y <- resp[ , 1]

# check the response
  if(any(n == 0)){
    err <- rownam[n == 0]
    rownam <- rownam[n > 0]
    y <- y[n > 0]
    dfr <- dfr[n > 0, ]
    n <- n[n > 0]   ### beware: condition must be applied on "y" and "dfr" BEFORE "n"
    warning("Lines with row names ",
            paste(sapply(err, function(x) dQuote(x)), collapse = ", "),
            " were discarded (0 observation).\n")
    }
  if(any(y > n))
    stop("Some values of ", dQuote(f[1]), " were > ", dQuote(f[2]), ".")
  nc <- length(f) - 2
  if(nc > 0)
    dfr <- dfr[ , -1, drop = FALSE]
  List <- vector(mode = "list", length = length(n))
  for(i in seq(length(n))){
    dat <- data.frame(id = rep(rownam[i], each = n[i]),
                      resp = rep(c(0, 1), times = c(n[i] - y[i], y[i])))
    names(dat) <- c(id, f[1])
    List[[i]] <- if(nc == 0) dat else merge(dat, dfr[i, , drop = FALSE], all = TRUE)
    }
  res <- do.call("rbind", List)
  rownames(res) <- seq(nrow(res))
  res
  }
