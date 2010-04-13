icc <- function(n, y, data, method = c("REML", "ML"), R = NULL){
  CALL <- match.call()
  method <- match.arg(method)
  require(nlme)
  resp <- c(deparse(substitute(n)), deparse(substitute(y)))
  datan <- data[ , resp]
  if(any(datan[ , 2] > datan[ , 1]))
    stop("Some ", deparse(substitute(y)), " were > ", deparse(substitute(n)), ".")
  if(any(datan[ , 1] <= 0))
    warning("Data with ", deparse(substitute(n)), " <= 0 were discarded.")
  names(datan) <- c("n", "y")
  datan <- datan[datan$n > 0, ]                                                          
  databin <- splitbin(cbind(y, n - y) ~ 1, datan)$tab
# REML/ML method
  resp <- names(databin)[2]
  f <- as.formula(paste(resp, "~ 1"))
  fm <- lme(fixed = f, random = ~ 1 | idbin, data = databin, method = method)
  nu <- attr(fm$apVar, "Pars")
  varnu <- fm$apVar
  rho <- as.numeric(exp(2 * nu[1]) / (exp(2 * nu[1]) + exp(2 * nu[2])))
  # delta method approximate of Var[rho]
  drho <- c(2 * nu[1] * nu[2]^2 / (nu[1]^2 + nu[2]^2)^2, -2 * nu[1]^2 * nu[2] / (nu[1]^2 + nu[2]^2)^2)
  varrho <- as.numeric(t(drho) %*% varnu %*% drho)
# ANOVA method
  k <- nrow(datan) ; n <- datan$n ; N <- sum(n) ; nA <- (N - sum(n^2) / N) / (k - 1) ; p <- sum(datan$y) / N
  f <- as.formula(paste(resp, "~ factor(idbin)"))
  lfm <- lm(f, data = databin)
  a <- anova(lfm)
  f <- c(F.value = a[1, 4], df.num = a[1, 1], df.denom = a[2, 1], P = a[1, 5])
  rho.anova <- (f[1] - 1) / (f[1] + nA - 1)
  varrho.anova <- NA  
  rho <- c(rho, rho.anova)
  names(rho) <- c(method, "ANOVA")
  varrho <- c(varrho, varrho.anova)
  names(varrho) <- c(method, "ANOVA")
# MC estimation
  rho.MC <- numeric(0)
  if(!is.null(R)){
    require(MASS)
    nu.MC <- mvrnorm(n = R, mu = attr(fm$apVar, "Pars"), Sigma = fm$apVar)
    rho.MC <- exp(2 * nu.MC[, 1]) / (exp(2 * nu.MC[, 1]) + exp(2 * nu.MC[, 2]))
    }
# outputs
  class(fm) <- "list"
  features <- c(N = k, n = sum(n), y = sum(datan$y))
  new(Class = "icc", CALL = CALL, rho = rho, f = f, rho.MC = rho.MC, varrho = varrho, method = method, features = features)
  }    

# show
setMethod(f = "show", signature = "icc", definition = function(object){                   
  cat("\nIntra-cluster correlation\n")
  cat("-------------------------\n")
  print(object@CALL)
  feat <- object@features
  cat("N = ", feat["N"], " clusters, n = ", feat["n"], " subjects, y = ", feat["y"], " cases.\n\n", sep = "")
  summry <- data.frame(rho = object@rho, se = object@varrho^0.5)
  R <- length(object@rho.MC)
  if(R > 0){
    summry$CI.low <- c(quantile(object@rho.MC, probs = 0.025), NA)
    summry$CI.up  <- c(quantile(object@rho.MC, probs = 0.975), NA)
    }
  
  List <- lapply(summry, function(x) ifelse(is.na(x), "", format(round(x, digits = 4), nsmall = 3)))
  summ <- as.data.frame(t(do.call("rbind", List)))
  rownames(summ) <- c(object@method, "ANOVA")
  print(summ)
  if(R > 0)
    cat("\n(Monte Carlo 95%CI based on", R, "replications)\n")
  cat("\nF test for ML/REML method (H0: rho = 0):\n")
  cat("F = ", round(object@f[1], digits = 1), ", ", sep = "")
  cat("df num. = ", object@f[2], ", ", sep = "")
  cat("df denom. = ", object@f[3], ", ", sep = "")
  cat("P(> F) = ", round(object@f[4], digits = 4), "\n", sep = "")
  invisible(summry)
  })
