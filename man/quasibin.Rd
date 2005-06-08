\name{quasibin}
\alias{quasibin}
\title{Quasi-Likelihood Model for Proportions}
\description{The function fits the generalized linear model \dQuote{II} proposed by Williams (1982) accounting 
  for overdispersion in clustered binomial data \eqn{(n, y)}.}
\usage{quasibin(formula, data, link = "logit", phi = NULL, tol =  0.001)}
\arguments{
  \item{formula}{A formula for the fixed effects. The left-hand side of the formula must be of the form 
    \code{cbind(y, n - y)} where the modelled probability is \code{y/n}.}
  \item{link}{The link function for the mean \eqn{p}: \dQuote{logit} or \dQuote{cloglog}.}
  \item{data}{A data frame containing the response (\code{n} and \code{y}) and explanatory variable(s).}
  \item{phi}{When \code{phi} is NULL (the default), the overdispersion parameter \eqn{\phi} is estimated from the data. 
    Otherwise, its value is considered as fixed.}
  \item{tol}{A positive scalar (default to 0.001). The algorithm stops at iteration \eqn{r + 1} when the condition
  \eqn{\chi{^2}[r+1] - \chi{^2}[r] <= tol} is met by the \eqn{\chi^2}{chi-squared} statistics .}
  }
\details{
  For a given cluster \eqn{(n, y)}, the model is:
  \deqn{y~|~\lambda \sim Binomial(n,~\lambda)}{y | \lambda ~ Binomial(n, \lambda)}
  with \eqn{\lambda} a random variable of mean \eqn{E[\lambda] = p}
  and variance \eqn{Var[\lambda] = \phi * p * (1 - p)}.\cr
  The marginal mean and variance are:
  \deqn{E[y] = p}
  \deqn{Var[y] = p * (1 - p) * [1 + (n - 1) * \phi]}
  The overdispersion parameter \eqn{\phi} corresponds to the intra-cluster correlation coefficient,
  which is here restricted to be positive.\cr
  The function uses the function \code{glm} and the parameterization: \eqn{p = h(X b) = h(\eta)}, where \eqn{h} is the
  inverse of a given link function, \eqn{X} is a design-matrix, \eqn{b} is a vector of fixed effects and \eqn{\eta = X b}
  is the linear predictor. \cr
  The estimate of \eqn{b} maximizes the quasi log-likelihood of the marginal model.
  The parameter \eqn{\phi} is estimated with the moment method or can be set to a constant
  (a regular \emph{glim} is fitted when \eqn{\phi} is set to zero). The literature recommends to estimate \eqn{\phi}
  from the saturated model. Several explanatory variables are allowed in \eqn{b}. None is allowed in \eqn{\phi}.
  }
\value{An object of formal class \dQuote{glimQL}: see \code{\link{glimQL-class}} for details.}
\references{
  Moore, D.F., 1987, \emph{Modelling the extraneous variance in the presence of extra-binomial variation}. 
    Appl. Statist. 36, 8-14.\cr
  Williams, D.A., 1982, \emph{Extra-binomial variation in logistic linear models}. Appl. Statist. 31, 144-148.
  }
\author{Matthieu Lesnoff \email{matthieu.lesnoff@cirad.fr}, Renaud Lancelot \email{renaud.lancelot@cirad.fr} }
\seealso{\code{\link{glm}}, \code{\link[geepack]{geese}} in the contributed package \pkg{geepack}, 
  \code{\link[dispmod]{glm.binomial.disp}} in the contributed package \pkg{dispmod}.}

 }
\examples{
  data(orob2) 
  fm1 <- glm(cbind(y, n - y) ~ seed * root,
             family = binomial, data = orob2)
  fm2 <- quasibin(cbind(y, n - y) ~ seed * root,
                  data = orob2, phi = 0)
  fm3 <- quasibin(cbind(y, n - y) ~ seed * root,
                  data = orob2)
  rbind(fm1 = coef(fm1), fm2 = coef(fm2), fm3 = coef(fm3))
  # show the model
  fm3
  # dispersion parameter and goodness-of-fit statistic
  c(phi = fm3@phi, 
    X2 = sum(residuals(fm3, type = "pearson")^2))
  # model predictions
  predfm1 <- predict(fm1, type = "response", se = TRUE)
  predfm3 <- predict(fm3, type = "response", se = TRUE)
  New <- expand.grid(seed = levels(orob2$seed),
                     root = levels(orob2$root))
  predict(fm3, New, se = TRUE, type = "response")
  data.frame(orob2, p1 = predfm1$fit, 
                    se.p1 = predfm1$se.fit,
                    p3 = predfm3$fit,
                    se.p3 = predfm3$se.fit)
  fm4 <- quasibin(cbind(y, n - y) ~ seed + root,
                  data = orob2, phi = fm3@phi)
  # Pearson's chi-squared goodness-of-fit statistic
  # compare with fm3's X2
  sum(residuals(fm4, type = "pearson")^2)
  }
\keyword{regression}
