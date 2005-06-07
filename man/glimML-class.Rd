\name{glimML-class}
\docType{class}
\alias{glimML-class}
\alias{show,glimML-method}
\title{Representation of Models of Formal Class "glimML"}
\description{Representation of models of formal class "glimML" fitted by maximum-likelihood method.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("glimML", ...)} or, 
  more commonly, via the \code{betabin} or \code{negbin} functions.}}
\section{Slots}{
  \describe{
    \item{\code{CALL}}{The call of the function.}
    \item{\code{link}}{The link function used to transform the mean: \dQuote{logit}, \dQuote{cloglog} or \dQuote{log}.}
    \item{\code{method}}{The type of fitted model: \dQuote{BB} for beta-binomial and \dQuote{NB} for negative-binomial models.}
    \item{\code{formula}}{The formula used to model the mean.}
    \item{\code{random}}{The formula used to model the overdispersion parameter \eqn{\phi}.}
    \item{\code{data}}{Data set to which model was fitted. Different from the original data in case of missing value(s).}
    \item{\code{param}}{The vector of the ML estimated parameters \eqn{b} and \eqn{\phi}.}
    \item{\code{varparam}}{The variance-covariance matrix of the ML estimated parameters \eqn{b} and \eqn{\phi}.}
    \item{\code{logL}}{The log-likelihood of the fitted model.}
    \item{\code{logL.max}}{The log-likelihood of the maximal model (data).}
    \item{\code{dev}}{The deviance of the model, i.e., \code{- 2 * (logL - logL.max)}.}
    \item{\code{df.residual}}{The residual degrees of freedom of the fitted model.}
    \item{\code{nbpar}}{The number of \bold{estimated} parameters, i.e., nbpar = total number of parameters - number 
      of fixed parameters. See argument \code{fixpar} in \code{\link{betabin}} or \code{\link{negbin}}.}
    \item{\code{iterations}}{The number of iterations performed in \code{optim}.}
    \item{\code{code}}{An integer (returned by \code{optim}) indicating why the optimization process terminated.
      \describe{
        \item{1}{Relative gradient is close to 0, current iterate is probably solution.}
        \item{2}{Successive iterates within tolerance, current iterate is probably solution.}
        \item{3}{Last global step failed to locate a point lower than estimate. Either estimate is an approximate 
                 local minimum of the function or \code{steptol} is too small.}
        \item{4}{Iteration limit exceeded.}
        \item{5}{Maximum step size \code{stepmax} exceeded 5 consecutive times. Either the function is unbounded below, 
                 becomes asymptotic to a finite value from above in some direction or \code{stepmax} is too small.}
        }
      }
    \item{\code{param.ini}}{The initial values provided to the ML algorithm.}
    \item{\code{na.action}}{A function defining the action taken when missing values are encountered.}
    }
  }

\keyword{classes}
