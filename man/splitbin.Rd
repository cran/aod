\name{splitbin}
\alias{splitbin}
\title{Splits Binomial Data into Bernoulli Data}
\description{
  Splits binomial data (n, y) and optional covariates into Bernoulli data (0 or 1) and the corresponding covariates.
  }
\usage{
  splitbin(formula, data, id = "id")
  }
\arguments{
  \item{formula}{A formula: its left-hand side must be of the form \code{cbind(y, n - y)} where the modelled 
    probability is \code{y/n}. The right-hand side indicates the covariates which are considered in the ungrouped 
    data set.}
  \item{data}{A data frame where all the variables described in \code{n}, \code{y} and \code{formula} are found.}
  \item{id}{An optional character string naming the identifier (= grouping factor). Default to \dQuote{id}.}
  }
\details{
  Any variable of mode character is converted into a factor before data processing.\cr
  Each line (with \eqn{n > 0}) of the initial data set generates \eqn{n} lines in the final data set with the same 
  covariate pattern, plus an identifier (named \dQuote{id}, by default). The value of this identifier is the row 
  name from which it was expanded in the initial data set.
  }
\value{
  A data frame with response and covariates as indicated in the formula, plus the identifier.
  }
\examples{
  data(orob1)
  splitbin(cbind(y, n - y) ~ 1, orob1)
  splitbin(cbind(y, n - y) ~ dilution, orob1)
  }
\keyword{datagen}
