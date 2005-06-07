\name{lizards}
\alias{lizards}
\docType{data}
\title{A Comparison of Site Preferences of Two Species of Lizard}
\description{
  \dQuote{These data describe the daytime habits of two species of lizards, \emph{grahami} and \emph{opalinus}. 
  They were collected by observing occupied sites or perches and recording the appropriate description, namely 
  species involved, time of the day, height and diameter of the perch and whether the site was sunny or shaded. 
  Time of the day is recorded as early, mid-day or late.} (McCullagh and Nelder, 1989, p.129).
  }
\usage{data(lizards)}
\format{A data frame with 24 observations on the following 6 variables.
  \describe{
    \item{Site}{A factor with levels \code{Sun} and \code{Shade}.}
    \item{Diameter}{A factor with levels \code{D <= 2} and \code{D > 2} (inches). }
    \item{Height}{A factor with levels \code{H < 5} and \code{H >= 5} (feet).}
    \item{Time}{A factor with levels \code{Early}, \code{Mid-day} and \code{Late}.}
    \item{grahami}{A numeric vector giving the observed sample size for \emph{grahami} lizards.}
    \item{opalinus}{A numeric vector giving the observed sample size for \emph{opalinus} lizards.}
    }
  }
\details{The data were originally published in Fienberg (1970).}
\source{McCullagh, P., Nelder, J.A., 1989. \emph{Generalized linear models}. London, Chapman & Hall, 511 p.}
\references{Fienberg, S.E., 1970. \emph{The analysis of multidimensional contingency tables}. Ecology 51: 419-433.}
\examples{data(lizards)}
\keyword{datasets}
