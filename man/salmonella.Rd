\name{salmonella}
\alias{salmonella}
\docType{data}
\title{Salmonella Reverse Mutagenicity Assay}
\description{
  \dQuote{Data for our third example were compiled by Margolin et al. (1981) from an Ames \emph{Salmonella}
  reverse mutagenicity assay. Table 1 shows the number of revertant colonies observed on each of 3 replicate 
  plates tested at each of 6 dose levels of quinoline.} (Breslow, 1984, Table 1).
  }
\usage{data(salmonella)}
\format{
  A data frame with 18 observations on the following 2 variables.
  \describe{
    \item{dose}{a numeric vector: the dose level of quinoline (microgram per plate).}
    \item{y}{a numeric vector: the number of revertant colonies of TA98 \emph{Salmonella}.}
    }
  }
\source{
  Breslow, N.E., 1984. \emph{Extra-Poisson variation in log-linear models}. Applied Statistics 33(1), 38-44.
  }
\references{
  Margolin, B.H., Kaplan, N., Zeiger, E., 1981. \emph{Statistical analysis of the Ames \emph{Salmonella} / microsome test}.
  Proc. Natl Acad. Sci. USA 76, 3779-3783.
  }
\keyword{datasets}
