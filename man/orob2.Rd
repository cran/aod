\name{orob2}
\alias{orob2}
\docType{data}
\title{ Germination Data }

\description{
  \dQuote{A 2 x 2 factorial experiment comparing 2 types of seed and 2 root extracts. There are 5 or 6 replicates in each 
  of the 4 treatment groups, and each replicate comprises a number of seeds varying between 4 and 81. The response 
  variable is the proportion of seeds germinating in each replicate.} (Crowder, 1978, Table 3).
  }

\usage{data(orob2)}

\format{
  A data frame with 21 observations on the following 4 variables.
  \describe{
    \item{seed}{a factor with 2 levels: \code{O73} and \code{O75}.}
    \item{root}{a factor with 2 levels \code{BEAN} and \code{CUCUMBER}.}
    \item{n}{a numeric vector: the number of seeds exposed to germination.}
    \item{y}{a numeric vector: the number of seeds which actually germinated.}
  }
}

\references{
  Crowder, M.J., 1978. \emph{Beta-binomial anova for proportions}. Appl. Statist. 27, 34-37.
  }

\keyword{datasets}
