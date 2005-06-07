\name{rats}
\alias{rats}
\docType{data}
\title{ Rats Diet Experiment }
\description{\dQuote{Weil (1970) in Table 1 gives the results from an experiment comprising two treatments. One group of 
  16 pregnant female rats was fed a control diet during pregnancy and lactation, the diet of a second group of 16 pregnant 
  females was treated with a chemical. For each litter the number \eqn{n} of pups alive at 4 days and the number \eqn{x} 
  of pups that survived the 21 day lactation period were recorded.} (Williams, 1975, p. 951).}
\usage{data(rats)}
\format{A data frame with 32 observations on the following 3 variables.
  \describe{
    \item{group}{A factor with levels \code{CTRL} and \code{TREAT}}
    \item{n}{A numeric vector: the number of pups alive at 4 days.}
    \item{y}{A numeric vector: the number of pups that survived the 21 day lactation.}
    }
  }
\source{
  Williams, D.A., 1975. \emph{The analysis of binary responses from toxicological experiments
  involving reproduction and teratogenicity}. Biometrics 31, 949-952.
  }
\references{
  Weil, C.S., 1970. \emph{Selection of the valid number of sampling units and a consideration of their combination 
  in toxicological studies involving reproduction, teratogenesis or carcinogenesis}. Fd. Cosmet. Toxicol. 8, 177-182.
  }
\keyword{datasets}
