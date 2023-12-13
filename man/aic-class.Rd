\name{aic-class}
\docType{class}
\alias{aic-class}
\title{Representation of Objects of Formal Class "aic"}
\description{Representation of the output of function \code{AIC}.}
\section{Slots}{%
  \describe{%
    \item{istats}{A data frame with 3 columns describing the models indicated by the row names:
      \describe{
        \item{\bold{df}}{number of parameters in the model},
        \item{\bold{AIC}}{Akaike information criterion for the model (see \code{\link[aod:AIC]{AIC}})},
        \item{\bold{AICc}}{small-sample corrected Akaike information criterion for the model (see \code{\link[aod:AIC]{AIC}}).}
        }
      }
    }
  }

\section{Methods}{%
  \describe{%
    \item{summary}{\code{signature(object = "aic")}}
    \item{show}{\code{signature(object = "aic")}}
    }
  }

\keyword{classes}
