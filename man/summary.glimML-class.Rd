\name{summary.glimML-class}
\docType{class}
\alias{summary.glimML-class}
\alias{show,glimML-class}
\alias{summary,glimML-method}
\alias{show,summary.glimML-method}

\title{Summary of Objects of Class "summary.glimML"}

\description{Representation of the summary of an object of formal class \dQuote{glimML}.}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("summary.glimML", ...)} or, 
  more commonly, via the \code{summary} or \code{show} method for objects of formal class \dQuote{glimML}.}

\section{Slots}{
  \describe{
    \item{\code{object}}{An object of formal class \dQuote{glimML}.}
    \item{\code{Coef}}{A data frame containing the estimates, standard error, z value and P value for the 
                       fixed-effect coefficients which were \emph{estimated} by the fitting function.}
    \item{\code{FixedCoef}}{A data frame containing the values of the fixed-effect coefficients which were \emph{set} 
                            to a fixed value.}
    \item{\code{Phi}}{A data frame containing the estimates, standard error, z values and P values for the overdispersion 
                      coefficients which were \emph{estimated} by the fitting function. Because the overdispersion 
                      coefficients are \eqn{> 0}, P values correspond to unilateral tests.}
    \item{\code{FixedPhi}}{A data frame containing the values of the overdispersion coefficients which were \emph{set} 
                            to a fixed value.}
    }
  }

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.glimML")}}
    \item{show}{\code{signature(object = "glimML")}}
    \item{summary}{\code{signature(object = "glimML")}}
    }
  }

\examples{
  data(orob2)
  fm1 <- betabin(cbind(y, n - y) ~ seed, ~ 1, data = orob2)
  # show method for objects of class "glimML"
  fm1
  # summary method for objects of class "glimML"
  res <- summary(fm1)
  res@Coef
  # show method for objects of class "summary.glimML"
  res
  }

\keyword{classes}
