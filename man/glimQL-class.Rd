\name{glimQL-class}
\docType{class}
\alias{glimQL-class}
\alias{show,glimQL-method}
\title{Representation of Models of Formal Class "glimQL"}
\description{Representation of models of formal class "glimQL" fitted by quasi-likelihood method.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("glimQL", ...)} or, 
  more commonly, via the \code{quasibin} or \code{quasipois} functions.}
\section{Slots}{
  \describe{
    \item{\code{CALL}}{The call of the function.}
    \item{\code{fm}}{A fitted model of class \dQuote{glm}.}
    \item{\code{phi}}{The overdispersion parameter.}
    }
  }
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "glimQL")}: Main results of \dQuote{glimQL} models.}
    }
  }
\keyword{classes}
