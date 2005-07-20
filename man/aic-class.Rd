\name{aic-class}
\docType{class}
\alias{aic-class}

\title{Representation of Objects of Formal Class "aic"}

\description{Representation of the output of function \code{AIC}.}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("aic", ...)} or, more commonly,
  via the function \code{AIC}.}}

\section{Slots}{
  \describe{
    \item{\code{istats}}{A data frame with 3 columns describing the models indicated by the row names:
    
      \bold{df} the number of parameters in the model,
      
      \bold{AIC} the Akaike Information Criterion for the model (see \code{\link[aod:AIC]{AIC}}),
      
      \bold{AICc} the small-sample corrected Akaike Information Criterion for the model 
        (see \code{\link[aod:AIC]{AIC}}).
        
      }
    }
  }

\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "aic")}}
    \item{show}{\code{signature(object = "aic")}}
    }
  }

\keyword{classes}
