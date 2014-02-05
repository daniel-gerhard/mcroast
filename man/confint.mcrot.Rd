\name{confint.mcrot}
\alias{confint.mcrot}


\title{Simultaneous Confidence Intervals Inverting the Rotation Test}
\description{
Computes simultaneous confidence intervals with a quantile based on the simulated maxT-statistics.
}
\usage{
\method{confint}{mcrot}(object, parm, level=0.95, ...)
}

\arguments{
  \item{object}{an object of class mcrot}
  \item{parm}{ignored here}
  \item{level}{the confidence level required}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
A list with
\item{B}{A matrix with contrast parameter estimates}
\item{lower}{A matrix with estimated lower confidence bounds}
\item{upper}{A matrix with estimated upper confidence bounds}
}

\author{
Daniel Gerhard
}


\keyword{htest}
