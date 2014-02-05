\name{mcroastyx}
\alias{mcroastyx}

\title{Multiple Contrast Rotation Tests}
\description{Adaption of the roast function in package limma to perform a maxT-test in the multivariate linear model, controlling FWER for multiple contrasts.}
\usage{
mcroastyx(Y, X, K, nrot = 9999, adjusted = TRUE, moderated = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A matrix of response variables with individuals in rows and responses in columns}
  \item{X}{A design matrix with the same number of rows as Y}
  \item{K}{A contrast matrix, with different contrasts in rows, and the same number of columns as X}
  \item{nrot}{Number of Monte Carlo simulation steps}
  \item{adjusted}{logical, if TRUE single-step maxT FWER control is performed; if FALSE the comparison-wise error rate is controlled}
  \item{moderated}{logical, if TRUE, moderated statistics (Smyth 2004) are used.}
}
\details{
Most of the code was taken from function roast in package limma (Bioconductor).
}
\value{
A list with
\item{Y}{The matrix of response variables}
\item{X}{The design matrix}
\item{statistic}{A matrix of test statistics}
\item{pvalue}{A matrix of (adjusted) p-values}
\item{B}{A matrix with contrast parameter estimates}
\item{SD}{A matrix with estimated standard deviations}
\item{mrot}{A vector with simulated maxT statistics under H0}
}
\references{
Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments. Bioinformatics 26, 2176-2182. 

Langsrud, O (2005). Rotation tests. Statistics and Computing 15, 53-60.
}
\author{
Daniel Gerhard
}


\keyword{multivariate}
\keyword{htest}
