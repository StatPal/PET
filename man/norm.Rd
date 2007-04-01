\name{norm}
\alias{norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ L1 and L2 norm }
\description{
  Computes the L1 or L2 norm between two vectors or matrices.
}
\usage{
norm(orgImage, testImage, mode = "L2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orgImage}{ A vector or a matrix.}
  \item{testImage}{ A vector or a matrix of same size as \code{orgImage}. }
  \item{mode}{ Defaults to \code{mode="L2"} to compute a L2 norm between \code{orgImage} and \code{testImage}. Alternative choice is \code{"L1"}. }
}
\details{
  If \eqn{X} and \eqn{Y} are vectors of length \eqn{n} or \eqn{(n_1,n_2)}{(n1,n2)}-matrices with \eqn{n=n_1*n_2}{n=n1*n2} then the L1-norm and the L2-norm will be compute as follow:
  
  \deqn{L1 = \frac{1}{n} \sum_{i=1}^{n}\left|X_i-Y_i\right| }{%
        L1 = 1/n*(\sum_{i=1}^{n} |X_i-Y_i|) }
  \deqn{L2 = \left(\frac{1}{n} \sum_{i=1}^{n}\left(X_i-Y_i\right)^2\right)^{\frac{1}{2}} }{%
        L2 = (1/n*(\sum_{i=1}^{n} (X_i-Y_i)^2))^0.5 }  
}
\value{
  Returns the L1 or L2 norm.
}

\author{ Joern Schulz, \email{jschulz78@web.de}. }

\examples{
P <- phantom(n=101)
P.pois <- markPoisson(P)$Data
cat("The L1 is:",norm(P, P.pois, mode="L1"),"\n")
rm(P,P.pois)
}

\keyword{math}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line
