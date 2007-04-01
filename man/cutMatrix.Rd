\name{cutMatrix}
\alias{cutMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Cut Matrix }
\description{
  Scale the minimum of a matrix to 0 and cut off the first and the last rows or columns where the added up values are less than \eqn{1e-10}.
}
\usage{
cutMatrix(A, mode = "col")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ Describe the matrix. }
  \item{mode}{ The default is set to \code{mode="col"} to cut off the first and the last columns of \code{A} where the added up columns values are less than \eqn{1e-10}. Alternative choice is \code{mode="row"}. }
}

\value{
  \item{A}{ The cut off matrix. }
  \item{dimOrg}{ The dimension of the original image. }
  \item{pattern}{ A pattern, contains the first and last column or row of \code{A} where the added up values are greater than \eqn{1e-10}. }
}

\author{ Joern Schulz \email{jschulz78@web.de}. }

\examples{
A <- matrix(c(rep(0,6),0:2,0,3,4,rep(0,4),5,6,rep(0,9)),nrow=3)
A
cutMatrix(A)
cutMatrix(A,mode="row")
rm(A)
}
\keyword{math}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line
