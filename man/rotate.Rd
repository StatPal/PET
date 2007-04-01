\name{rotate}
\alias{rotate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rotates Data}
\description{
  The function implements the rotation of a vector, a matrix or a 3-dimensional array.
}
\usage{
rotate(A, grad = 90)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{Is a vector, a matrix or a 3-dimensional array.}
  \item{grad}{Can be set to 0, 90, 180, 270, -90, -180 or -270. Defaults to \code{grad = 90}.}
}
\value{
 Returns the rotated data.
}
\author{ Joern Schulz, \email{jschulz78@web.de}. }

\examples{
# Rotation of a vector
(x <- c(1:5))
rotate(x,270)
rm(x)

# Rotation of a matrix
(A <- matrix(1:15, nrow=3))
rotate(A, -90)
rm(A)

# Rotation of an array
(A <- array(1:20, c(2,5,2)))
rotate(A, 180)
rm(A)
}
\keyword{math}% at least one, from doc/KEYWORDS
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
