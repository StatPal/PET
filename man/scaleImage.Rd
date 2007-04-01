\name{scaleImage}
\alias{scaleImage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scales an Image}
\description{
  Scales the minimum of an image to \eqn{0} and the mean or the maximum to \eqn{1}.
}
\usage{
scaleImage(image, mode = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{ The raw image data. }
  \item{mode}{ \code{mode} specifies if mean or maximum will be scaled to \eqn{1}. Defaults to \code{mode="mean"}. Alternative choice is \code{"max"}. }
}

\value{
  Returns the scaled data.
}

\author{ Joern Schulz, \email{jschulz78@web.de}. }

\examples{
(A <- matrix(-9:10,nrow=4))
scaleImage(A)
rm(A)
}
\keyword{math}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line
