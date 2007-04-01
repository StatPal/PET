\name{phantom}
\alias{phantom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creation of a Phantom}
\description{
  This function creates two dimensional phantom data. There are different options, e.g. to generate a big ellipse that represents the head, and several smaller ellipses, that represent pathological areas to be located in the space of the bigger ellipse.
}
\usage{
phantom(n = 257, design = "A", addIm = "none", DebugLevel = "Normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ (integer)  Is the number of columns and rows in the generated phantom. It is assumed that the number of columns is equal to the number of rows. Defaults to \code{n=257}. }
  \item{design}{ (character)  \code{design} characterizes the phantom data. It is possible to define different ellipses, with different intensities. There are four default-designs, these are \code{design = "A"}, \code{"B"}, \code{"C"} and \code{"D"}. To define an own design of the phantom see below to details. Defaults to \code{design = "A"}. }
  \item{addIm}{ (character)  Adds an additional image to the phantom. There are six default-designs for \code{addIm}. These are \code{addIm="blurred1"}, \code{"blurred2"}, \code{"keen1"}, \code{"keen2"}, \code{"simple1"} and \code{"simple2"}. The defaults of \code{addIm} will be generated with the function \code{partEllipse}. If \code{addIm="none"}, no image will be added to the phantom. A further possibility is a matrix \code{A} (\code{addIm=A}) with the same size as the phantom (i.e. \code{dim(A)==c(n,n)}). Defaults to \code{addIm = "none"}.}
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="Normal"} for a standard level output. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"HardCore"} for no information at all. }
}
\details{
  \code{design}:
  
  To define an own design of the phantom (e.g. \code{design=P}) you have to note the following conditions:
  \enumerate{
    \item \code{P} is a \eqn{(n,5)} or \eqn{(n,6)} matrix, whereas \eqn{n > 0}.
    \item All elements \eqn{x} of \code{P} have to be between -1 and 1.
    \item Each row from \code{P} define a ellipse on \eqn{[-1,1]^2}, where
  }
  \tabular{rll}{
    [,1] \tab \eqn{A} \tab Is the additive intensity of the corresponding ellipse. \cr
    [,2] \tab \eqn{x_0} \tab Is the \eqn{x}-coordinate of the center of the ellipse. \cr
    [,3] \tab \eqn{y_0} \tab Is the \eqn{y}-coordinate of the center of the ellipse. \cr
    [,4] \tab \eqn{a} \tab Is the half length of the horizontal axis. \cr
    [,5] \tab \eqn{b} \tab Is the half length of the vertical axis. \cr
    [,6] \tab \eqn{\alpha} \tab Is the angle in degree between the \eqn{x}-axis of the ellipse \cr 
    \tab \tab and the \eqn{x}-axis of the grid. This parameter is optional, \cr
    \tab \tab but the number of columns have to be the same in all rows. \cr 
  }

}
\value{
 Returns a phantom image.
}

\author{ Joern Schulz \email{jschulz78@web.de}. }

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{partEllipse}}, \code{\link{markPoisson}} }

\examples{
P1 <- phantom()
P2 <- phantom(addIm="blurred1")
PhPa1 <- c(0.5, 0, 0, 0.4, 0.6)
PhPa2 <- matrix( c(0.6, -0.35, 0, 0.4, 0.6,
         0.3,  0.5,  0, 0.2, 0.35), nrow=2, byrow=TRUE )
P3 <- phantom(design=PhPa1)
P4 <- phantom(design=PhPa2)
viewData(list(P1, P2, P3, P4), list("Default Phantom",  
         "addIm='blurred1'", "First new design",
         "Second new design"))
rm(P1,P2,P3,P4,PhPa1,PhPa2)
}

\keyword{math}
\keyword{smooth}
\concept{iradonIT}
\concept{iradon}
\concept{radon}
\concept{pet}

