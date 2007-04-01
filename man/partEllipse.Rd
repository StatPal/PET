\name{partEllipse}
\alias{partEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generates a Fragment of an Ellipse}
\description{
  The function offers the possibility to generate and to manipulate an ellipse. The data will be generated in \eqn{[-1,1]^2}.
}
\usage{
partEllipse(mod = "hcirc1", x = 0, y = 0, intensity = 0.1, 
            n = 257, re1 = 0.2, re2 = 0.2, ring.wide = 0.05, 
            no.xax1 = -1, no.xax2 = 1, no.yax1 = -1, no.yax2 = 1, 
            in.r = 0.01, DebugLevel = "HardCore")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{ Defines the type of computation technique. Default is \code{mod = "hcirc1"}, i.e. the ellipse is generated only by the parameters from \code{mod} to \code{no.yax2}. Also implemented are \code{mod = "hcirc2"} and \code{mod = "hcirc3"} where additionally in each point a circle with radius \code{in.r} will be generated in the area with an intensity. The effect is that the border of the ellipse is rounded. If \code{mod = "hcirc2"} all values in an accepted area will be set to \code{intensity}, but if \code{mod = "hcirc3"} then values in an accepted area will be cumulated. The result is more blurred. Note that in both cases the costs for computation increase obviously, all the more value exist with \code{intensity!=0}. }
  \item{x}{ Is the \eqn{x}-coordinate of the center of the ellipse. Defaults to \code{x = 0}. }
  \item{y}{ Is the \eqn{y}-coordinate of the center of the ellipse. Defaults to \code{y = 0}. }
  \item{intensity}{ The intensity of the values, those are in the defined area. Defaults to \code{intensity = 0.1}. }
  \item{n}{ (integer)  Is the number of columns and rows in the generated image. It is assumed that the number of columns are equal to the number of rows. Defaults to \code{n = 257}. }
  \item{re1}{ Is the half length of the horizontal axis. Defaults to \code{re1 = 0.2}. }
  \item{re2}{ Is the half length of the vertical axis. Defaults to \code{re2 = 0.2}. }
  \item{ring.wide}{ Defines the wide of the ellipse. This area is set to \code{intensity}. Note \eqn{0 < \code{ring.wide} <= 1}. Defaults to \code{ring.wide = 0.05}. }
  \item{no.xax1}{ The area of \eqn{x}-axis between \eqn{-1} and \code{no.xax1} is set to \eqn{0}. Defaults to \code{no.xax1 = -1}. }
  \item{no.xax2}{ The area of \eqn{x}-axis between \eqn{1} and \code{no.xax2} is set to \eqn{0}. Defaults to \code{no.xax2 = 1}. }
  \item{no.yax1}{ The area of \eqn{y}-axis between \eqn{-1} and \code{no.yax1} is set to \eqn{0}. Defaults to \code{no.yax1 = -1}. }
  \item{no.yax2}{ The area of \eqn{y}-axis between \eqn{1} and \code{no.yax2} is set to \eqn{0}. Defaults to \code{no.yax2 = 1}. }
  \item{in.r}{ If \code{mod = "hcirc2"} or \code{mod = "hcirc3"} \code{in.r} will define the radius of circles that will be generated in each point of the areas with \eqn{\code{itensity}>0}. If \code{mod="hcirc1"} \code{in.r} will be ignored. Defaults to \code{in.r = 0.01}. }
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="HardCore"} for no information at all. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"Normal"} for a standard level output. }
}

\value{
  Returns a matrix.
}

\author{ Joern Schulz \email{jschulz78@web.de}. }

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{phantom}} }

\examples{
\dontrun{
ellip1 <- partEllipse()
ellip2 <- partEllipse(x=0.1, y=-0.1, ring.wide=0.4)
ellip3 <- partEllipse(mod="hcirc2", y=0.2, intensity=0.4, 
                  re1=0.5, no.xax1=0, DebugLevel="Normal")
ellip4 <- partEllipse(mod="hcirc3", y=0.2, intensity=0.4, 
                  re1=0.5, no.xax2=0, DebugLevel="Normal")
viewData(list(ellip1, ellip2, ellip3, ellip4,),
         list("Part of an ellipse 1", "Part of an ellipse 2",
              "Part of an ellipse 3", "Part of an ellipse 4"))
rm(ellip1,ellip2,ellip3,ellip4)
}
}

\keyword{math}
\keyword{smooth}
\concept{iradonIT}
\concept{iradon}
\concept{radon}
\concept{pet}
