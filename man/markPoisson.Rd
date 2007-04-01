\name{markPoisson}
\alias{markPoisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Marked Poisson Process}
\description{
  The function implements the Acceptance Rejection (AR) procedure to simulate a marked Poisson process with spatial varying intensity. Therewith it is possible to create noisy data from a phantom, as generated in Positron Emission Tomography (PET).
}
\usage{
markPoisson(DataInt, nSample = 200000, ThetaSamples = 181, 
            RhoSamples = 2*round(sqrt(sum((dim(DataInt))^2))/2)+1,
            RhoMin = -0.5*sqrt(2),
            vect.length = 100000, image = TRUE, DebugLevel = "Normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataInt}{ (matrix)  A two dimensional image in which the matrix elements comply with intensities. That means the intensity of the decay of positrons in a certain tissue. }
  \item{nSample}{ (integer)  \code{nSample} determines the number of accepted events that will be generated with AR method. Defaults to \code{nSample = 200000}. }
  \item{ThetaSamples}{ (integer)  Specifies the number of samples in the angular parameter \eqn{\theta}{theta} in the sinogram. The sinogram is sampled linearly from \code{0} to (approximately) \eqn{\pi}{pi} radians. Defaults to \code{ThetaSamples=181}. }
  \item{RhoSamples}{ (integer)  Specifies the number of samples in the distance parameter \eqn{\rho}{rho} in the sinogram. Defaults to \code{RhoSamples=2*round(sqrt(sum((dim(DataInt))^2))/2)+1}. }
  \item{RhoMin}{ (double)  Specifies the minimum sample position in the sinogram on the second axis. It is assumed that \code{DataInt} is scaled to \eqn{[-0.5,0.5]^2}. Defaults to \code{RhoMin=-0.5*sqrt(2)}. }
  \item{vect.length}{ (integer)  Determines a bound of number of generated accepted events in each iteration. That means, if \code{nSample > vect.length} then in each iteration \code{vect.length/(0.1217)} events will be generated and roughly \code{vect.length} accepted events. If you scale down \code{vect.length} the computation will be more time-consuming and if you increase the value of \code{vect.length}, more memory will be necessary. Defaults to \code{vect.length = 100000}. }
  \item{image}{ (logical)  If \code{image=FALSE}, only a generated sinogram will be returned, i.e. a noisy Radon transformation of the phantom will be returned. Otherwise, if \code{image=TRUE} an additional image will be returned, where each pixel of the image corresponds to a number of accepted events generated with the AR method. Defaults to \code{image=TRUE}. }
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Available are: The default \code{DebugLevel="Normal"} for standard level output or alternative \code{"Detail"} if it desirable to logged almost all output to screen or \code{"HardCore"} for no information at all. }
}
\details{
The function implements the Acceptance Rejection (AR) method to simulate a marked Poisson process with spatial varying intensity. The function was developed to simulate data of a PET scanner. Therefore, new random events are generated with the AR method from \code{DataInt}. An angle \eqn{\theta}{theta} will be generated to each event (uniform and iid. in \eqn{[0,\pi]}{[0,pi]}), because a PET scanner detects two photons, who will travel in (nearly) opposite directions. The line between the two detectors has the line parameters \eqn{(\rho,\theta)}{(rho,theta)}, whereas \eqn{\rho}{rho} is the shortest distance from the origin of the coordinate system to the line, and \eqn{\theta}{theta} an angle corresponding to the angular orientation of the line.  The only obtainable information from the two photons is the fact, that the photon emission took place somewhere along that line.
}
\value{
  \item{rData}{ A matrix, that contains the Radon transformed data. }
  \item{Data}{ Will be returned only in case of \code{image = TRUE}. A matrix, where each pixel of the image corresponds to a number of accepted events generated with the AR method. }
  \item{Header}{ A list of following values:
  \describe{ 
    \item{SignalDim}{ The dimension of the \code{rData}.}
    \item{XYmin}{ The minimum x- and y-position in \code{rData}.}
    \item{DeltaXY}{ Sampling distance on the x- and y-axis in \code{rData}.}
  } }
  \item{call}{Arguments of the call to \code{markPoisson}. }
}
\references{Schulz, Joern, \emph{Diploma Thesis: Analyse von PET Daten unter Einsatz adaptiver Glaettungsverfahren}, Humboldt-Universitaet zu Berlin, Institut fuer Mathematik, 2006. \cr
 
Gentle, J.E., \emph{Elements of Computational Statistics}, Springer-Verlag New York Berlin Heidelberg, 2002. }

\author{ Joern Schulz \email{jschulz78@web.de}.}

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{radon}}, \code{\link{iradon}}, \code{\link{iradonIT}} }

\examples{
\dontrun{
P <- phantom()
rP <- radon(P)
mP1 <- markPoisson(P)
mP2 <- markPoisson(P, nSample = 1000000)
viewData(list(P, mP1$Data, mP2$Data, rP$rData, mP1$rData, mP2$rData),
         list("Phantom", "nSample = 200000", "nSample = 1000000",
              "Radon Transfom of Phantom",  "nSample = 200000", 
              "nSample = 1000000"))
cat("Number of generated accepted events for mP2:",sum(mP2$Data),"\n")
rm(mP1,mP2,P,rP)
}
}

\keyword{math}
\keyword{smooth}
\concept{distribution}
\concept{radon}
\concept{positronen}
\concept{emission}
\concept{tomographie}
\concept{pet}
\concept{sinogram}
