\name{hough}
\alias{hough}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hough Transformation }
\description{
  The function implements a Hough transformation for an image.
}
\usage{
hough(oData, mode=1, XYSamples=nrow(oData), DeltaXY=1.0, 
      XYmin=-0.5*DeltaXY*(XYSamples-1), ThetaSamples=181, 
      RhoSamples=2*round(sqrt(sum((dim(oData))^2))/2)+1,
      ThetaMin=0, RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1),
      DeltaTheta=pi/ThetaSamples, DeltaRho=(2*abs(RhoMin)+1)/RhoSamples,       
      DebugLevel = "Normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oData}{  (matrix)  A matrix that contains the image (for the Hough transformation). }
  \item{mode}{  (integer) The optimization strategy for the Hough transformation. Default is \code{mode=1} for no optimization. Also implemented are \code{mode=2} for 'vector mapping', \code{mode=3} for 'avoid limit check and use vector mapping' and \code{mode=4} for 'avoid limit check and use matrix mapping'. }
  \item{XYSamples}{ (integer)  Specifies the number of samples on the x-axis (rows) and y-axis (columns) of \code{oData}. Defaults to \code{XYSamples=nrow(oData)}. }
  \item{DeltaXY}{ (double)  Specifies the sampling distance of both axes in the image. Defaults to \code{DeltaXY=1}. }
  \item{XYmin}{ (double)  Specifies the minimum sample position in the image on the first and second axis. If not given, the image is centered around the middle. Defaults to \code{XYmin=-0.5*DeltaXY*(XYSamples-1)}. }
  \item{ThetaSamples}{ (integer)  Specifies the number of samples in the angular parameter \eqn{\theta}{theta} in the Hough transformation image. It is sampled linearly from \eqn{0} to (approximately) \eqn{\pi}{pi} radians. Defaults to \code{ThetaSamples=181}. }
  \item{RhoSamples}{ Specifies the number of samples in the distance parameter \eqn{\rho}{rho} in the Hough transformation image. Defaults to \code{RhoSamples=2*round(sqrt(sum((dim(oData))^2))/2)+1}. }
  \item{ThetaMin}{ (double)  Specifies the minimum sample position in the Hough transformation image on the first axis. Defaults to \code{ThetaMin=0}. }
  \item{RhoMin}{ (double)  Specifies the minimum sample position in the Hough transformation image on the second axis. Defaults to \code{RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1)}. }
  \item{DeltaTheta}{ (double)  Angular sampling distance. Defaults to \code{DeltaTheta=pi/ThetaSamples}. }
  \item{DeltaRho}{ (double)  Specifies the sampling distance in \eqn{\rho}{rho}. The program will center the sampling points around \eqn{0}. Defaults to \code{DeltaRho=(2*abs(RhoMin)+1)/RhoSamples} }
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="Normal"} for a standard level output. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"HardCore"} for no information at all. }
}
\details{
  It is shown in the paper of P.Toft that the Hough transformation can be defined in a way that gives exactly the same discrete parameter domain as found with the nearest neighbour approximation of the discrete Radon transformation. The Hough transformation does not have the same property. It is also shown that the Hough transformation behaves very differently when changing the sampling intervals in the discrete parameter domain, compared to the discrete Radon transformation. Furthermore different optimization schemes for the Hough transformation are described there.
}
\value{
  \item{hData}{ A matrix, that contains the Hough transformation of \code{oData}. }
  \item{Header}{ A list of following values:
  \describe{ 
    \item{SignalDim}{ The dimension of the \code{hData}.}
    \item{XYmin}{ The minimum x- and y-position in \code{hData}.}
    \item{DeltaXY}{ Sampling distance on the x- and y-axis in \code{hData}.}
  } }
  \item{call}{Arguments of the call to \code{hough}. }  
}
\references{  Toft, Peter, \emph{Ph.D. Thesis, The Radon Transform - Theory and Implementation}, Department of Mathematical Modelling Section for Digital Signal Processing, Technical University of Denmark, 1996. \cr
\url{http://eivind.imm.dtu.dk/staff/ptoft/ptoft_papers.html} }

\author{ Peter Toft, Joern Schulz \email{jschulz78@web.de}. }

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{radon}}, \code{\link{markPoisson}} }

\examples{
P <- phantom()
hP <- hough(P)
viewData(list(P, hP$hData), list("Phantom", "Hough transformed phantom"))
rm(P,hP)
}
\keyword{math}
\keyword{smooth}
\concept{radon}
\concept{hough}
