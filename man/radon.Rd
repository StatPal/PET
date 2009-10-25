\name{radon}
\alias{radon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Radon Transformation}
\description{
  This function will forward project (normal Radon transformation) a quadratic image into a sinogram.
}
\usage{

radon(oData, mode="NN", XYSamples=nrow(oData), 
      XYmin=-0.5*(nrow(oData)-1), DeltaXY=1, ThetaSamples=181, 
      RhoSamples=2*round(sqrt(sum((dim(oData))^2))/2)+1, ThetaMin=0,
      RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1),
      DeltaTheta=pi/ThetaSamples, DeltaRho=(2*abs(RhoMin)+1)/RhoSamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oData}{ (matrix)  A matrix that contains the image (for the Radon transformation). }
  \item{mode}{ (character)  The interpolation procedure. Currently supported functions are \code{mode="NN"} (Nearest Neighbour Interpolation), "LI" (Linear Interpolation) and "SINC" (Sinc Interpolation). Defaults to \code{mode="NN"}. }
  \item{XYSamples}{ (integer)  Specifies the number of samples on the \eqn{x}-axis (rows) and \eqn{y}-axis (columns) in \code{oData}. Defaults to \code{XYSamples=nrow(oData)}. }
  \item{XYmin}{ (double)  Specifies the minimum sample position in the image on the first and second axis. If not given, the image is centered around the middle. Defaults to \code{XYmin=-0.5*(nrow(oData)-1)}. }
  \item{DeltaXY}{ (double)  Specifies the sampling distance of both axes in the image. Defaults to \code{DeltaXY=1}. }
  \item{ThetaSamples}{ (integer)  Specifies the number of samples in the angular parameter \eqn{\theta}{theta} in the sinogram. The sinogram is sampled linearly from \eqn{0} to (approximately) \eqn{\pi}{pi} radians. Defaults to \code{ThetaSamples=181}. }
  \item{RhoSamples}{ (integer)  Specifies the number of samples in the distance parameter \eqn{\rho}{rho} in the sinogram. Should be an odd number. \cr
  Defaults to \code{RhoSamples = 2*round(sqrt(sum((dim(oData))^2))/2)+1}. }
  \item{ThetaMin}{ (double)  Specifies the minimum sample position in the sinogram on the first axis. Defaults to \code{ThetaMin=0}. }
  \item{RhoMin}{ (double)  Specifies the minimum sample position in the sinogram on the second axis. Defaults to \code{RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1)}. }
  \item{DeltaTheta}{ (double)  Angular sampling distance. Defaults to \code{DeltaTheta=pi/ThetaSamples}. }
  \item{DeltaRho}{ (double)  Specifies the sampling distance in \eqn{\rho}{rho}. The program will center the sampling points around \eqn{0}. Defaults to \code{DeltaRho=(2*abs(RhoMin)+1)/RhoSamples}. }
}

\details{
The function implements a forward projection of a quadratic image. The used projection is the Radon transformation, named after Johann Radon (1917). The normal form of a line is given by \deqn{\rho=x\cos(\theta)+y\sin(\theta).}{rho=x*cos(theta)+y*sin(theta).} There are different ways to define a Radon transformation of a two dimensional continuous function \eqn{f}. One possibility is to integrate values of \eqn{f} along lines, whereas the location of the line is determined by the parameter \eqn{\rho}{rho}, the shortest distance from the origin of the coordinate system to the line, and \eqn{\theta}{theta} an angle corresponding to the angular orientation of the line. If negative values of \eqn{\rho}{rho} are introduced the parameter domain is bounded by \eqn{0 \leq \theta < \pi}{0 <= theta < pi} and \eqn{-\rho_{max} \leq \rho \leq \rho_{max}}{RhoMin <= rho <= RhoMax, whereas RhoMin=-RhoMax}. Numerical computations use the discrete Radon transformation. Therefore all parameters are sampled linearly\cr
\tabular{l}{
\eqn{x=x_m=XMin+m*DeltaX, m=0,\ldots,XSamples-1},\cr
\eqn{y=y_n=YMin+n*DeltaY, n=0,\ldots,YSamples-1},\cr
\eqn{\rho=\rho_r=RhoMin+r*DeltaRho, r=0,\ldots,RhoSamples-1} and\cr
\eqn{\theta=\theta_t=ThetaMin+t*DeltaTheta, t=0,\ldots,ThetaSamples-1}.\cr
}
Given the digitally sampled image a fundamental problem arises. The discrete Radon transformation requires samples that are not found in the digital image, because linear sampling of all variables implies that sampling parameter in the Radon transformation in general never coincides with the sampling parameter of the digital image. That should be clear and it is necessary to approximate the \eqn{x} and \eqn{y} by using \code{"NN"} (nearest neighbour approximation), \code{"LI"} (linear interpolation) or \code{"SINC"} (sinc interpolation). \cr
Several things must be fulfilled to ensure a reasonable performance. Firstly sampling must be adequate in all parameters (see to references to get detail information). This will imply bounds on the sampling intervals. Secondly it is assumed that the fundamental function \eqn{f(x,y)} to be reconstructed have compact support, or more precisely is zero if \eqn{\sqrt{x^2+y^2} > |\rho_{max}|}{sqrt(x^2+y^2) > |rho_{max}|}. This demand will ensure that \eqn{Rf(\rho,\theta)=0}{Rf(rho,theta)=0}, if \eqn{|\rho|>|\rho_{max}|}{|rho|>|rho_{max}|}. If this cannot be fulfilled, numerical problems must be expected. \cr
Assuming that \eqn{f(x,y)} has compact support, then \eqn{(x,y)=(0,0)} should be placed to minimize \eqn{|\rho_{max}|}{|rho_{max}|}. This will reduce the size of the data array used for the discrete Radon transform. \cr

The forward projection methods that are implemented in this function are developed by P. Toft (1996) and implemented in R by J. Schulz (2006). For example see to references to get more detailed theoretical information about the Radon transformation.
}

\value{
  \item{rData}{ A matrix, that contains the Radon transformed data. }
  \item{Header}{ A list of following values:
  \describe{ 
    \item{SignalDim}{ The dimension of the \code{rData}.}
    \item{XYmin}{ The minimum x- and y-position in \code{rData}.}
    \item{DeltaXY}{ Sampling distance on the x- and y-axis in \code{rData}.}
  } }
  \item{call}{Arguments of the call to \code{radon}. }
}

\references{ Toft, Peter, \emph{Ph.D. Thesis, The Radon Transform - Theory and Implementation}, Department of Mathematical Modelling Section for Digital Signal Processing, Technical University of Denmark, 1996. \cr
\url{http://eivind.imm.dtu.dk/staff/ptoft/ptoft_papers.html}

Schulz, Joern, \emph{Diploma Thesis: Analyse von PET Daten unter Einsatz adaptiver Glaettungsverfahren}, Humboldt-Universitaet zu Berlin, Institut fuer Mathematik, 2006.

Radon, Johann, \emph{Ueber die Bestimmung von Funktionen durch ihre Integralwerte laengs gewisser Mannigfaltigkeiten}, Berichte ueber die Verhandlungen der koeniglich saechsischen Gesellschaft der Wissenschaften zu Leipzig, Mathematisch-Physische Klasse, Band 69, 1917.}

\author{ Joern Schulz \email{jschulz78@web.de}, Peter Toft.}

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{markPoisson}}, \code{\link{iradon}}, \code{\link{iradonIT}} }

\examples{
P <- phantom()
R1 <- radon(P)
R2 <- radon(P, DeltaRho=0.5)
R3 <- radon(P, RhoSamples=501)
R4 <- radon(P, RhoSamples=251, RhoMin=-125, DeltaRho=1)
R5 <- radon(P, RhoSamples=451, RhoMin=-225, DeltaRho=1)
viewData(list(P, R1$rData, R2$rData, R3$rData, R4$rData, R5$rData), 
         list("Phantom P with 256x256 voxels", "Default sinogram of P", 
              "DeltaRho=0.5", "RhoSamples=501", 
              "RhoSamples=251, RhoMin=-125, DeltaRho=1",
              "RhoSamples=451, RhoMin=-225, DeltaRho=1"))
rm(P,R1,R2,R3,R4,R5)
}

\keyword{math}
\keyword{smooth}
\concept{radon}
\concept{iradon}
\concept{positronen}
\concept{emission}
\concept{tomographie}
\concept{pet}
\concept{sinogram}
