\name{iradon}
\alias{iradon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Direct Inverse Radon Transformation}
\description{
The function implements different direct reconstruction methods for a radon transformed two-dimensional image.
}

\usage{
iradon(rData, XSamples, YSamples, mode = "FB", 
       Xmin = -sqrt(0.5)*(ncol(rData)/XSamples)*0.5*(XSamples-1), 
       Ymin = -sqrt(0.5)*(ncol(rData)/YSamples)*0.5*(YSamples-1), 
       DeltaX = sqrt(0.5)*(ncol(rData)/XSamples), 
       DeltaY = sqrt(0.5)*(ncol(rData)/YSamples), 
       InterPol=1, FilterTyp="Hamming1", oData=NULL, 
       DebugLevel="Normal", iniFile=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rData}{ (matrix) A matrix that contains the sinogram image (for the reconstruction functions). }
  \item{XSamples}{ (integer) Number of samples on the x-axis (rows) in the reconstructed image. }
  \item{YSamples}{ (integer) Number of samples on the y-axis (columns) in the reconstructed image. }
  \item{mode}{ (character) The reconstruction method. Default is \code{mode="FB"} for Filtered Backprojection. Also implemented are \code{"BF"} for Filtering After Backprojection, \code{"CNF"} for Central Slice (FFT based with Nearest Neighbour approximation), \code{"CBF"} for Central Slice (FFT based with Bilinear Interpolation), \code{"CC"}  for using a variation of the Central Slice (by the use of nonlinear sampling of the Radon domain with the use of the Chirp-z algorithm, the interpolation can be reduced to a simple one  dimensional linear interpolation compared to normal CS), \code{"CNC"} for Central Slice (Chirp-z based with Nearest Neighbour approximation),  \code{"CBC"} for Central Slice (Chirp-z based with Bilinear Interpolation) and \code{"Test"} to run all available reconstruction  routines on the same image and to compare with a L1-Norm and L2-Norm. In case of \code{"Test"} the surrender value is \code{NULL}. }
  \item{Xmin}{ (double) The minimum x-position of the reconstructed image. Defaults to \code{Xmin=-sqrt(0.5)*(ncol(rData)/XSamples)*0.5*(XSamples-1)}. }
  \item{Ymin}{ (double) The minimum y-position of the reconstructed image. Defaults to \code{Ymin=-sqrt(0.5)*(ncol(rData)/YSamples)*0.5*(YSamples-1)}. }
  \item{DeltaX}{ (double) Sampling distance on the x-axis. \cr
  Defaults set to \code{DeltaX = sqrt(0.5)*(ncol(rData)/XSamples)}. }
  \item{DeltaY}{ (double) Sampling distance on the y-axis. \cr
  Defaults set to \code{DeltaY = sqrt(0.5)*(ncol(rData)/YSamples)}. }
  \item{InterPol}{ (integer) Interpolation level. Used by Filtered Backprojection. Defaults \code{InterPol=1}. }
  \item{FilterTyp}{ (character) Filter type is only used by Filtered Backprojection (\code{"FB"}). Defaults to \code{FilterTyp="Hamming1"} for generalized Hamming Filter with \eqn{\alpha=0.5}. The alternative choice is \code{"Hamming2"} for generalized Hamming Filter with \eqn{\alpha=0.54} or a \code{"Ramp"} to get very good results of data without noise. }
  \item{oData}{ (matrix) If \code{mode="Test"}, measures of misfit can be made between \code{rData} and \code{oData}. \code{oData} is the original image of the radon transformed data. Defaults to \code{oData=NULL}. }
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="Normal"} for a standard level output. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"HardCore"} for no information at all. }
  \item{iniFile}{ (character)  If \code{iniFile!=NULL}, \code{iniFile} has to be the name of an ini-file including a pathname to the file. In the case of a specified \code{iniFile} all parameters are read from the file. Note that in this case contingently setting parameters (except for \code{DebugLevel}) in R are ignored when calling \code{iradon}. The parameters which are not specified in \code{iniFile} are set to defaults. Defaults to \code{iniFile=NULL}. }
}

\details{
The function implements different direct reconstruction methods for a radon transformed two-dimensional image. The different methods will be specified with parameter \code{mode}. The reconstruction methods are developed by P. Toft (1996) and implemented in R by J. Schulz (2006). For detail theoretical information about the radon-transformation see references. \cr
Several things must be fulfilled to ensure a reasonable performance. Firstly sampling must be adequate in all parameters (see to references to get detail information). This will imply bounds on the sampling intervals. Secondly it is assumed that the fundamental function \eqn{f(x,y)} to be reconstructed have compact support, or more precisely is zero if \eqn{\sqrt{x^2+y^2} > |\rho_{max}|}{sqrt(x^2+y^2) > |rho_{max}|}. This demand will ensure that \eqn{Rf(\rho,\theta)=0}{Rf(rho,theta)=0}, if \eqn{|\rho|>|\rho_{max}|}{|rho|>|rho_{max}|}. If this cannot be fulfilled, numerical problems must be expected. \cr
Assuming that \eqn{f(x,y)} has compact support, then \eqn{(x,y)=(0,0)} should be placed to minimize \eqn{|\rho_{max}|}{|rho_{max}|}. This will reduce the size of the data array used for the discrete Radon transform. \cr

Another variation to determine the parameter is the offer with \code{iniFile}. \code{iniFile} has to be a name of a file (e.g. \code{iniFile="/home/work/sino.ini"}) with the following structure. Each line begins with the name of parameter, then an equal sign follows and the value of the parameter. The first line must contain the parameter \code{mode}. Characters are not written in \code{""}, e.g. \code{RadonKernel=NN} and not \code{RadonKernel="NN"}. Furthermore note that in an ini-file \code{rData} and \code{oData} both have to be of type character, videlicet the name of the corresponding file. Supported file formats are ".txt", ".dat", ".fif", ".pet", ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif, ".jpg, ".jpeg. See \code{?readData} to get detailed information about supported formats. \cr
If a file of the type ".fif", ".pet" or ".dat" is specified for \code{oData} (see \code{?readData} or \code{?writeData} in R to get more information about these formats), the parameters \code{XSamples}, \code{YSamples}, \code{XMin}, \code{YMin}, \code{DeltaX} and \code{DeltaY} will be read from the file-header, but only in the case of unspecified corresponding parameters. \cr
Note that in the case of an ini-file contingently setting parameters (except for \code{DebugLevel}) in R are ignored with the calling of \code{iradon}. Parameters that are not specified in the \code{iniFile} are set to default.
}
\value{
  \item{irData}{ A matrix, that contains the reconstructed image (matrix) of \code{rData}. }
  \item{Header}{ A list of following values:
  \describe{ 
    \item{SignalDim}{ The dimension of the \code{irData}.}
    \item{XYmin}{ The minimum x- and y-position in \code{irData}.}
    \item{DeltaXY}{ Sampling distance on the x- and y-axis in \code{irData}.}
  } }
  \item{call}{Arguments of the call to \code{iradon}. }
}
\references{ Toft, Peter, \emph{Ph.D. Thesis, The Radon Transform - Theory and Implementation}, Department of Mathematical Modelling Section for Digital Signal Processing, Technical University of Denmark, 1996. \cr
\url{http://eivind.imm.dtu.dk/staff/ptoft/ptoft_papers.html}

Schulz, Joern, \emph{Diploma Thesis: Analyse von PET Daten unter Einsatz adaptiver Glaettungsverfahren}, Humboldt-Universitaet zu Berlin, Institut fuer Mathematik, 2006.}

\author{ Joern Schulz \email{jschulz78@web.de}, Peter Toft. }

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{radon}}, \code{\link{iradonIT}}}

\examples{
#    
# Compare the results of different direct reconstruction methods
#
\dontrun{
P <- phantom(design="B")
rP <- markPoisson(P, nSample=3000000 )
irP1 <- iradon(rP$rData , nrow(P), ncol(P))
irP2 <- iradon(rP$rData , nrow(P), ncol(P), 
        mode="BF", DebugLevel="HardCore")
irP3 <- iradon(rP$rData , nrow(P), ncol(P), 
        mode="CBF", DebugLevel="HardCore")
viewData(list(rP$rData, irP1$irData, irP2$irData, irP3$irData),
         list("Generated PET Data", "Reconstruction: mode='FB'", 
         "Reconstruction: mode='BF'", "Reconstruction: mode='CBF'"))
rm(irP1,irP2,irP3,P,rP)
}

#    
# Compare the results of different values for RhoSamples in 'markPoisson'
#
\dontrun{
P <- phantom()
rP1 <- markPoisson(P, nSample=1000000, RhoSamples=101, image=FALSE)
rP2 <- markPoisson(P, nSample=1000000, RhoSamples=256, image=FALSE)
rP3 <- markPoisson(P, nSample=1000000, RhoSamples=501, image=FALSE)
rP4 <- markPoisson(P, nSample=1000000, RhoSamples=801, image=FALSE)
irP1 <- iradon(rP1$rData, 257, 257)
irP2 <- iradon(rP2$rData, 257, 257, DebugLevel="HardCore")
irP3 <- iradon(rP3$rData, 257, 257, DebugLevel="HardCore")
irP4 <- iradon(rP4$rData, 257, 257, DebugLevel="HardCore")
viewData(list(irP1$irData, irP2$irData, irP3$irData, irP4$irData,),
         list("RhoSamples=101",  "RhoSamples=256", "RhoSamples=501",
         "RhoSamples=801"))
rm(P,rP1,rP2,rP3,rP4,irP1,irP2,irP3,irP4)
}

# 
# mode="Test"
#
P <- phantom()
R <- radon(P)
iradon(R$rData, XSamples=257, YSamples=257, mode="Test", oData=P)
rm(P,R)
}
\keyword{math}
\keyword{smooth}
\concept{iradon}
\concept{radon}
\concept{positronen}
\concept{emission}
\concept{tomographie}
\concept{pet}
\concept{sinogram}
