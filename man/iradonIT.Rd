\name{iradonIT}
\alias{iradonIT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inverse Radon Transformation }
\description{
The function implements three of the major iterative reconstruction techniques: ART (Algebraic Reconstruction Technique), EM (Likelihood Reconstruction using Expectation Maximization) and LSCG (Least Squares Conjugate Method). 
}
\usage{
iradonIT(rData, XSamples, YSamples, StartImage = "None", 
         mode = "EM", UseFast = 1, RadonKernel = "NN", 
         Iterations = 20, IterationsType = "random", 
         SaveIterations = 0, SaveIterationsName = "", 
         LowestALevel = 0, ConstrainMin = -1, 
         ConstrainMax = -1, Alpha = 1, Beta = 1, 
         Regularization = 0, KernelFileSave = 0, 
         KernelFileName = "", RefFileName = "None", 
         ThetaSamples = nrow(rData), RhoSamples = ncol(rData), 
         ThetaMin = 0, RhoMin = -0.5*((2*round(sqrt(XSamples^2+
         YSamples^2)/2)+1)-1), DeltaTheta = pi/ThetaSamples, 
         DeltaRho = (2*abs(RhoMin)+1)/RhoSamples, 
         Xmin = -0.5*(XSamples-1), Ymin = -0.5*(YSamples-1), 
         DeltaX = 1, DeltaY = 1, OverSamp = 0, 
         DebugLevel = "Normal", iniFile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rData}{ (matrix)  A matrix that contains the sinogram. The rows of the image correspond to the sampled angles \eqn{\theta} and the columns correspond to the sampled distance \eqn{\rho}. These both parameters determine the lines. }
  \item{XSamples}{ (integer)  \code{XSamples} is the number of samples on the x-axis (rows) in the reconstructed image. }
  \item{YSamples}{ (integer)  \code{YSamples} is the number of samples on the y-axis (columns) in the reconstructed image. }
  \item{StartImage}{ (matrix)  If provided, then the initial guess on the reconstructed image is taken from this image, otherwise a constant solution will be assumed from the start. Default is set to \code{StartImage="None"}. Note, that the dimension of the \code{StartImage} has to be \code{dim(StartImage)==c(XSamples,YSamples)}. }
  \item{mode}{ (character)  The iterative reconstruction method. Default is \code{mode="EM"} for Likelihood Reconstruction using Expectation Maximization. Additional implementations are \code{"ART"} (Algebraic Reconstruction Technique) and \code{"CG"} (Least Squares Conjugate Method).}
  \item{UseFast}{ (logical)  If \code{UseFast=1} fast reconstruction is used, where the system matrix is stored (in case of \code{KernelFileSave=1}) using sparse techniques, otherwise a slower but more efficient memory reconstruction is used. Default is \code{UseFast=1}. }
  \item{RadonKernel}{ (character)  The type of kernel is used to model the system matrix. Default is \code{RadonKernel="NN"} for two-level Nearest Neighbour approximation (Memory consuming). Additional implementations are  \code{"RNN"} for ray driven Nearest Neighbour approximation discrete Radon transformation based (very fast with small system matrix), \code{"RL"} for ray driven Linear Interpolation discrete Radon transformation based (fast with small system matrix), \code{"P1"} for method based on Radon transformation of square with pre-guidance  (slow but good) and \code{"P2"} for method based on Radon transformation of square with no pre-guidance (slower but better). }
  \item{Iterations}{ (integer)  Using \code{mode="EM"} or \code{"CG"}, it is the number of iterations before the  iteration ends. Using \code{mode="ART"}, it is the number of full iterations, i.e., divided by the number of rows in the system matrix.  Defaults to \code{Iterations=20}. }
  \item{IterationsType}{ (character)  Using \code{mode="ART"}, two iterations types are possible:  \cr
  \code{"cyclic"} selection of the row index or \code{"random"} selection. Defaults to \code{IterationsType = "random"}. }
  \item{SaveIterations}{ (integer)  If it is set to \code{SaveIterations=1}, the current solution will be saved after each iteration under the name SaveIterationsName + CurrentIteration (also possible 2,3 ...). Defaults to \code{SaveIterations=0}. \cr
NOTE: When \code{Iterations/itINI.SaveIterations>300} then too many new files are requested. In this case default is used. }
  \item{SaveIterationsName}{ (character)  In case of \code{SaveIteration} is non-zero,  the parameter determines the name of the currently saved iteration. It contains the path to the file and the filename. The path contains the path relatively to the working-directory of your R-session or it contains the full path of the file. Defaults to  \code{SaveIterationsName=""}. }
  \item{LowestALevel}{ (double)  If fast reconstruction is used, the matrix elements will be truncated to this level relatively to the sampling distance of \eqn{x}. Defaults to  \code{LowestALevel=0.0}. }
  \item{ConstrainMin}{ (double)  After each iteration the solution in each sample will be forced above this limit. Defaults to \code{ConstrainMin=0}. }
  \item{ConstrainMax}{ (double)  After each iteration the solution in each sample will be forced below this limit. Defaults to \code{ConstrainMax=-1}. \cr
NOTE: For both limit it is assumed that negative limits imply that the feature is not used. }
  \item{Alpha}{ (double)  If \code{mode="ART"}, it will be the initial update weight. Defaults to \code{Alpha=1.0}. }
  \item{Beta}{ (double)  If \code{mode="ART"}, it will be the multiplicative change to the weight factor, which should be less than one. Default is set to \code{Beta=1}. }
  \item{Regularization}{ (double)  If the Regularization is not zero and positive and using fast reconstruction, rows will be appended to the system matrix with a simple Laplace operator. A weight factor should also be incorporated. If negative, identity matrix will be used. Defaults to \code{Regularization=0}. }
  \item{KernelFileSave}{ (logical)  If \code{KernelFileSave=1} and fast reconstruction is used, the system matrix will be saved under the name \code{KernelFileName}. Defaults to \code{KernelFileSave=0}. }
  \item{KernelFileName}{ (character)  If fast reconstruction is used, the system matrix will be saved and restored with this sif-name. If the system matrix is incompatible to the sampling parameters, a new one will be generated. KernelFileName contains the path to the filename and the filename. The path contains the path relatively to the working-directory of your R-session or it contains the full path to the file.  Defaults to \code{KernelFileName=""}. }
  \item{RefFileName}{ (character)  Error measures can be made between the image which is contained in \code{RefFileName} and the reconstructed image from \code{rData}. The file has to be either a fif-file or a pet-file (type \code{?writeData} or \code{?readData} in R to get more information about these formats). If \code{RefFileName} exists, the measures will be logged in a file with the name of \code{RefFileName} and extension 'dif'. Defaults to \code{RefFileName="None"}. }
  \item{ThetaSamples}{ (integer)  Number of angular samples \eqn{\theta} in the sinogram.  Defaults to \code{ThetaSamples=nrow(rData)}. }
  \item{RhoSamples}{ (integer)  Number of samples in the sinogram \eqn{\rho} in the \eqn{\rho}-direction. Defaults to \code{RhoSamples=ncol(rData)}. }
  \item{ThetaMin}{ (double)  Start of the angular sampling (should be set to zero).  Defaults to \code{ThetaMin=0}. }
  \item{RhoMin}{ (double)  Start of sampling positions in \eqn{\rho}{rho}.  Defaults to \code{RhoMin=-0.5*} \code{((2*round(sqrt(XSamples^2+YSamples^2)/2)+1)-1)}. }
  \item{DeltaTheta}{ (double)  Angular sampling distance. Should be set to the default \code{DeltaTheta = pi/ThetaSamples}. }
  \item{DeltaRho}{ (double)  Sampling distance in \eqn{\rho}{rho}. Defaults to \code{DeltaRho=(2*abs(RhoMin)+1)/RhoSamples}. }
  \item{Xmin}{ (double)  The minimum \eqn{x}-position of the reconstructed image. Defaults to \code{Xmin=-0.5*(XSamples-1)}. }
  \item{Ymin}{ (double)  The minimum \eqn{y}-position of the reconstructed image. Defaults to \code{Ymin=-0.5*(YSamples-1)}. }
  \item{DeltaX}{ (double)  Sampling distance on the \eqn{x}-axis. Defaults to \code{DeltaX=1}. }
  \item{DeltaY}{ (double)  Sampling distance on the \eqn{y}-axis. Defaults to \code{DeltaY=1}. }
  \item{OverSamp}{ (integer)  Uses oversampling - squared number of samples are used. Defaults to \code{OverSamp=0}. }
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="Normal"} for a standard level output. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"HardCore"} for no information at all. }
  \item{iniFile}{ (character)  If \code{iniFile!=NULL}, \code{iniFile} has to be the name of an ini-file including a pathname to the file. In the case of a specified \code{iniFile} all parameters are read from the file. Note that in this case contingently setting parameters (except for \code{DebugLevel}) in R are ignored when calling \code{iradonIT}. The parameters which are not specified in \code{iniFile} are set to defaults. Default of \code{iniFile} is \code{iniFile=NULL}. }
}
\details{
The function implements three of the major iterative reconstruction techniques: ART (Algebraic Reconstruction Technique), EM (Likelihood Reconstruction using Expectation Maximization) and LSCG (Least Squares Conjugate Method). The reconstruction methods were developed by P. Toft (1996) and implemented in R by J. Schulz (2006). For detailed theoretical information about the Radon-transformation see references. \cr

There are a lot of different parameters which influence the reconstruction method. In common cases, it should be enough to determine \code{rData}, \code{XSamples} and \code{YSamples}. For the default parameters, it is assumed that the sinogram is standardised to a coordinate system with \eqn{\theta}{theta} from \eqn{0} to \eqn{\pi} and \eqn{\rho}{rho} from \eqn{-0.5*((2*round(sqrt(M^2+N^2)/2)+1)-1)} to \eqn{0.5*((2*round(sqrt(M^2+N^2)/2)+1)-1)}. It is also assumed that the reconstructed image is standardised to a coordinate system with \eqn{x} from \eqn{-0.5*(M-1)} to \eqn{0.5*(M-1)} and \eqn{y} from \eqn{-0.5*(N-1)} to \eqn{0.5*(N-1)}. \eqn{R}, \eqn{M} and \eqn{N} are the number of sampling parameters in \eqn{\rho}, \eqn{x} and \eqn{y}. If that is not true, you will have to adapt the parameters. \cr
Several things must be fulfilled to ensure a reasonable performance. Firstly sampling must be adequate in all parameters (see to references to get detail information). This will imply bounds on the sampling intervals. Secondly it is assumed that the fundamental function \eqn{f(x,y)} to be reconstructed have compact support, or more precisely is zero if \eqn{\sqrt{x^2+y^2} > |\rho_{max}|}{sqrt(x^2+y^2) > |rho_{max}|}. This demand will ensure that \eqn{Rf(\rho,\theta)=0}{Rf(rho,theta)=0}, if \eqn{|\rho|>|\rho_{max}|}{|rho|>|rho_{max}|}. If this cannot be fulfilled, numerical problems must be expected. \cr
Assuming that \eqn{f(x,y)} has compact support, then \eqn{(x,y)=(0,0)} should be placed to minimize \eqn{|\rho_{max}|}{|rho_{max}|}. This will reduce the size of the data array used for the discrete Radon transform. \cr


The numerical complexity of the procedure is mainly determined by \code{Iterations}, the size of \code{rData}, \code{XSamples}, \code{YSamples} and by the use of fast reconstruction method. \cr
Fast reconstruction is used in case of \code{UseFast=1}. It is possible to store the computed system matrix (\code{KernelFileSave=1}) under the name \code{KernelFileName}. The speed of the reconstruction can obviously be increased by using of an existing system matrix. NOTE: This is only possible, if the sampling parameters are compatible with the system matrix. That means, that the fifteen parameters \code{XSamples}, \code{YSamples}, \code{DeltaX}, \code{DeltaY}, \code{Xmin}, \code{Ymin}, \code{ThetaSamples}, \code{RhoSamples}, \code{DeltaRho}, \code{DeltaTheta}, \code{ThetaMin}, \code{RhoMin}, \code{LowestALevel}, \code{RadonKernel} and \code{OverSamp} are equal to the parameter from which the system matrix is generated. If the system matrix is incompatible to the sampling parameters, a new one will be generated. \cr

Another variation to determine the parameter is the offer with \code{iniFile}. \code{iniFile} has to be the name of a file (e.g. \code{iniFile="/home/work/sino.ini"}) with the following structure. Each line contains at first the name of a parameter, then an equal sign follows and the value of the parameter. The first line must contain the parameter \code{mode}. Characters are not written in \code{""}, e.g. \code{RadonKernel=NN} and not \code{RadonKernel="NN"}. Furthermore note that in an ini-file \code{rData} and \code{StartImage} are to be of type character, videlicet the name of the corresponding file. Supported file formats are ".txt", ".dat", ".fif", ".pet", ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg", ".jpeg". See \code{?readData} to get detailed information about supported formats. \cr
If a file of the type ".fif", ".pet" or ".dat" is specified for \code{rData} (see \code{?writeData} or \code{?readData} in R to get more information about these formats), the parameters \code{ThetaSamples}, \code{RhoSamples}, \code{ThetaMin}, \code{RhoMin}, \code{DeltaTheta} and \code{DeltaRho} will be read from the file-header of \code{rData}, but only in case of unspecified corresponding parameters. Just as, if a file of type ".fif", ".pet" or ".dat" is specified for \code{StartImage}, then the parameters \code{XSamples}, \code{YSamples}, \code{XMin}, \code{YMin}, \code{DeltaX} and \code{DeltaY} will be read from the file-header. \cr
Note that in case of an ini-file, contingently setting parameters (except for \code{DebugLevel}) in R are ignored when calling \code{iradonIT}. Parameters that are not specified in the \code{iniFile} are set to defaults.
}

\value{
  \item{irData}{ A matrix, that contains the reconstructed image (matrix) of \code{rData}. }
  \item{Header}{ A list of following values:
  \describe{ 
    \item{SignalDim}{ The dimension of the \code{irData}.}
    \item{XYmin}{ The minimum x- and y-position in \code{irData}.}
    \item{DeltaXY}{ Sampling distance on the x- and y-axis in \code{irData}.}
  } }
  \item{call}{Arguments of the call to \code{iradonIT}. }
}

\references{ Toft, Peter, \emph{Ph.D. Thesis, The Radon Transform - Theory and Implementation}, Department of Mathematical Modelling Section for Digital Signal Processing, Technical University of Denmark, 1996. \cr
\url{http://eivind.imm.dtu.dk/staff/ptoft/ptoft_papers.html}

Schulz, Joern, \emph{Diploma Thesis: Analyse von PET Daten unter Einsatz adaptiver Glaettungsverfahren}, Humboldt-Universitaet zu Berlin, Institut fuer Mathematik, 2006.}

\author{ Joern Schulz \email{jschulz78@web.de}, Peter Toft.}

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{radon}}, \code{\link{iradon}}}

\examples{
#    
# Compare the results of iterative reconstruction method "EM" and 
# direct reconstruction method "FB"
#
\dontrun{
P <- phantom(design="B")
rP <- markPoisson(P, nSample=1600000 )
irP1 <- iradon(rP$rData , nrow(P), ncol(P))
irP2 <- iradonIT(rP$rData , nrow(P), ncol(P))
viewData(list(P, rP$rData, irP1$irData, irP2$irData),
         list("Generated unnoisy Phantom", "Generated PET Data",
         "Direct rec.: mode='FB'", "Iterative rec.: mode='EM'"))
rm(irP1,irP2,P,rP)
}

#    
# Compare the results from the iterative reconstruction methods "EM"
# "CG" and "ART"
#
\dontrun{
P <- phantom(n=151, addIm="blurred1")
rP <- markPoisson(P, nSample=2000000, RhoSamples=401)
irP1 <- iradonIT(rP$rData , nrow(P), ncol(P), Iterations=20, 
        ConstrainMin=0, ConstrainMax=10)
irP2 <- iradonIT(rP$rData , nrow(P), ncol(P), mode="CG", 
        Iterations=4, ConstrainMin=0, ConstrainMax=10)
irP3 <- iradonIT(rP$rData , nrow(P), ncol(P), mode="ART", 
        Iterations=5, ConstrainMin=0, ConstrainMax=10, 
        Alpha=0.1, Beta=0.5)
viewData(list(P,irP1$irData,irP2$irData,irP3$irData),
         list("Generated unnoisy Phantom",
         "mode='EM', Iterations=20","mode='CG', Iterations=4", 
         "mode='ART', Iter.=20, Alpha=0.1, Beta=0.5"))
rm(irP1,irP2,irP3,P,rP)
}

}
\keyword{math}
\keyword{smooth}
\concept{iradonIT}
\concept{iradon}
\concept{radon}
\concept{positronen}
\concept{emission} 
\concept{tomographie}
\concept{pet}
\concept{sinogram}
