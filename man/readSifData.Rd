\name{readSifData}
\alias{readSifData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reloaded System Matrix}
\description{
  The function reloaded a system matrix stored by \code{iradonIT}.
}
\usage{
readSifData(filename, DebugLevel = "Normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ (character)  \code{filename} gives the name of the ".sif"-file to read, including the pathname to the file. }
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="Normal"} for a standard level output. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"HardCore"} for no information at all. }
}
\details{
  If fast iterative reconstruction methods (\code{iradonIT} with parameters \code{UseFast = 1} and \code{KernelFileSave = 1}) are used, the system matrix will be saved with the name \code{KernelFileName}. See to help of \code{iradonIT} to get detail information about these parameters. The speed of the reconstruction can obviously be increased by using of an existing system matrix. \cr
  
  The system matrix will be stored as a binary ".sif"-file with the following structure: At first a header is defined, then the raw data follows. The form of the header is shown below. The raw data will be stored in a ".sif"-file in the following form: \cr
  At first \code{M} integers \eqn{N_m, m=1,\ldots,M} are stored, which specify the number of value in each row \eqn{m}. Then the \eqn{N_m} indicies (integers) follow alternatingly for the appropriate columns and \eqn{N_m} accordingly (float) values for each row. Thus, a sparse matrix is defined, i.e. a \eqn{(M,N)}-matrix without zero-elements. \cr

Header:    
\tabular{llll}{
  [1] \tab XSamples (integer) \tab [18] \tab mode (integer) \cr
 [2] \tab YSamples (integer) \tab [19] \tab UseFast (integer) \cr
 [3] \tab DeltaX (float) \tab [20] \tab ConstrainMin (float) \cr
 [4] \tab DeltaY (float) \tab [21] \tab ConstrainMax (float) \cr
 [5] \tab Xmin (float) \tab [22] \tab Alpha (float) \cr
 [6] \tab Ymin (float) \tab [23] \tab Beta (float) \cr
 [7] \tab ThetaSamples (integer) \tab [24] \tab IterationType (integer) \cr
 [8] \tab RhoSamples (integer) \tab [25] \tab KernelFileSave (integer) \cr
 [9] \tab DeltaRho (float) \tab [26] \tab SaveIterations (integer) \cr
 [10] \tab DeltaTheta (float) \tab [27] \tab RefFileName (200-byte character) \cr
 [11] \tab ThetaMin (float) \tab [28] \tab KernelFileName (200-byte character) \cr
 [12] \tab RhoMin (float) \tab [29] \tab InFileName (200-byte character) \cr
 [13] \tab LowestALevel (float) \tab [30] \tab OutFileName (200-byte character) \cr
 [14] \tab RadonKernel (integer) \tab [31] \tab SaveIterationsName (200-byte character) \cr
 [15] \tab OverSamp (integer) \tab [32] \tab M (integer) \cr
 [16] \tab Regularization (float) \tab [33] \tab N (integer) \cr
 [17] \tab Iterations (integer) \tab \tab \cr
}

\code{M = RhoSamples*ThetaSamples}, \code{N = XSamples*YSamples}, \code{InFileName = "InputData"} and \code{OutFileName = "OutputData"}. \cr
The other parameters are described in detail in the function \code{iradonIT}, but note that \code{RadonKernel} will be stored as an integer ("NN" = 1, "RNN" = 2, "RL" = 3, "P1" = 4 and "P2" =  5), \code{mode} will be stored as an integer ("ART" = 0, "EM" = 1, "CG" = 2) and \code{IterationType} will be stored as an integer ("random" = 0, "cyclic" = 1). 

}

\value{
  \item{MIndexNumber}{A vector of \code{M} integers \eqn{N_m, m=1,\ldots,M}, that specify the number of values in each row \eqn{m}.}
  \item{SparseIndex}{A list of \eqn{N_m} indicies (integers) for the appropriate columns. }
  \item{SparseSignal}{A list of \eqn{N_m} accordingly (float) values for each row.}
  \item{Header}{A list of all values defined above in Header.}
}

\references{ Toft, Peter, \emph{Ph.D. Thesis, The Radon Transform - Theory and Implementation}, Department of Mathematical Modelling Section for Digital Signal Processing, Technical University of Denmark, 1996. \cr

Schulz, Joern, \emph{Diploma Thesis: Analyse von PET Daten unter Einsatz adaptiver Glaettungsverfahren}, Humboldt-Universitaet zu Berlin, Institut fuer Mathematik, 2006.}

\author{ Joern Schulz, \email{jschulz78@web.de}. }

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
P <- phantom(n=101)
rP <- markPoisson(P, nSample=1800000)
irP <- iradonIT(rP$rData, 101, 101, KernelFileSave=1, 
              KernelFileName="SystemMatrix.sif")
rm(irP,P,rP)

# reading the sif-matrix:
SysMat <- readSifData("SystemMatrix.sif")
names(SysMat)
SysMat$Header
rm(SysMat)
}

}
\keyword{IO}% at least one, from doc/KEYWORDS
\keyword{file}% __ONLY ONE__ keyword per line
\concept{iradonIT}
\concept{radon}
%\concept{}
