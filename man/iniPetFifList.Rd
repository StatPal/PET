\name{iniPetFifList}
\alias{iniPetFifList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initializes a Header for PET or FIF File}
\description{
  The function initializes a list with values that corresponds to a header of a ".pet" or a ".fif" file. The file structures of these formats are explicitly described in the functions \code{readData} and \code{writeData} of this package.
}
\usage{
iniPetFifList(listType = "pet", data = NULL, imType = "normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listType}{ Can be chosen between \code{"pet"} and \code{"fif"}. Default is \code{listType = "pet"} to  generate a default header for ".pet"-files. Choose \code{listType = "fif"} to  generate a default header for ".fif"-files. }
  \item{data}{ Using \code{data}, it is possible to consign the matrix corresponding to the list and to generate the parameters \code{SignalDim}, \code{XYmin}, \code{DeltaXY} and \code{SignalMinMax} with default values. Otherwise, if \code{data = NULL} these values will be set to \code{NULL} and you will have to specify these values. }
  \item{imType}{ (character)  \code{imType} specifies the type of image and thus the default values of the list. Default is \code{imType="normal"}. It is also implemented \code{imType="radon"}, i.e. that the image is a sinogram (radon transformed image). }
}
\details{
Set \code{listType="pet"} to generate a list with the following structure. It is assumed that \code{data} is a matrix and not equal to \code{NULL}.
  \tabular{rl}{
  Description \tab Contains a short description of the image. The maximal number  \cr
              \tab of characters has to be 80. Defaults to \code{Description="}\eqn{\backslash 0}{\0}\code{"}. \cr
  SignalDim   \tab Contains the number of rows and columns of the matrix. \cr
              \tab Defaults to \code{SignalDim=dim(data)}. \cr
  XYmin       \tab Leftmost coordinate and lowest coordinate in the original \cr
              \tab image. If \code{imType="normal"} then defaults to \code{XYmin=} \cr
              \tab \code{-0.5*(dim(data)-1)} or if \code{imType="radon"} then \cr
              \tab  defaults to \code{XYmin=c(0, -0.5*(ncol(data)-1))}. \cr
  DeltaXY     \tab  Quantization steps in the original image in \eqn{x} and \eqn{y} direction. \cr
    \tab If \code{imType="normal"} then defaults to \code{DeltaXY=c(1,1)} \cr 
    \tab or if \code{imType="radon"} then defaults to \code{DeltaXY=} \cr
    \tab \code{c(pi/(nrow(data)),1)}. \cr
  }
          
Set \code{listType="fif"} to generate a list, with the following structure. It is assumed that \code{data} is a matrix and not equal to \code{NULL}.
\tabular{rl}{
  FIFIdType \tab ID used to restore FIF:17737:'\eqn{\backslash 0}{\0}''\eqn{\backslash 0}{\0}''E''I'. Defaults to \cr
            \tab \code{FIFIdType=17737}. \cr
  FileName  \tab   Name used for saving/restoring this image. Defaults to \cr
            \tab  \code{FileName="}\eqn{\backslash 0}{\0}\code{"}. \cr
  Description \tab See above.  \cr
  Date      \tab Date (YYYY-MM-DD). Defaults to \code{Date=} \cr 
            \tab \code{as.character(Sys.Date())}. \cr
  SignalDim \tab See above.  \cr
  ArrayType \tab Defines the format number. 1 for Real and 2 for Complex. \cr
            \tab Complex matrices are determined by \eqn{A=a_{i,j}} \cr
            \tab where \eqn{a_{i,2n}} is the imaginary part to the real value \cr
            \tab \eqn{a_{i,2n-1}}, \eqn{i=1,\ldots,M}, \eqn{j=1,\ldots,2N}, \eqn{n=1,\ldots,N}. \cr
            \tab Defaults to \code{ArrayType=1}. \cr
  XYmin     \tab See above.  \cr
  DeltaXY   \tab See above.  \cr
  SignalMinMax \tab The lowest and highest signal value in the matrix. Defaults \cr
            \tab to \code{SignalMinMax=c(min(data), max(data))}. \cr
  }
}
\value{
  Returns a list with above structure.
}

\author{ Joern Schulz \email{jschulz78@web.de}, Peter Toft.}

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{writeData}}, \code{\link{readData}} }

\examples{
A <- phantom()
PetList <- iniPetFifList(data = A, imType = "normal")
PetList
rm(A,PetList)
}

\keyword{math}
\keyword{smooth}
\concept{iradon}
\concept{radon}
