\name{writeData}
\alias{writeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write Datasets in different Formats}
\description{
  The function offers the possibility to write ".txt", ".dat", ".dat.gz", ".pet" or ".fif" files. Furthermore the graphic formats ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" and ".jpeg" are supported. The function was especially written to handle sinograms and datasets, that are processed with the functions \code{radon}, \code{iradon} and \code{iradonIT} of the R package 'PET'.
}
\usage{
writeData(data, outputfile,  fileHeader = NULL, imType = "normal",
          DebugLevel = "Normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ (matrix)  The data to be written out. }
  \item{outputfile}{ (character)  \code{outputfile} naming the file to write to, including the pathname. The path has to be given relatively to the working-directory of your R-session or it contains the full path of the file. The extension of the filename specifies the format. Currently, ".txt", ".dat", ".dat.gz", ".fif", ".pet" and the graphic formats ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" and ".jpeg" are supported. See below to details to get more information about these formats.}
  \item{fileHeader}{ (list)  If the file format "is .dat", ".dat.gz", ".pet" or ".fif", then \code{fileHeader} is necessary. If in these cases \code{fileHeader=NULL}, default values will be generated dependent on \code{imType}. See below to details to get more information. }
  \item{imType}{ (character)  \code{imType} specifies the type of \code{data} and thus the default values of \code{fileHeader}. Default is \code{imType="normal"}. Additional implementation is \code{imType="radon"}, i.e. that the image is a sinogram (Radon transformed image). }
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="Normal"} for a standard level output. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"HardCore"} for no information at all. }
}
\details{
In the following different supported formats of \code{inputfile} are explained. For ".dat", ".dat.gz", ".pet" and ".fif" the default values dependent on \code{imType} are defined.

\bold{".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" and ".jpeg"} \cr
For using these graphic formats the R package \code{adimpro} is necessary. See there, to get more information about these formats.

\bold{".txt"} \cr
Writes files in ASCII format without header, which contain only the raw data. The R-routine \code{write.table} is used.
 
\bold{".dat" and ".gz"} \cr
Write files in ASCII format where files contain a header with three or four rows followed by the data. Compared to ".txt" the advantage is that is possible to save and load higher-dimensional data. Furthermore important values are saved, that are needed for Radon and the inverse Radon transformation. If the extension ".gz" is specified, then files in ".dat" format will be stored  which are additionally compressed by 'gzip'. \cr
Assuming the data are two-dimensional, at first \code{fileHeader} will be written out. Default values will be defined (dependent on \code{imType}) for not specified parameters. The file has only the following form: \cr
 \cr
Description:  char (optional) \cr
SignalDim:  int int \cr
XYmin:  double  double \cr
DeltaXY:  double double \cr
... and after the header \code{data} will be written out \cr

Explanation of the parameters:
\tabular{rl}{
Description: \tab Can be a short description of the data. It is possible \cr
\tab to leave out this parameter. In this case the first \cr
\tab line of the file begins with 'SignalDim'. \cr
SignalDim: \tab In case of 2-dim data, 'SignalDim' is the number of \cr
\tab rows and columns of the array, but also a greater \cr
\tab dimension is permissible. \cr
XYmin: \tab Leftmost coordinate and lowest coordinate of the \cr
\tab data. \cr
DeltaXY: \tab Quantization steps in image, in x and y-direction. \cr
}
Defaults to:
\tabular{lcc}{
 \tab \code{imType="normal"} \tab \code{imType="radon"} \cr
 \tab  \tab \cr
SignalDim: \tab \code{dim(data)} \tab \code{dim(data)} \cr
XYmin: \tab \code{-0.5*(dim(data)-1)} \tab \code{c(0,-0.5*(ncol(data)-1))} \cr
DeltaXY: \tab \code{c(1,1)} \tab \code{c(pi/(nrow(data)),1)} \cr
}
If the dimension of data is greater than two, the number of 'SignalDim', 'XYmin' and 'DeltaXY' should be increased appropriately. Note in this case you must define the \code{fileHeader}.

\bold{".pet"} \cr
Writes a picture in ".pet" format (raw float with header) to the file e.g. \code{inputfile = "FileName.pet"} (binary format). The ".pet"-file has to be of the following structure: The first part contains a header and the second part the data, e.g. an image. 
Only the values of the parameters are saved in the header and not the parameter names.
\tabular{rl}{
Description: \tab (80-byte character) Should be a short description of \cr 
\tab image. Defaults to \code{Description="}\eqn{\backslash 0}{\0}\code{"}. \cr
XSamples: \tab (integer) Number of rows of array. Defaults \cr
 \tab to \code{XSamples=nrow(data)}. \cr
YSamples: \tab (integer) Number of columns of array. Defaults \cr
 \tab to \code{YSamples=ncol(data)}. \cr
Xmin: \tab (float) Leftmost coordinate in original image. \cr
 \tab If \code{imType="normal"} then defaults to \code{Xmin=}\cr
 \tab \code{-0.5*(nrow(data)-1)} or if \code{imType="radon"} \cr
 \tab then \code{Xmin=0}. \cr
Ymin: \tab (float) Lowest coordinate in original image. \cr
 \tab Defaults to \code{Ymin=-0.5*(ncol(data)-1)}. \cr
DeltaX: \tab (float) Quantization steps in original image (x). \cr
 \tab If \code{imType="normal"} then defaults to \code{DeltaX=1} \cr
 \tab or if \code{imType="radon"} then \code{DeltaX=pi/(nrow(data))}. \cr
DeltaY: \tab (float) Quantization steps in original image (y). \cr
 \tab Defaults to \code{DeltaY=1}. \cr
}
The dimension of the data should be not greater than two. \cr
The ".pet" and the sequencing ".fif" format will be used in the \code{iradonIT} method to store iteration-steps and to read a references image.

\bold{".fif"} \cr
Writes a picture in ".fif" format (raw float with header) to the file e.g. \code{inputfile = "FileName.fif"} (binary format). The ".fif"-file has to be of the following structure: The first part contains a header and the second part the data, e.g. an image. At first the values of the following parameters are saved in the file:

\tabular{rl}{
FIFIdType: \tab (integer) ID used to restore FIF:17737:'\eqn{\backslash 0}{\0}''\eqn{\backslash 0}{\0}''E''I'. \cr
 \tab Defaults to \code{FIFIdType=17737}. \cr
FileName: \tab (100-byte character) Name used for saving this image. \cr
 \tab Defaults to \code{FileName=outputfile}. \cr
Description: \tab See above. \cr
Date: \tab (10-byte character) Date (YYYY-MM-DD). Defaults to \cr 
 \tab \code{Date=as.character(Sys.Date())}. \cr
XSamples: \tab See above. \cr
YSamples: \tab See above. \cr
ArrayType: \tab (integer) Defines the format number: \eqn{1} for Real or \cr
\tab \eqn{2} for Complex. Complex matrices are determined by \eqn{A=a_{i,j}} \cr
\tab where \eqn{a_{i,2n}} is the imaginary part to the real value \cr
\tab \eqn{a_{i,2n-1}}, \eqn{i=1,\ldots,M}, \eqn{j=1,\ldots,2N}, \eqn{n=1,\ldots,N}. \cr \tab Defaults to \code{ArrayType=1}. \cr
Xmin: \tab (float) See above. \cr
Ymin: \tab (float) See above. \cr
DeltaX: \tab (float) See above. \cr
DeltaY: \tab (float) See above. \cr
SignalMin: \tab (float) Lowest signalvalue in array. Defaults \cr
 \tab to \code{SignalMax=max(data)}. \cr
SignalMax: \tab (float) Highest signalvalue in array. Defaults \cr
 \tab to \code{SignalMax=max(data)}. \cr
}

After this header \code{data} are following.
}

\value{
  \code{writeData} creates a file, but returns nothing.
}

\references{Schulz, Joern, \emph{Diploma Thesis: Analyse von PET Daten unter Einsatz adaptiver Glaettungsverfahren}, Humboldt-Universitaet zu Berlin, Institut fuer Mathematik, 2006.}

\author{ Joern Schulz, \email{jschulz78@web.de}. }

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{readData}}, \code{\link{radon}}}

\examples{
\dontrun{
P <- phantom(addIm="blurred1")
rP <- radon(P, RhoSamples=401)
irP <- iradon(rP$rData, 257, 257)
# Saving 'P' as a 'jpeg'-file
writeData(P, "Phantom.jpeg")
# Saving 'rP' as a 'pet'-file
writeData(rP$rData, "RadonPhantom.pet", fileHeader=rP$Header)
# Saving 'irP' as a 'dat.gz'-file'
writeData(irP$irData, "RecPhantom.dat.gz", fileHeader=irP$Header)
rm(irP,P,rP)

Image1 <- readData("Phantom.jpeg")
Image2 <- readData("RadonPhantom.pet")
Image3 <- readData("RecPhantom.dat.gz")
viewData(list(Image1,Image2$Signal,Image3$Signal))
rm(Image1,Image2,Image3)
}

}
\keyword{IO}   % at least one, from doc/KEYWORDS
\keyword{file} % __ONLY ONE__ keyword per line
%\concept{}
