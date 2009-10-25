\name{readData}
\alias{readData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reload Saved Datasets in differnt Formats.}
\description{
  The function offers the possibility to read a ".txt", ".dat", ".dat.gz", ".pet" or ".fif" file. Furthermore the graphics formats ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" and ".jpeg" are supported.
}
\usage{
readData(inputfile, DebugLevel = "Normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfile}{ (character)  \code{inputfile} gives the name of the file to read, including the pathname to the file. The extension of the name of the file specifies the format. Currently ".txt", ".dat", ".dat.gz", ".fif", ".pet" and the graphic formats ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" and ".jpeg" are supported. See below to details to get more information about the formats.}
  \item{DebugLevel}{ (character)  This parameter controls the level of output. Defaults to \code{DebugLevel="Normal"} for a standard level output. Alternative implementations are \code{"Detail"} if it is desirable to show almost all output on screen or \code{"HardCore"} for no information at all. }
}

\details{
In the following different supported formats of \code{inputfile} are explained.

\bold{".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" and ".jpeg"} \cr
For using these graphic formats the R package \code{adimpro} is necessary. See there, to get more information about these formats.

\bold{".txt"} \cr
Reads files in ASCII format without header, which contain only the raw data. The R-routine \code{read.table} is used. Returns a matrix.
 
\bold{".dat" and ".gz"} \cr
Read files in ASCII format that contain a header with three or four rows followed by the data. Compared to ".txt" the advantage is that is possible to save and load higher-dimensional data. Furthermore important values are saved, that are needed for Radon and the inverse Radon Transformation. If the extension ".gz" is specified, then files in ".dat" format will be read which were additionally compressed by 'gzip'. \cr
Assuming the data are two-dimensional, then the header has to be only of the following form: \cr
 \cr
Description:  char (optional) \cr
SignalDim:  int int \cr
XYmin:  double  double \cr
DeltaXY:  double double \cr
... and after the header follow the raw data \cr

Explanation of the parameters:
\tabular{rl}{
Description: \tab Can be a short description of the data. It is possible \cr
\tab to leave out this parameter. In this case the first \cr
\tab line of the file begins with 'SignalDim'. \cr
SignalDim: \tab In case of 2-dim data 'SignalDim' is the number of \cr
\tab rows and columns of the array, but also a greater \cr
\tab dimension is permissible. \cr
XYmin: \tab Leftmost coordinate and lowest coordinate of the \cr
\tab data. \cr
DeltaXY: \tab Quantization steps in image, in \eqn{x} and \eqn{y}-direction. \cr
}
If the dimension of data is greater than two, the number of 'SignalDim', 'XYmin' and 'DeltaXY' should be increase appropriately.

\bold{.pet} \cr
Reads a picture in ".pet" format (raw float with header) from the file e.g. \code{inputfile = "FileName.pet"} (binary format). The ".pet"-file has to be of the following structure: The first part contains a header and the second part the data, e.g. an image. 
Only the values of the parameter are saved in the header and not the parameter names.
\tabular{rl}{
Description: \tab (80-byte character) Should be a short description \cr 
\tab of image. \cr
XSamples: \tab (integer) Number of rows of array. \cr
YSamples: \tab (integer) Number of columns of array. \cr
Xmin: \tab (float) Leftmost coordinate in original image. \cr
Ymin: \tab (float) Lowest coordinate in original image. \cr
DeltaX: \tab (float) Quantization steps in original image (x). \cr
DeltaY: \tab (float) Quantization steps in original image (y). \cr
}
The dimension of the data should be not greater than two. \cr
The ".pet" and the sequencing ".fif" format will be used in the \code{iradonIT} method to store iteration-steps and to read a references image.

\bold{.fif} \cr
Reads a picture in ".fif" format (raw float with header) from the file e.g. \code{inputfile = "FileName.fif"} (binary format). The ".fif"-file has to be of the following structure: The first part contains a header and the second part the data, e.g. an image. At first, the values of the following parameters are saved in the file:
\tabular{rl}{
FIFIdType: \tab (integer) ID used to restore FIF:17737:'\eqn{\backslash 0}{\0}''\eqn{\backslash 0}{\0}''E''I'. \cr
FileName: \tab (100-byte character) Name used for saving this image. \cr
Description: \tab See above. \cr
Date: \tab (10-byte character) Date (YYYY-MM-DD).  \cr
XSamples: \tab See above. \cr
YSamples: \tab See above. \cr
ArrayType: \tab (integer) Defines the format number: \eqn{1} for Real or \cr
\tab \eqn{2} for Complex. Complex matrices are determined by \eqn{A=a_{i,j}} \cr
\tab where \eqn{a_{i,2n}} is the imaginary part to the real value \cr
\tab \eqn{a_{i,2n-1}}, \eqn{i=1,\ldots,M}, \eqn{j=1,\ldots,2N}, \eqn{n=1,\ldots,N}. \cr
Xmin: \tab (float) See above. \cr
Ymin: \tab (float) See above. \cr
DeltaX: \tab (float) See above. \cr
DeltaY: \tab (float) See above. \cr
SignalMin: \tab (float) Lowest signal value in the array. \cr
SignalMax: \tab (float) Highest signal value in the array. \cr
}
After this header the data are following.

}
\value{
  In case of \code{inputfile="/.../FileName.txt"}, ".tif", ".tiff", ".pgm", ".ppm", ".png", ".pnm", ".gif", ".jpg" and ".jpeg" only a matrix (image) will be returned. \cr
  
  In case of \code{inputfile="/.../FileName.dat"}, ".dat.gz", ".pet" and ".fif" a list will be returned, with
  \item{Signal}{ The image or the array that was read. }
  \item{Header}{ Is a list of values specified through the file format. }
  In case of \code{inputfile="/.../FileName.dat"}, ".dat.gz" and ".pet" \code{$Header} contain the following list elements: Description (if ".dat" or ".dat.gz" is used then optional), SignalDim, XYmin and DeltaXY. \cr
  In case of \code{inputfile="/.../FileName.fif"} \code{$Header} contain the following list elements: FIFIdType, FileName, Description, Date,  SignalDim, ArrayType, XYmin, DeltaXY and SignalMinMax.

}
\references{Schulz, Joern, \emph{Diploma Thesis: Analyse von PET Daten unter Einsatz adaptiver Glaettungsverfahren}, Humboldt-Universitaet zu Berlin, Institut fuer Mathematik, 2006.}

\author{ Joern Schulz, \email{jschulz78@web.de}. }

%\note{ ~~further notes~~ }
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{writeData}}, \code{\link{iradon}}}

\examples{
P <- phantom()
writeData(P, "Phantom.pet")
P.new <- readData("Phantom.pet", DebugLevel = "Normal")
viewData(P.new$Signal)
rm(P, P.new)
}
\keyword{IO}% at least one, from doc/KEYWORDS
\keyword{file}% __ONLY ONE__ keyword per line
\concept{radon}
%\concept{}
