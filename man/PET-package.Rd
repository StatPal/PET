\name{PET-package}
\alias{PET-package}
\docType{package}
\title{Simulation and Reconstruction of PET Images}
\description{This package implements different analytic or direct and 
iterative reconstruction methods of Peter Toft. It also
offer the possibility to simulate PET data.}
\details{
\tabular{ll}{
Package: \tab PET\cr
Version: \tab 0.4.7\cr
Date: \tab 2009-10-25\cr
License: \tab GPL (>=2)\cr
Copyright: \tab 2006-2009 Joern Schulz and Peter Toft \cr
}

Index:
\preformatted{
cutMatrix               Cut matrix.
hough                   Hough transformation.
iniPetFifList           Initializes a header for PET or FIF file.
iradonIT                Iterative inverse radon transformation.
iradon                  Direct Inverse radon transformation.
markPoisson             Marked poisson process.
norm                    L1 and L2 norm.
partEllipse             Generates a fragment of an ellipse.
phantom                 Creation of a phantom.
radon                   Radon transformation.
readData                Reload saved datasets in differnt formats.
readSifData             Reloaded system matrix.
rotate                  Rotates data.
scaleImage              Scales an image.
viewData                Display several images.
writeData               Write datasets in different formats.
}
}
\author{
Joern Schulz <jschulz78@web.de>,
Peter Toft <pto@imm.dtu.dk>,
Jesper James Jensen <jjj@oedan.dk>,
Peter Philipsen <pap@imm.dtu.dk>

Maintainer: Joern Schulz <jschulz78@web.de>
}
\keyword{package}
