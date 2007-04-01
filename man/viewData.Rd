\name{viewData}
\alias{viewData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Display Several Images.}
\description{
  This function prints 1 up to 8 images in a new window on the display or in a current window \code{curWindow}. The function uses the R function \code{image} to display the images. Note the conditions of this function.
}
\usage{
viewData(Data, Title = NULL, curWindow = TRUE, 
         colors = gray((0:255)/255), s = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{ Has to be a matrix or a list of matrices. One, two, up to eight image are possible. }
  \item{Title}{ A character or a list of titles, that specify the title of the images. Note if \code{Data} and \code{Title} are of type list then \code{length(Data) == length(Title)}. Defaults to \code{Title=NULL} (no title is uses).}
  \item{curWindow}{ If \code{curWindow=TRUE} or \code{curWindow=i} then the current active window is used or the specified window \code{i}. In these cases the size of the corresponding window is used. Otherwise if \code{curWindow=FALSE} then a new window is open. Defaults to \code{curWindow=TRUE}. }
  \item{colors}{ Corresponds to the parameter \code{col} of the function \code{image}. Defaults to \code{colors = gray((0:255)/255)}. }
  \item{s}{ Scaling parameter for the size of a new window. Defaults to \code{s=1}. }
}

\value{
  Returns the current window.
}

\author{ Joern Schulz, \email{jschulz78@web.de}. }

\examples{
P <- phantom()
R <- markPoisson(P)
viewData(list(P, R$Data, R$rData), list("Phantom", 
         "Marked Poisson Data", "Simulated PET Data"))
rm(P,R)
}

\keyword{print}
\concept{PET}
\concept{radon}
