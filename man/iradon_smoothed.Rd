\name{iradon_smoothed}
\alias{iradon_smoothed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Iterative Inverse Radon Transformation with Smoothing}
\description{
The function implements iterative reconstruction techniques. 
}
\usage{
iradon_smoothed(rData, XSamples, YSamples, StartImage = "None", 
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

\examples{
#    
# Compare the results of iterative reconstruction method "EM" and 
# direct reconstruction method "FB"
#
\dontrun{
P <- phantom(design="B")
rP <- markPoisson(P, nSample=1600000 )
irP1 <- iradon_smoothed(rP$rData , nrow(P), ncol(P))
viewData(list(P, rP$rData, irP1$irData),
         list("Generated unnoisy Phantom", "Generated PET Data",
         "Direct rec.: mode='FB'"))
rm(irP1,P,rP)
}
}
\keyword{math}
\keyword{smooth}
\concept{iradonIT}
\concept{iradon}
\concept{iradon_smoothed}
\concept{radon}
\concept{positronen}
\concept{emission} 
\concept{tomographie}
\concept{pet}
\concept{sinogram}
