## New poisson data, directly susing the radon function.
# Author: Subrata Pal
# Here, the total number of case would be the nSample, 
# not the final number of cases.
## 
newPoisson <- function(oData, nSample=200000,
                    mode="NN", 
                    XYSamples=nrow(oData),
                    XYmin=-0.5*(nrow(oData)-1),
                    DeltaXY=1,
                    ThetaSamples=181,
                    RhoSamples=2*round(sqrt(sum((dim(oData))^2))/2)+1,
                    ThetaMin=0,
                    RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1),  
                    # Defaiult does not match with markPoisson
                    DeltaTheta=pi/ThetaSamples,
                    DeltaRho=(2*abs(RhoMin)+1)/RhoSamples){

    Forward_transform <- radon(oData, mode, XYSamples, XYmin, DeltaXY, ThetaSamples, RhoSamples, 
                                ThetaMin, RhoMin, DeltaTheta, DeltaRho)$rData
    Forward_transform[Forward_transform<0] <- 0
    Forward_transform <- Forward_transform/sum(Forward_transform)
    param <- nSample * Forward_transform

    array(rpois(n=length(param), lambda=param), dim=dim(param))
}
