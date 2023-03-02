Forward_R <- function(oData, 
                    mode="NN", 
					XYSamples=nrow(oData),
					XYmin=-0.5*(nrow(oData)-1),
					DeltaXY=1,
					ThetaSamples=181,
					RhoSamples=2*round(sqrt(sum((dim(oData))^2))/2)+1,
					ThetaMin=0,
				          RhoMin=-0.5*((2*round(sqrt(sum((dim(oData))^2))/2)+1)-1),
					DeltaTheta=pi/ThetaSamples,
					DeltaRho=(2*abs(RhoMin)+1)/RhoSamples)
{
    args <- match.call()
	 
	 ######################################### 
	 # checking parameter
	  #if (DeltaRho>(DeltaXY/sqrt(2)))
	  #    stop("DeltaRho should be less than DeltaXY/sqrt(2)")
	  if (!(is.matrix(oData)))
		  stop("'oData' has to be of type 'matrix'.")

	  setpar      <- matrix(0, nrow=9, ncol=1)
	  setpar[1,1] <- XYSamples
	  setpar[2,1] <- XYmin
	  setpar[3,1] <- DeltaXY
	  setpar[4,1] <- ThetaSamples 
	  setpar[5,1] <- DeltaTheta
	  setpar[6,1] <- RhoSamples
	  setpar[7,1] <- RhoMin
	  setpar[8,1] <- DeltaRho
	  setpar[9,1] <- ThetaMin
      
      
    irData <-.C("Forward_C_orig_dim", 
                as.double(oData), 
                dat=double(ThetaSamples*RhoSamples), 
                as.character(mode),
                as.integer(1),      ## What is interpol?
                as.character("Ramp"),  ## What is FilterTypC 
                as.character("Normal"), 
                as.double(ThetaMin), 
                as.double(RhoMin), 
                as.double(DeltaTheta), 
                as.double(DeltaRho), 
                as.integer(rDataDim[1]), 
                as.integer(rDataDim[2]), 
                as.integer(ThetaSamples), 
                as.integer(RhoSamples),
                PACKAGE="PET")

    irdat = irData$dat

    z <- list(rData=irdat,
            Header=list(SignalDim=c(ThetaSamples, RhoSamples), 
                        XYmin=c(ThetaMin, RhoMin), 
                        DeltaXY=c(DeltaTheta, DeltaRho)), 
        call=args)
    class(z) <- "pet"

    return(z)
}
