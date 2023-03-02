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

      print(setpar)
      


      rdata <- .C("Forward_C_orig_dim", 
                  rdata=matrix(0, nrow=ThetaSamples, ncol=RhoSamples), 
                  as.double(oData), 
                  as.double(setpar),
                  PACKAGE="PET")$rdata  

    z <- list(rData=rdata, 
                Header=list(SignalDim=c(ThetaSamples, RhoSamples), 
                            XYmin=c(ThetaMin, RhoMin), 
                            DeltaXY=c(DeltaTheta, DeltaRho)), 
                call=args )
    class(z) <- "pet"
    return(z)
}
