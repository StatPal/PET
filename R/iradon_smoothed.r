## Calling from C

iradon_smoothed <- function(rData, XSamples, YSamples, 
                   mode="FB", 
                   Xmin=-sqrt(0.5)*(ncol(rData)/XSamples)*0.5*(XSamples-1), 
                   Ymin=-sqrt(0.5)*(ncol(rData)/YSamples)*0.5*(YSamples-1), 
                   DeltaX=sqrt(0.5)*(ncol(rData)/XSamples), 
                   DeltaY=sqrt(0.5)*(ncol(rData)/YSamples),
                   InterPol=1,
                   FilterTyp="Hamming1",
                   oData=NULL,
                   DebugLevel="Normal", 
                   iniFile=NULL )
{
  args <- match.call()
  if (is.character(iniFile)){
      iniList <- readIni(iniFile, DebugLevel=DebugLevel)
  
      mode <- iniList$mode
      rData <- iniList$rData
      XSamples <- iniList$XSamples
      YSamples <- iniList$YSamples
      Xmin <- iniList$Xmin
      Ymin <- iniList$Ymin
      DeltaX <- iniList$DeltaX
      DeltaY <- iniList$DeltaY
      InterPol <- iniList$InterPol
      FilterTyp <- iniList$FilterTyp
      oData <- iniList$oData
      DebugLevel <- iniList$DebugLevel


  } else if (!is.null(iniFile))
      stop("'iniFile' must specified an INI-file or set to FALSE.")

 ################################################################
 #              Checking the input-parameter
 #
  if (!(is.matrix(rData)))
      stop("'rData' has to be of type 'matrix'.")

  if (as.integer(XSamples)!=XSamples){
      cat("WARNING : XSamples is not of type integer and is truncated to", (as.integer(XSamples)), "\n")
      }
  if (as.integer(YSamples)!=YSamples){
      cat("WARNING : YSamples is not of type integer and is truncated to", (as.integer(YSamples)), "\n")
      }
  if (as.integer(InterPol)!=InterPol){
      cat("WARNING : InterPol is not of type integer and is truncated to", (as.integer(InterPol)), "\n")
      }
  if (!(FilterTyp=="Ramp" || FilterTyp=="Hamming1" || FilterTyp=="Hamming2")){
      cat("WARNING: FilterTyp=''",FilterTyp,"'' is not supported. \n", sep="")
      cat("Default is used: FilterTyp=''Hamming1'' \n")
      FilterTyp<-"Hamming1"
      }
  if (!(DebugLevel=="Normal" || DebugLevel=="Detail" || DebugLevel=="HardCore")){
      cat("WARNING: DebugLevel=''",DebugLevel,"'' is not supported. \n", sep="")
      cat("Default is used: DebugLevel=''Normal'' \n")
      DebugLevel<-"Normal"
      }
    
  if (FilterTyp=="Ramp")
      FilterTypC <- "Ramp"
  else if (FilterTyp=="Hamming1")
      FilterTypC <- "Hanning"
  else if (FilterTyp=="Hamming2")
      FilterTypC <- "Hamming"

  if (mode=="CC" || mode=="CNC" || mode=="CBC" || mode=="CNF" || 
      mode=="CBF" || mode=="FB" || mode=="BF" ){

      rDataDim  <- dim(rData)
    
     # calling the C-routine "iradon" from iradon.c
      irData <-.C("iradon_smoothed_C", 
                 as.double(rData), 
                 irData=double(XSamples*YSamples), 
                 as.character(mode),
                 as.integer(InterPol), 
                 as.character(FilterTypC),
                 as.character(DebugLevel), 
                 as.double(Xmin), 
                 as.double(Ymin), 
                 as.double(DeltaX), 
                 as.double(DeltaY), 
                 as.integer(rDataDim[1]), 
                 as.integer(rDataDim[2]), 
                 as.integer(XSamples), 
                 as.integer(YSamples),
                 PACKAGE="PET")$irData
    
      irData <- scaleImage(matrix(irData,nrow=XSamples,ncol=YSamples, byrow=TRUE))
      z <- list(irData=irData, 
            Header=list(SignalDim=c(XSamples,YSamples), 
                        XYmin=c(Xmin, Ymin), 
                        DeltaXY=c(DeltaX,DeltaY)), 
            call=args)
      class(z) <- "pet"

  } else
      stop("The mode=",mode," is not supported.")

  return(z)
}
