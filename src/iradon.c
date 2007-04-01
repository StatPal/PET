/***************************************************************************
[HEADER]

This program is the 'main' program, which contains routines to call the
different direct reconstruction routines. The mode of operation
is determined by the parameter sending from R.

The ordering of the following items are not important, but each require
on line.

double *InImage
\item[{\tt InImage}] Contain the sinogram image (for the reconstruction 
functions).
\item[{\tt OutImage}] If mode not equal 'Test' than OutImage contains the
reconstructed image, that is sending to R.
If mode='Test', measures of misfit can be made between the rData
and the oData. 'oData' contain a matrix with orginal data of the
radon tranformed data and in this case 'OutImage' countain this matrix.
\item[{\tt mode}] Specify the function the program should do. Currently
supported functions are
\begin{description}
\item[{\tt FB}] Filtered Backprojection.
\item[{\tt BF}] Filtering After Backprojection.
\item[{\tt CNF}] Central Slice. FFT based with Nearest Neighbor approximation.
\item[{\tt CBF}] Central Slice. FFT based with Bilinear Interpolation.
\item[{\tt CC}]  Using a variation of the Central Slice. By the use of nonlinear 
sampling of the Radon domain with the use of the Chirp-z algorithm, the 
interpolation can be reduced to simple one dimensional linear interpolation 
compared to normal CS.
\item[{\tt CNC}] Central Slice. Chirp-z based with Nearest Neighbor approximation.
\item[{\tt CBC}] Central Slice. Chirp-z based with Bilinear Interpolation.
\item[{\tt Test}] Mode to run all the avaliable reconstruction routines on the 
same image and compare with a L1Norm and L2Norm. Surrender value is "NULL".
\end{description}
\item[{\tt DebugLevel}] This parameter controls the level of output. The
parameter is mixed and overrules with the one used in the
Print-statements.
\begin{description}
\item[{\tt Detail}] Allmost all output is logged to screen.
\item[{\tt Normal}] Standard level of output to screen.
\item[{\tt HardCore}] No information at all.
\end{description}
\item[{\tt InterPol}] Interpolation level. Used by Filtered Backprojection.
\item[{\tt FilterTyp}] Filter type. Only used by Filtered Backprojection (FB).
\begin{description}
\item[{\tt Ramp}] Very good results by data without noise.
\item[{\tt Hamming1}] Generalized Hamming Filter with alpha=0.5 
\item[{\tt Hamming2}] Generalized Hamming Filter with alpha=0.54 
\end{description}
\item[{\tt Xmin}] The minimum x-position of the reconstructed image.
\item[{\tt Ymin}] The minimum y-position of the reconstructed image.
\item[{\tt DeltaX}] Sampling distance on the x-axis.
\item[{\tt DeltaY}] Sampling distance on the y-axis.
\item[{\tt XSamples}] Number of samples on the x-axis.
\item[{\tt YSamples}] Number of samples on the y-axis.
\end{description}


PT April 96
JS July 2006
***************************************************************************/

#include "iradon.h"
#include <R.h>
#include <Rdefines.h>


/***************************************************************************
[NAME]
DoTest

[SYNOPSIS]
DoTest();

[DESCRIPTION] 

Function to run all the avaliable reconstruction routines on the same
image. The $L_1$ and $L_2$ measures are calculated for each routine.

[USAGE]

{\tt DoTest();} 

[REVISION]
Nov. 94, JJJ and PAP
March 2006 Jörn Schulz.
***************************************************************************/
void DoTest(Image *Image1, Image *ref)
{
  Image *Image2,*ImageFB,*ImageBF,*ImageCNC,
        *ImageCBC,*ImageCC,*spectrum;
  float Tid;
  
  Image2=CopyImage(Image1);
   
  Print(_DNormal,"1: \n");
  Tid=clock();
  spectrum=CentralSliceNN(Image2);
  ImageCNC=IFFTSpectrum(spectrum);
  RenameImage(ImageCNC, "CNC reconstructed image.");
  Tid=(clock()-Tid)/(float)CLOCKS_PER_SEC;
  Print(_DNormal,"IRadon: Reconstruction was active for %.2f seconds\n",Tid) ;
  FreeImage(Image2);
  FreeImage(spectrum);
 
  Print(_DNormal,"2:                    \n");
  Image2=CopyImage(Image1);
  Tid=clock();
  spectrum=CentralSliceBL(Image2);
  ImageCBC=IFFTSpectrum(spectrum);
  RenameImage(ImageCBC, "CBC reconstructed image.");
  Tid=(clock()-Tid)/(float)CLOCKS_PER_SEC;
  Print(_DNormal,"IRadon: Reconstruction was active for %.2f seconds\n",Tid) ;
  FreeImage(Image2);
  FreeImage(spectrum);

  Print(_DNormal,"3:                    \n");
  Image2=CopyImage(Image1);
  Tid=clock();
  spectrum=CentralSliceCZ(Image2);
  ImageCC=IChirpSpectrum(spectrum);
  RenameImage(ImageCC, "CC reconstructed image.");
  Tid=(clock()-Tid)/(float)CLOCKS_PER_SEC;
  Print(_DNormal,"IRadon: Reconstruction was active for %.2f seconds\n",Tid) ;
  FreeImage(Image2); 
  FreeImage(spectrum);

  Print(_DNormal,"4:                    \n");
  Image2=CopyImage(Image1);
  Tid=clock();
  ImageBF=BackFilter(Image2);
  RenameImage(ImageBF, "BF reconstructed image.");
  NormImage(ImageBF,1.0,MeanValue(ref));
  Tid=(clock()-Tid)/(float)CLOCKS_PER_SEC;
  Print(_DNormal,"IRadon: Reconstruction was active for %.2f seconds\n",Tid) ;
  FreeImage(Image2);
  
  Print(_DNormal,"5:                    \n");
  Image2=CopyImage(Image1);
  Tid=clock();
  ImageFB=FilteredBack(Image2);
  RenameImage(ImageFB, "FB reconstructed image.");
  Tid=(clock()-Tid)/(float)CLOCKS_PER_SEC;
  Print(_DNormal,"IRadon: Reconstruction was active for %.2f seconds\n",Tid) ;
  FreeImage(Image2);
  
  ScaleImage(ImageFB);
  ScaleImage(ImageBF);
  ScaleImage(ImageCNC);
  ScaleImage(ImageCBC);
  ScaleImage(ImageCC);
  
  PrintStats(_DDetail, ImageFB);
  PrintStats(_DDetail, ImageBF);
  PrintStats(_DDetail, ImageCNC);
  PrintStats(_DDetail, ImageCBC);
  PrintStats(_DDetail, ImageCC);
  
  Print(_DNormal,"FB:  L1=%9.6f, L2=%9.6f \n",L1Norm(ref,ImageFB),L2Norm(ref,ImageFB));
  Print(_DNormal,"BF:  L1=%9.6f, L2=%9.6f \n",L1Norm(ref,ImageBF),L2Norm(ref,ImageBF));
  Print(_DNormal,"CNC: L1=%9.6f, L2=%9.6f \n",L1Norm(ref,ImageCNC),
	L2Norm(ref,ImageCNC));
  Print(_DNormal,"CBC: L1=%9.6f, L2=%9.6f \n",L1Norm(ref,ImageCBC),
        L2Norm(ref,ImageCBC));
  Print(_DNormal,"CC:  L1=%9.6f, L2=%9.6f \n",L1Norm(ref,ImageCC),
        L2Norm(ref,ImageCC));	  

  FreeImage(ImageFB);
  FreeImage(ImageBF);
  FreeImage(ImageCNC);
  FreeImage(ImageCBC);
  FreeImage(ImageCC);
  
  FreeImage(ref);
}


/***************************************************************************
[NAME]
Main Subroutine

[DESCRIPTION] 

This is the main subroutine function for the direct reconstruction
methods. The routine have to call from R through iradon. 
See above for all possibilities. 

It reads the specified parameters getting from R
and executes one of the following functions:

\begin{itemize}
\item CentralSlice.
\item FilteredBack.
\item BackFiltering.
\item Convert.
\item Trace.
\item ImageInfo.
\item Test.
\end{itemize}

[REVISION]
Oct. 94, JJJ and PAP\\
April 5, 96 PT More info\\
April 9, 96 PT Better check on Value\\
April 13, 96 PT Error in percent print if 0 sec.\\
March 2006, Jörn Schulz, Modification of the complete routine
***************************************************************************/
void iradon(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
{
  Image *NewImage, *InvNewImage, *spectrum, *RefImage;
  
  if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
  else DebugNiveau=_DNormal;
   
  Print(_DNormal,"\n");
  Print(_DNormal,"Start of inverse Radontransformation\n");
  Print(_DNormal,"------------------------------------\n");
  Print(_DNormal,"          iradon (ver 2.0)          \n");
  Print(_DNormal,"      Made by Jesper J. Jensen      \n");
  Print(_DNormal,"   Peter Philipsen and Peter Toft   \n");
  Print(_DNormal,"  Implemented in R by Joern Schulz  \n");
  Print(_DNormal,"------------------------------------\n");

  ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 
      
  // ==================================================================
  // initialization of radon-image
  NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
  RDoubleToImage(NewImage, InImage, *M, *N );
  InitImage(NewImage);
  
  
  // ==================================================================
  
  if (strequal(*mode,"CC")||strequal(*mode,"CNC")||strequal(*mode,"CBC")){
    // Performs reconstruction using Central Slice with either 
    // nearest neighbor or chirp-z interpolation, renames the 
    // image and saves the result. Calls \tc{CentralSliceNN}, 
    // \tc{CentralSliceBL} or \tc{CentralSliceCZ}. Inverts the 
    // spectrum using inverse two dimensional chirp-z transformation.
    Print(_DNormal,"Using Central Slice with Chirp-z interpolation\n");
    if (strequal(*mode,"CNC"))
      spectrum=CentralSliceNN(NewImage);
    else if (strequal(*mode,"CBC"))
      spectrum=CentralSliceBL(NewImage);
    else
      spectrum=CentralSliceCZ(NewImage);
        
    InvNewImage=IChirpSpectrum(spectrum);
    FreeImage(spectrum);
    ScaleImage(InvNewImage);
    PrintStats(_DDetail,InvNewImage);
    ImageToFloat(OutImage, InvNewImage);
    FreeImage(InvNewImage);
  }
  else if (strequal(*mode,"CNF")||strequal(*mode,"CBF")){
    // Performes reconstruction using Central Slice with either 
    // nearest neighbor or bi-linear interpolation. Calls 
    // \tc{CentralSliceNN}, \tc{CentralSliceBL} or \tc{CentralSliceCZ}. 
    // Uses the two dimensional FFT to invert the spectrum.
    if (strequal(*mode,"CNF")){
      Print(_DNormal,"Using Central Slice with FFT - Nearest Neighbor interpolation\n");
      spectrum=CentralSliceNN(NewImage);
      }
    else if (strequal(*mode,"CBF")){
      Print(_DNormal,"Using Central Slice with FFT - Bilinear interpolation\n");
      spectrum=CentralSliceBL(NewImage);
      }
    
    InvNewImage=IFFTSpectrum(spectrum);
    FreeImage(spectrum);
    ScaleImage(InvNewImage);
    PrintStats(_DDetail,InvNewImage);
    ImageToFloat(OutImage, InvNewImage);
    FreeImage(InvNewImage);
  }
  else if (strequal(*mode,"FB")){ 
    //Performes an reconstruction using Filtered Backprojection.
    Print(_DNormal,"Using Filtered Backprojection\n");
    InvNewImage=FilteredBack(NewImage);
    ScaleImage(InvNewImage);
    PrintStats(_DDetail,InvNewImage);
    ImageToFloat(OutImage, InvNewImage);
    FreeImage(InvNewImage);
    }
  else if (strequal(*mode,"BF")){ 
    // Performes reconstruction using Filtering after Backprojection.
    Print(_DNormal,"Using Filtering after Backprojection\n");
    InvNewImage=BackFilter(NewImage);
    ScaleImage(InvNewImage);
    PrintStats(_DDetail,InvNewImage);
    ImageToFloat(OutImage, InvNewImage);
    FreeImage(InvNewImage);
    }
  else if (strequal(*mode,"Test")){
    // initialization of the reference-image
    RefImage=NewFloatImage("Reference Image", *XSamples, *YSamples,_RealArray);
    RDoubleToImage(RefImage, OutImage, *XSamples, *YSamples );
    // run all the avaliable reconstruction routines on the same image.
    DoTest(NewImage, RefImage);
    }
  else {
    Print(_DHardCore,"Warning: Function not recognized: '%s'\n",*mode);
  }
  
  FreeImage(NewImage);
  Print(_DNormal,"return to R.          \n");

}

/***************************************************************************
[NAME]
loadFile Subroutine

[DESCRIPTION] 
loadFile(SEXP fileR, SEXP DebugLevelR)

The routione have to be calling from R and SEXP is an specific format to
handle data between R and C.  It's makes possible to read data with
{\tt .fif}, {\tt .pet} or {\tt .mat} format.

\begin{description}
\item[{\tt fileR}] Name of the file to write the image. The
program recognizes the following file extentions: 
\begin{description}
\item[{\tt fif}] Float image format which includes sampling parameters. 
\item[{\tt pet}] raw and outdated imageformat. 
\item[{\tt mat}] Matlab files with one matrix. 
\end{description}
\item[{\tt dataR}] The image, that should written under fileR.
\item[{\tt DebugLevelR}] This parameter controls the level of output. The
parameter is mixed and overrules with the one used in the
Print-statements.
\begin{description}
\item[{\tt Detaill}]  Allmost all output is logged to screen.
\item[{\tt Normal}] Standard level of output to screen.
\item[{\tt HardCore}] No information at all.
\end{description}
\end{description}

[REVISION]
July 2006, Jörn Schulz
***************************************************************************/

SEXP loadFile(SEXP fileR, SEXP DebugLevelR){

  int typeF, m, n, llist;
  int *p_fifiInt, *p_dimInt, *p_arraytypeInt;
  double *p_xyminDouble, *p_deltaxyDouble, *p_minmaxDouble;
  char *fileName, *DebugLevel, fileType[100], *chrprt;
  char *namesMAT[3]  = {"Signal","SignalDim","SignalMinMax"};
  char *namesPET[5]  = {"Signal", "Description", "SignalDim", "XYmin", "DeltaXY"};
  char *namesFIF[10] = {"Signal", "FIFIdType", "FileName", "Description", "Date", 
                        "SignalDim", "ArrayType", "XYmin", "DeltaXY", "SignalMinMax"};
  Image *OutImage;
  SEXP ans, list, list_names, 
       fifiInt, fileNameChar, desrciptionChar, dateChar, dimInt, arraytypeInt,
       xyminDouble, deltaxyDouble, minmaxDouble;

  // ---------------------------------------------------------------------
  // Evaluating DegbugLevelR
  if(isString(DebugLevelR)){
    PROTECT(DebugLevelR = AS_CHARACTER(DebugLevelR));
    // allocate memory
    DebugLevel = R_alloc(strlen(CHAR(STRING_ELT(DebugLevelR, 0))),sizeof(char));
    // ... and copy the character value
    strcpy(DebugLevel, CHAR(STRING_ELT(DebugLevelR, 0)));
    UNPROTECT(1);
           
	   if (strstr(DebugLevel,"Normal"))   DebugNiveau=_DNormal;
    else if (strstr(DebugLevel,"Detail"))   DebugNiveau=_DDetail;
    else if (strstr(DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
    else {
      Print(_DNormal, "Unknown DebugLevel. Default is used. \n");
      DebugNiveau=_DNormal;
    }
  } else {
      Print(_DNormal, "DebugLevel has to be of type character. Default is used. \n");
	DebugNiveau=_DNormal;
  }
  // ---------------------------------------------------------------------
      
  // ---------------------------------------------------------------------
  // Evaluating fileR
  if(isString(fileR)){
    PROTECT(fileR = AS_CHARACTER(fileR));
    // allocate memory
    fileName = R_alloc(strlen(CHAR(STRING_ELT(fileR, 0))),sizeof(char));
    // ... and copy the character value
    strcpy(fileName, CHAR(STRING_ELT(fileR, 0)));
    UNPROTECT(1);
    
    chrprt=strrchr(fileName,'.');
    if(chrprt) {
      strcpy(fileType, chrprt+1);
      *chrprt='\0';
    }
        
         if (strstr(fileType,"fif")) typeF = _FIF;
    else if (strstr(fileType,"pet")) typeF = _PET;
    else if (strstr(fileType,"mat")) typeF = _MAT;
    else Error("Unknown file type: '%s'", fileType);
    
    OutImage = ReadImage(fileName, typeF);
    
    if (strstr(fileType,"mat")){
    // ***************************************************
    // mat Files
      llist = (int) 3;
      
      PROTECT(ans = allocMatrix(REALSXP, OutImage->M, OutImage->N));
      for(n = 0; n < OutImage->N; n++){
        for(m = 0; m < OutImage->M; m++){
          REAL(ans)[m + n*(OutImage->M)] = (double) OutImage->Signal[m][n];
        }
      }
    
      // Creating a character string vector of the names attribute of the 
      // objects in out list
      PROTECT(list_names = allocVector(STRSXP,llist));
      for(n=0; n < llist; n++)
        SET_STRING_ELT(list_names,n,mkChar(namesMAT[n]));
            
	// Creating a list with 1 vector elements
      PROTECT(list = allocVector(VECSXP,llist));
      // attaching ans vector to list:
      SET_VECTOR_ELT(list, 0, ans);
      // attaching the vector names:
      setAttrib(list, R_NamesSymbol, list_names);
    
      UNPROTECT(3);
    
    } else if (strstr(fileType,"pet") ){
    // ***************************************************
    // fif - files
    // l <- readData("MegaTest.fif", fformat="fif", DebugLevel="Detail")
      llist = (int) 5;
      
      PROTECT(ans = allocMatrix(REALSXP, OutImage->M, OutImage->N));
      for(n = 0; n < OutImage->N; n++){
        for(m = 0; m < OutImage->M; m++){
          REAL(ans)[m + n*(OutImage->M)] = (double) OutImage->Signal[m][n];
        }
      }
    
      // Creating a character string vector of the names attribute of the 
      // objects in out list
      PROTECT(list_names = allocVector(STRSXP,llist));
      for(n=0; n < llist; n++)
        SET_STRING_ELT(list_names,n,mkChar(namesPET[n]));

      PROTECT(desrciptionChar = allocVector(STRSXP, 1));
      SET_STRING_ELT(desrciptionChar, 0, mkChar(OutImage->Description));            

      PROTECT(dimInt = NEW_INTEGER(2));
      p_dimInt = INTEGER_POINTER(dimInt);
      p_dimInt[0] = OutImage->M;
      p_dimInt[1] = OutImage->N;
    
      PROTECT(xyminDouble = NEW_NUMERIC(2));
      p_xyminDouble = NUMERIC_POINTER(xyminDouble);
      p_xyminDouble[0] = (double) OutImage->Xmin;
      p_xyminDouble[1] = (double) OutImage->Ymin;

      PROTECT(deltaxyDouble = NEW_NUMERIC(2));
      p_deltaxyDouble = NUMERIC_POINTER(deltaxyDouble);
      p_deltaxyDouble[0] = (double) OutImage->DeltaX;
      p_deltaxyDouble[1] = (double) OutImage->DeltaY;

      // Creating a list with 5 vector elements
      PROTECT(list = allocVector(VECSXP,llist));
      // attaching ans vector to list:
      SET_VECTOR_ELT(list, 0, ans);
      SET_VECTOR_ELT(list, 1, desrciptionChar);
      SET_VECTOR_ELT(list, 2, dimInt);
      SET_VECTOR_ELT(list, 3, xyminDouble);
      SET_VECTOR_ELT(list, 4, deltaxyDouble);
      // attaching the vector names:
      setAttrib(list, R_NamesSymbol, list_names);
    
      UNPROTECT(7);

    } else if (strstr(fileType,"fif") ){
    // ***************************************************
    // fif - files
    // l <- readData("MegaTest.fif", fformat="fif", DebugLevel="Detail")
      llist = (int) 10;
      
	PROTECT(ans = allocMatrix(REALSXP, OutImage->M, OutImage->ArrayType*OutImage->N));
	for(n = 0; n < OutImage->ArrayType*OutImage->N; n++){
	  for(m = 0; m < OutImage->M; m++){
	    REAL(ans)[m + n*(OutImage->M)] = (double) OutImage->Signal[m][n];
        }
      }
    
      // Creating a character string vector of the names attribute of the 
      // objects in out list
      PROTECT(list_names = allocVector(STRSXP,llist));
      for(n=0; n < llist; n++)
        SET_STRING_ELT(list_names,n,mkChar(namesFIF[n]));

      PROTECT(fifiInt = NEW_INTEGER(1));
      p_fifiInt = INTEGER_POINTER(fifiInt);
      p_fifiInt[0] = OutImage->FIFIdType;

      PROTECT(fileNameChar = allocVector(STRSXP, 1));
      SET_STRING_ELT(fileNameChar, 0, mkChar(OutImage->FileName));

      PROTECT(desrciptionChar = allocVector(STRSXP, 1));
      SET_STRING_ELT(desrciptionChar, 0, mkChar(OutImage->Description));            

      PROTECT(dateChar = allocVector(STRSXP, 1));
      SET_STRING_ELT(dateChar, 0, mkChar(OutImage->Date));

      PROTECT(dimInt = NEW_INTEGER(2));
      p_dimInt = INTEGER_POINTER(dimInt);
      p_dimInt[0] = OutImage->M;
      p_dimInt[1] = OutImage->N;
    
      PROTECT(arraytypeInt = NEW_INTEGER(1));
      p_arraytypeInt = INTEGER_POINTER(arraytypeInt);
      p_arraytypeInt[0] = OutImage->ArrayType;

      PROTECT(xyminDouble = NEW_NUMERIC(2));
      p_xyminDouble = NUMERIC_POINTER(xyminDouble);
      p_xyminDouble[0] = (double) OutImage->Xmin;
      p_xyminDouble[1] = (double) OutImage->Ymin;

      PROTECT(deltaxyDouble = NEW_NUMERIC(2));
      p_deltaxyDouble = NUMERIC_POINTER(deltaxyDouble);
      p_deltaxyDouble[0] = (double) OutImage->DeltaX;
      p_deltaxyDouble[1] = (double) OutImage->DeltaY;

      PROTECT(minmaxDouble = NEW_NUMERIC(2));
      p_minmaxDouble = NUMERIC_POINTER(minmaxDouble);
      p_minmaxDouble[0] = (double) OutImage->SignalMin;
      p_minmaxDouble[1] = (double) OutImage->SignalMax;

      // Creating a list with 6 vector elements
      PROTECT(list = allocVector(VECSXP,llist));
      // attaching ans vector to list:
      SET_VECTOR_ELT(list, 0, ans);
      SET_VECTOR_ELT(list, 1, fifiInt);
      SET_VECTOR_ELT(list, 2, fileNameChar);
      SET_VECTOR_ELT(list, 3, desrciptionChar);
      SET_VECTOR_ELT(list, 4, dateChar);
      SET_VECTOR_ELT(list, 5, dimInt);
      SET_VECTOR_ELT(list, 6, arraytypeInt);
      SET_VECTOR_ELT(list, 7, xyminDouble);
      SET_VECTOR_ELT(list, 8, deltaxyDouble);
      SET_VECTOR_ELT(list, 9, minmaxDouble);
      // attaching the vector names:
      setAttrib(list, R_NamesSymbol, list_names);
    
      UNPROTECT(12);
	
    }
     
  } else {
      Print(_DHardCore, "Warning: Input has to be of type character. \n");
	ans = R_NilValue;
	
	
	
  }
  // ---------------------------------------------------------------------
 
  return(list);

}

/***************************************************************************
[NAME]
writeFile Subroutine

[DESCRIPTION] 
writeFile(SEXP fileR, SEXP dataR, SEXP headerR, SEXP DebugLevelR)

The routione have to be calling from R and SEXP is an specific format to
handle data between R and C. It's makes possible to save data under
{\tt .fif}, {\tt .pet} or {\tt .mat} format.

\begin{description}
\item[{\tt fileR}] Name of the file to write the image. The
program recognizes the following file extentions: 
\begin{description}
\item[{\tt fif}] Float image format which includes sampling parameters. 
\item[{\tt pet}] raw and outdated imageformat. 
\item[{\tt mat}] Matlab files with one matrix. 
\end{description}
\item[{\tt dataR}] The image, that should written under fileR.
\item[{\tt headerR}] A list that contains the necessary data for header of the specific file format.
\item[{\tt DebugLevelR}] This parameter controls the level of output. The
parameter is mixed and overrules with the one used in the
Print-statements.
\begin{description}
\item[{\tt Detaill}]  Allmost all output is logged to screen.
\item[{\tt Normal}] Standard level of output to screen.
\item[{\tt HardCore}] No information at all.
\end{description}
\end{description}

[REVISION]
March 2006, Joern Schulz
July 2006, Joern Schulz
***************************************************************************/

SEXP writeFile(SEXP fileR, SEXP dataR, SEXP headerR, SEXP DebugLevelR){

  int typeF, M, N;
  char *fileName, *DebugLevel, fileType[100], *chrprt;
  Image *SaveImage;
  double *dInImage;
  //SEXP names = getAttrib(headerR, R_NamesSymbol);

  // ---------------------------------------------------------------------
  // Evaluating DegbugLevelR
  if(isString(DebugLevelR)){
    PROTECT(DebugLevelR = AS_CHARACTER(DebugLevelR));
    // allocate memory
    DebugLevel = R_alloc(strlen(CHAR(STRING_ELT(DebugLevelR, 0))),sizeof(char));
    // ... and copy the character value
    strcpy(DebugLevel, CHAR(STRING_ELT(DebugLevelR, 0)));
    UNPROTECT(1);
           
	   if (strstr(DebugLevel,"Normal"))   DebugNiveau=_DNormal;
    else if (strstr(DebugLevel,"Detail"))   DebugNiveau=_DDetail;
    else if (strstr(DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
    else {
      Print(_DNormal,"Unknown DebugLevel. Default is used. \n");
      DebugNiveau=_DNormal;
    }
  } else
      Print(_DNormal,"DebugLevel has to be of type character. Default is used. \n");
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  // Evaluating InImage and fileR
  if(isString(fileR)){
    PROTECT(fileR = AS_CHARACTER(fileR));
    // allocate memory
    fileName = R_alloc(strlen(CHAR(STRING_ELT(fileR, 0))),sizeof(char));
    // ... and copy the character value
    strcpy(fileName, CHAR(STRING_ELT(fileR, 0)));
    UNPROTECT(1);
    
    chrprt=strrchr(fileName,'.');
    if(chrprt) {
      strcpy(fileType, chrprt+1);
      *chrprt='\0';
    }
        
         if (strstr(fileType,"fif")) typeF = _FIF;
    else if (strstr(fileType,"pet")) typeF = _PET;
    else if (strstr(fileType,"mat")) typeF = _MAT;
    else Error("Unknown file type: '%s'", fileType);
    
    Print(_DDetail,"`OutFileName' entry is set to `%s'\n", fileName);        
    
    // protecting the image-data
    PROTECT(dataR = AS_NUMERIC(dataR));
    dInImage = NUMERIC_POINTER(dataR);
    UNPROTECT(1);

    // List from R with the values (the header)
    PROTECT(headerR = AS_VECTOR(headerR));
    
    if(isInteger(VECTOR_ELT(headerR, 0))){
	  M = INTEGER(VECTOR_ELT(headerR, 0))[0];
        N = INTEGER(VECTOR_ELT(headerR, 0))[1];
    }
    if(isReal(VECTOR_ELT(headerR, 0))){
        M = (int) REAL(VECTOR_ELT(headerR, 0))[0];
        N = (int) REAL(VECTOR_ELT(headerR, 0))[1];
    }
    
    if(strstr(fileType,"fif")){
        if(isInteger(VECTOR_ELT(headerR, 7))){
	      SaveImage=NewFloatImage(fileName, M, N, INTEGER(VECTOR_ELT(headerR, 7))[0]);
        } else {
	     SaveImage=NewFloatImage(fileName, M, N, (int) REAL(VECTOR_ELT(headerR, 7))[0]);
        }
    } else
        SaveImage=NewFloatImage(fileName, M, N, _RealArray);    
    RDoubleToImage(SaveImage, dInImage, M, N);
    
    
    //if (strstr(fileType,"mat")){
    // ***************************************************
    // mat Files
    //  SaveImage->SignalMin = (float) REAL(VECTOR_ELT(listR, 2))[0];
    //  SaveImage->SignalMax = (float) REAL(VECTOR_ELT(listR, 2))[1];
    //
    //} else 
    if (strstr(fileType,"pet") || strstr(fileType,"fif")) {
    // ***************************************************
    // pet & fif - files
      strcpy(SaveImage->Description, CHAR(STRING_ELT(VECTOR_ELT(headerR, 1), 0)));
      SaveImage->Xmin = (float) REAL(VECTOR_ELT(headerR, 2))[0];
      SaveImage->Ymin = (float) REAL(VECTOR_ELT(headerR, 2))[1];
      SaveImage->DeltaX = (float) REAL(VECTOR_ELT(headerR, 3))[0];
      SaveImage->DeltaY = (float) REAL(VECTOR_ELT(headerR, 3))[1];

    if (strstr(fileType,"fif")) {
    // ***************************************************
    // fif - files
      Print(_DDetail,"Special fif-parameter. \n");
      SaveImage->SignalMin = (float) REAL(VECTOR_ELT(headerR, 4))[0];
      SaveImage->SignalMax = (float) REAL(VECTOR_ELT(headerR, 4))[1];
      strcpy(SaveImage->Date, CHAR(STRING_ELT(VECTOR_ELT(headerR, 5), 0)));
      // FIFIdType
	if(isInteger(VECTOR_ELT(headerR, 6))){
        SaveImage->FIFIdType = INTEGER(VECTOR_ELT(headerR, 6))[0];
      }
      if(isReal(VECTOR_ELT(headerR, 6))){
        SaveImage->FIFIdType = (int) REAL(VECTOR_ELT(headerR, 6))[0];
      }
     // ArrayType
     if(isInteger(VECTOR_ELT(headerR, 7))){
        SaveImage->ArrayType = INTEGER(VECTOR_ELT(headerR, 7))[0];
      }
      if(isReal(VECTOR_ELT(headerR, 7))){
        SaveImage->ArrayType = (int) REAL(VECTOR_ELT(headerR, 7))[0];
      }
     //strcpy(SaveImage->FileName, CHAR(STRING_ELT(VECTOR_ELT(headerR, 8), 0)));
     /* // Type
      if(isInteger(VECTOR_ELT(headerR, 9))){
        SaveImage->Type = INTEGER(VECTOR_ELT(headerR, 9));
      }
      if(isReal(VECTOR_ELT(headerR, 9))){
        SaveImage->Type = (int) REAL(VECTOR_ELT(headerR, 9));
      } */
    }
            
    UNPROTECT(1);
    
    if (DebugNiveau==_DDetail){
    if ( strstr(fileType,"fif")){
      Print(_DDetail,"`FIFIdType' entry is set to %i\n",SaveImage->FIFIdType);
      Print(_DDetail,"`FileName' entry is set to `%s'\n", SaveImage->FileName );
    }
    if ( strstr(fileType,"dat") || strstr(fileType,"fif")){
      Print(_DDetail,"`Description' entry is set to `%s'\n", SaveImage->Description );
    }
    if ( strstr(fileType,"fif")){
      Print(_DDetail,"`Date' entry is set to `%s'\n", SaveImage->Date );
      // Print(_DDetail,"`Type' entry is set to %i\n",SaveImage->Type);
    }
    Print(_DDetail,"`M' entry is set to %i\n",SaveImage->M);
    Print(_DDetail,"`N' entry is set to %i\n",SaveImage->N);
    if (strstr(fileType,"fif")){
      Print(_DDetail,"`ArrayType' entry is set to %i\n",SaveImage->ArrayType);
    }
    if ( strstr(fileType,"dat") || strstr(fileType,"fif")){
      Print(_DDetail,"`Xmin' entry is set to %f\n",SaveImage->Xmin);
      Print(_DDetail,"`Ymin' entry is set to %f\n",SaveImage->Ymin);
      Print(_DDetail,"`DeltaX' entry is set to %f\n",SaveImage->DeltaX);
      Print(_DDetail,"`DeltaY' entry is set to %f\n",SaveImage->DeltaY);
    }
    if (strstr(fileType,"fif")){
      Print(_DDetail,"`SignalMin' entry is set to %f\n",SaveImage->SignalMin);
      Print(_DDetail,"`SignalMax' entry is set to %f\n",SaveImage->SignalMax);
    }
    }
    }
    if (strstr(fileType,"mat")){
      Print(_DDetail,"`M' entry is set to %i\n",SaveImage->M);
      Print(_DDetail,"`N' entry is set to %i\n",SaveImage->N);
      UNPROTECT(1);
    }

    WriteImage(SaveImage, typeF);
    FreeImage(SaveImage);
      
  } else Print(_DHardCore,"Warning: Filename and path has to be of type character. \n");
  // ---------------------------------------------------------------------

  return(R_NilValue);
}
