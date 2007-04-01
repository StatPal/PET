/*****************************************************************
[HEADER]

The ``it'' program implements three of the major iterative
reconstruction techniques: ART, EM, and LSCG, but with the basic
framework  provided virtually all iterative reconstruction algorithms
based on matrix operations can easily be implemented. The program uses
many of the same elements as ``iradon'' and provides several new for
matrix operations relevant for iterative reconstruction.


\begin{description}
\item[{\tt InImage}] Contain the sinogram image (for the reconstruction 
functions).
\item[{\tt OutImage}] If mode not equal 'Test' than OutImage contains the
reconstructed image, that is sending to R.
\item[{\tt RefFileName}] Error measures can be made between the RefFile and 
the OutFile.
\item[{\tt StartImageTrue}] If 1 then StartImage is used.
\item[{\tt StartImage}] If provided the initial guess on the
reconstructed image is taken from this image, else a constant solution
will be assumed from start.
\item[{\tt mode}] The iterative reconstruction method used:
Currently three are available: ART, EM and CG (LSCG).
\item[{\tt UseFast}] An integer. If 1 then fast reconstruction is used,
where the system matrix is stored (if KernelFileSave 1) using sparse 
techniques, else slower but memory efficient reconstruction is used.
\item[{\tt RadonKernel}] The type of kernel used to model the system
matrix. Currently available are the following methods
\begin{description}
  \item[{\tt NN}] Two-level Nearest Neighbor approximation. (Memory consuming).
  \item[{\tt RNN}] Ray driven Nearest Neighbor discrete Radon transform based (Very fast with small system matrix).
    \item[{\tt RL}] Ray driven Linear Interpolation discrete Radon transform
    based (Fast with small system matrix).
    \item[{\tt P1}] Method based on Radon transformation of square with
    pre-guidance (slow but good).
    \item[{\tt P1}] Method based on Radon transformation of square with
    no pre-guidance (slower but better).
    \item[{\tt SINC}] Sinc interpolation methods in the image and analytically
    Radon of that. (Very slow).
\end{description}
\item[{\tt IterationType}] For ART: If set to 'cyclical' a cyclical
selection of the row index is used, and else if set to 'random' a random 
selection is chosen.
\item[{\tt Iterations}] For EM and CG the number of iterations before the
iteration ends. For ART the number of full iterations, i.e., divided
by the number of rows in the system matrix.
\item[{\tt SaveIterations}] If set to 1 the current solution will be saved
after each corresponding iteration.
\item[{\tt SaveIterationsName}] The name of the file where the iterations 
are saved. Supported are files of type pet and fif.
\item[{\tt LowestALevel}] If fast reconstruction is use this the matrix
elements are truncated to this level relative to the sampling distance
of x, $\Delta x=$DeltaX.
\item[{\tt ConstrainMin}] After each iteration the solution in each sample
will forced above this limit.
\item[{\tt ConstrainMax}] After each iteration the solution in each sample
will forced below this limit. For both limit it is assumed that
negative limits imply that the feature is not used.
\item[{\tt Alpha}] For ART the initial update weight.
\item[{\tt Beta}] For ART the multiplicative change to the weight factor,
which should be less than one.
\item[{\tt Regularization}] If set to 1 and using fast reconstruction rows
will be appended to the system matrix with a simple Laplace
operator. A weight factor should also be incorporated.
\item[{\tt KernelFileSave}] Specify, whether system is saved or not (1 or 0).
\item[{\tt KernelFileName}] If using fast reconstruction the system matrix
will be saved and restored with this sif-name. If the system matrix read
is incopatible with the sampling parameters, then a new will be
generated.
\item[{\tt ThetaSamples}] Number of angular samples $theta$ in the sinogram.
\item[{\tt ThetaMin}] Start of the angular sampling (should be set to zero.).
\item[{\tt DeltaTheta}] Angular sampling distance (should be set to
pi/ThetaSamples).
\item[{\tt RhoSamples}] Number of samples in the sinogram $R$ in the
$\rho$-direction.
\item[{\tt DeltaRho}] Sampling distance in $\rho$, i.e., $\Delta \rho$
\item[{\tt RhoMin}] Start of sampling positions in $\rho$ (should be set to
$-\Delta\rho\frac{R-1}{2}$.
\item[{\tt Xmin}] The minimum x-position of the reconstructed image.
\item[{\tt Ymin}] The minimum y-position of the reconstructed image.
\item[{\tt DeltaX}] Sampling distance on the x-axis.
\item[{\tt DeltaY}] Sampling distance on the y-axis.
\item[{\tt XSamples}] Number of samples on the x-axis.
\item[{\tt YSamples}] Number of samples on the y-axis.
\item[{\tt OverSamp}] Use oversampling - Using squared no. of samples
\item[{\tt DebugLevel}] This parameter controls the level of output. The
parameter is mixed and overrules with the one used in the
Print-statements.
\begin{description}
\item[{\tt Detail}] Allmost all output is logged to screen.
\item[{\tt Normal}] Standard level of output to screen.
\item[{\tt HardCore}] No information at all.
\end{description}
\end{description}

April 96, JJJ and PT
July 2006, J.Schulz
*****************************************************************/

#include "it.h"
#include <R.h>

int DebugNiveau;
FILE *LogFile;
char LogFileName[100];
float multtemp;
INI IniFile;
float start,slut;
SparseMatrix *AMatrix, *AMatrix2;
Image *MyImage;
itINItype itINI;

/***************************************************************
[NAME]
SaveIteration

[SYNOPSIS]
void SaveIteration(Vector *MyVector, 
                   int iteration, 
                   char *filename)

[DESCRIPTION] 

This function is used to save a series of images while iterating. The
function saves an image with the specified outputname + a number
indicating which iteration the picture represents.

[USAGE]
{\tt SaveIteration(TestVector, 123, "testImage");}

Saves the vector {\tt TestVector} as the image `{\tt testimage.123.pet}'.

[REVISION]
Dec. 94, JJJ
July 2006, Joern Schulz
***************************************************************/
void SaveIteration(Vector *MyVector, int iteration, char *filename)
{
  char outfilename[100];
  char itstr[10];
  Image *NewImage;

  strcpy(outfilename,filename);
  sprintf(itstr,".%d",iteration);  
  strcat(outfilename,itstr);

  NewImage=VectorToImage(MyVector,itINI.XSamples,itINI.YSamples);

  RenameImage(NewImage,outfilename);
  NewImage->DeltaX=itINI.DeltaX;
  NewImage->DeltaY=itINI.DeltaY;
  NewImage->Xmin=itINI.Xmin;
  NewImage->Ymin=itINI.Ymin;
  ScaleImage(NewImage);
  WritePET(NewImage);
  FreeImage(NewImage);
}


/***************************************************************
[NAME]
main

[SYNOPSIS]
void main(int argc, char *argv[])

[DESCRIPTION]
This is the main function which controls the program.

[USAGE]
{\tt it test.ini}

Starts the main program with the parameteres specified in {\tt test.ini}.

[REVISION]
March 96, JJJ and PT\\
April 2, 96 PT (Moved last call to clock - error in SGI CC)
July 2006,  Modification of the complete routine for the R calling.
***************************************************************/
void it(double *InImage, double *OutImage, int *StartImageTrue, double *StartImage, char **mode, int *UseFast, char **RadonKernel, char **IterationsType, int *Iterations, int *SaveIterations, char **SaveIterationsName, double *LowestALevel, double *ConstrainMin, double *ConstrainMax, double *Alpha, double *Beta, double *Regularization, int *KernelFileSave, char **KernelFileName, char **RefFileName, int *ThetaSamples, double *ThetaMin, double *DeltaTheta, int *RhoSamples, double *RhoMin, double *DeltaRho, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *XSamples, int *YSamples, int *OverSamp, char **DebugLevel )
{
  int RealTid1,RealTid2,n;
  float Tid, tempsuma, tempsumb, mean;
  char Value[100];
  Vector *xvector, *tempvector, *bvector;

  GetDateTime(Value,_RealTime);
  sscanf(Value,"%i",&RealTid1);
  Tid=clock();

  if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
  else DebugNiveau=_DNormal;
      
  Print(_DNormal,"\n********************************************\n\n");
  Print(_DNormal,"Iterative Reconstruction program version 2.0\n");
  Print(_DNormal,"    Peter Toft and Jesper James Jensen\n");
  Print(_DNormal,"     Implemented in R by Joern Schulz\n");
  Print(_DNormal,"\n********************************************\n");


  ReadItArgs("InputData", *mode, UseFast, *RadonKernel, *IterationsType, Iterations, SaveIterations, *SaveIterationsName, LowestALevel, ConstrainMin, ConstrainMax, Alpha, Beta, Regularization, KernelFileSave, *KernelFileName, *RefFileName, ThetaSamples, ThetaMin, DeltaTheta, RhoSamples,  RhoMin, DeltaRho, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples, OverSamp, *DebugLevel);  
  
  // ==================================================================
  // Initialization of radon-image and the System-Matrix A
  
  if (itINI.IsFast)  AMatrix=GenerateAMatrix();
  
  bvector=DoubleToFloatVector(InImage, itINI.ThetaSamples, itINI.RhoSamples);

  // ==================================================================
  // Regularisation and initialisation of a Start-Image

  /* If regularisation is used, concatenate the bvector with zeroes */
  if (itINI.Regularization>0) 
    VectorCat(bvector,InitVector(AMatrix->M-bvector->N));

  if (*StartImageTrue) {
    Print(_DDetail, "A StartImage was specified. \n");
    xvector=DoubleToFloatVector(StartImage, itINI.XSamples, itINI.YSamples);
  }
  else {
    xvector=InitVector(itINI.XSamples*itINI.YSamples);
    if (itINI.IsFast) {
      /* The fast version, we can use the a-matrix */
      tempsuma=0; tempsumb=0;
      tempvector=SumRowSparseMatrix(AMatrix);
      for(n=0; n<tempvector->N; n++)
        tempsuma+=tempvector->value[n];
      for(n=0; n<bvector->N; n++)
        tempsumb+=bvector->value[n];
      mean=tempsumb/tempsuma;
      for(n=0; n<xvector->N; n++)
        xvector->value[n]=mean;
    } else {
      /* we have to estimate the mean value */
      tempsumb=0;
      for(n=0; n<bvector->N; n++)
        tempsumb+=bvector->value[n];
      mean=tempsumb/(itINI.ThetaSamples*itINI.RhoSamples*itINI.XSamples*itINI.DeltaX);
      for(n=0; n<xvector->N; n++)
        xvector->value[n]=mean;
    }    
    //Print(_DNormal,"Estimated mean of output image: %f \n",mean);
    Print(_DNormal,"Image-Initialization with: %f \n",mean);
  }

  // ==================================================================
  // Choice the corresponding Algorithm
  
  switch (itINI.Algorithm){           /* Main choice of algorithm */
    case _CG:   
      if (itINI.IsFast)
        MyImage=FAST_CG(AMatrix,xvector,bvector);
      else
        MyImage=SLOW_CG(xvector,bvector);                       
      break;
    case _EM:
      if (itINI.IsFast)
        MyImage=FAST_EM(AMatrix,xvector, bvector);                       
      else
        MyImage=SLOW_EM(xvector, bvector);                       
      break;
    case _ART:   
      if (itINI.IsFast)
        MyImage=FAST_ART(AMatrix,xvector,bvector);
      else
        MyImage=SLOW_ART(xvector,bvector);
      break;
  }
  
  ScaleImage(MyImage);
  Print(_DDetail,"\n");
  PrintStats(_DDetail, MyImage);

  // ==================================================================
  // Storage of reconstructed data and empty out the memory.
  
  ImageToFloat(OutImage, MyImage);
  FreeImage(MyImage);
  FreeVector(xvector);
  FreeVector(bvector);
  if (itINI.IsFast)
    FreeSparseMatrix(AMatrix);
  
  // ==================================================================
  // Calculation of the time-consuming.
  
  Tid=(clock()-Tid)/(float)CLOCKS_PER_SEC;
  GetDateTime(Value,_RealTime);
  sscanf(Value,"%i",&RealTid2);
  Print(_DNormal,"IT: Program was active for %.2f seconds\n",Tid) ;
  Print(_DNormal,"    World time elapsed %d seconds\n",
	(RealTid2-RealTid1));
  Print(_DNormal,"    Program used %.2f %% cpu time\n",
         Tid/((float)RealTid2-RealTid1)*100);
  Print(_DNormal,"\n");

}
