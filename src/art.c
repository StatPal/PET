/*****************************************************************
[HEADER]

This unit contains the programs associated with iterative
reconstruction with the ART algorithm.

Dec. 94, JJJ and PT
July 06, J.Schulz
*****************************************************************/

#include "it.h"
#include<Rmath.h>

/****************************************************************************
[NAME]
ARTUpdateAddVector

[SYNOPSIS]
void ARTUpdateAddVector(Vector *MyV, 
                        SparseMatrix *MySm, 
                        float weigth, 
                        int rownr)

[DESCRIPTION] 
This function preforms the following ART specific update operation
$$\mathbf b =\mathbf b - \alpha \mathbf a_n,$$  where $\mathbf b$={\tt MyV},
$\mathbf A$={\tt MySm}, $\alpha$={\tt weigth} and $n$={\tt rownr}.

[USAGE]
Do not use, used internally.

[REVISION]
Jan. 95, JJJ and PT
March 2007, J.Schulz 
      tempV[tempI[n]]=tempSM[n]*weigth; modified to tempV[tempI[n]]+=tempSM[n]*weigth;  
February 2019, J.Schulz
      random number generation by R instead of rand()
****************************************************************************/
void ARTUpdateAddVector(Vector *MyV, SparseMatrix *MySm, float weigth, int rownr)
{
  int n, *tempI, tempNm;
  float *tempSM, *tempV, *tempVn;

  tempI=MySm->index[rownr];
  tempSM=MySm->value[rownr];
  tempV=MyV->value;
  tempNm=MySm->Nm[rownr];

  if ((itINI.ConstrainMin<0) || (itINI.ConstrainMax<0)) 
    for(n=0; n<tempNm; n++) 
      tempV[tempI[n]]+=tempSM[n]*weigth; 
  else {
    for(n=0; n<tempNm; n++) {
      tempVn=&tempV[tempI[n]];      *tempVn+=tempSM[n]*weigth; 
      if (*tempVn<itINI.ConstrainMin) *tempVn=itINI.ConstrainMin;
      else 
	if (*tempVn>itINI.ConstrainMax) *tempVn=itINI.ConstrainMax; 
    }
  }
}


/****************************************************************************
[NAME]
ARTUpdateAddVector2

[SYNOPSIS]
void ARTUpdateAddVector2(Vector *MyV, Vector *MyAVector, float weigth)

[DESCRIPTION] 
This function preforms the following ART specific update operation
$$\mathbf b =\mathbf b - \alpha \mathbf a_n,$$  where $\mathbf b$={\tt MyV},
$\mathbf A$={\tt MySm}, $\alpha$={\tt weigth} and $n$={\tt rownr}. Does not 
use sparse matrix.

[USAGE]
Do not use, used internally.

[REVISION]
April 96 PT JJJ
****************************************************************************/
void ARTUpdateAddVector2(Vector *MyV, Vector *MyAVector, float weigth)
{
  int n;
  float *TempA, *TempV;

  TempV=MyV->value;
  TempA=MyAVector->value;
  for (n=0;n<MyV->N;n++)
    TempV[n]+=weigth*TempA[n];
  
  if ((itINI.ConstrainMin>=0) && (itINI.ConstrainMax>=0)) 
    ConstrainVector(MyV,itINI.ConstrainMin,itINI.ConstrainMax);
}


/****************************************************************************
[NAME]
FAST\_ART

[SYNOPSIS]
Image *FAST_ART(SparseMatrix *AMatrix, 
                Vector *xvector,
                Vector *bvector)

[DESCRIPTION]
This function will iterate towards a solution for the sparse system of 
equations \mb{b}=\mb{A x}, where \mb{b} is the sinogram (Radon domain)
and \mb{x} is the reconstructed image to be found. The function uses a 
fast version of ART (Algebraric Reconstruction Techniques) where a 
pre-calculated transformation matrix is used.

[USAGE]
{\tt Image=ART(TestMatrix, TestSinogram);}

Reconstructs the sinogram {\tt TestSinogram}, returns it as an image.

[REVISION]
Jan. 95, JJJ and PT
July 06, J.Schulz, Modification of ReadRefImage and the condition to
                   SaveIterions
****************************************************************************/
Image *FAST_ART(SparseMatrix *AMatrix, Vector *xvector, Vector *bvector)
{
  int ARows,ACols,currentrow,currentiteration,TotalIterations,AntPrint,UseRefImage;
  float lambda,brk;
  //float refxdev=0.0;
  float *tempXv,*tempBv;
  char DiffFileName[200];
  FILE *DiffFile=NULL;
  Vector *refxvector=NULL, *L2Norm;
  Image *Recon,*RefImage=NULL;

  ARows=AMatrix->M;
  ACols=AMatrix->N;

  Print(_DNormal,"Using ART (fast) to solve %i equations with %i unknowns\n",ARows,ACols);

  UseRefImage=(strlen(itINI.RefFileName)!=0);
  if (UseRefImage!=0) {
    RefImage=ReadRefImage(itINI.RefFileName);
    refxvector=ImageToVector(RefImage);
    FreeImage(RefImage);
    //refxdev=DeviationVector(refxvector);
    strcpy(DiffFileName,itINI.RefFileName);
    strcat(DiffFileName,".dif");
    DiffFile=fopen(DiffFileName,"wt");
    Print(_DNormal,"Logging differences in `%s' \n", DiffFileName);
  }

  tempXv=xvector->value;
  tempBv=bvector->value;
  //srand((int)clock());
  lambda=itINI.Alpha/itINI.Beta;

  TotalIterations=itINI.Iterations*ARows;
  AntPrint=(int)(TotalIterations/93);
  L2Norm=SumSqRowSparseMatrix(AMatrix);
  for (currentiteration=0;currentiteration<TotalIterations;currentiteration++)
  {
    if (currentiteration%ARows==0) lambda*=itINI.Beta;
    if (currentiteration%AntPrint==0) 
      Print(_DNormal,"Iterating %6.2f %% done\r",
	    (currentiteration+1)*100.0/TotalIterations); 
    if (itINI.IterationType==1)
      currentrow=currentiteration%ARows;
    else {
      //currentrow=(int)(ARows*(float)rand()/(RAND_MAX+1.0));
      GetRNGstate();
      currentrow = (int)(ARows*runif(0, 1));
      //currentrow = (int)(ARows*runif(0, RAND_MAX)/(RAND_MAX+1.0));
      PutRNGstate();
    }
    if (AMatrix->Nm[currentrow]>0) {
      brk=lambda*(tempBv[currentrow]-
		  MultSparseMatrixRowVector(AMatrix,xvector,currentrow))/
		    L2Norm->value[currentrow];
      ARTUpdateAddVector(xvector, AMatrix, brk, currentrow);      
    }
    if ((itINI.SaveIterations) && (currentiteration != 0) && (!(currentiteration%(ARows*(itINI.SaveIterations)))))
      SaveIteration(xvector,(int)(currentiteration/ARows),itINI.SaveIterationsName);
    
    if (UseRefImage==1)
      if (currentiteration%AntPrint==0)
	  fprintf(DiffFile,"%f %f\n",(double)currentiteration/ARows,
	 	    (double)L2NormVector(refxvector,xvector));
	/*fprintf(DiffFile,"%f %f\n",(double)currentiteration/ARows,
		(double)L2NormVector(refxvector,xvector,refxdev));*/ 
  }
  Print(_DNormal,"                                                  \r");
  Recon=VectorToImage(xvector,itINI.XSamples,itINI.YSamples);
  if (UseRefImage==1){
    Print(_DNormal,"L2 = %9.6f \n",L2NormVector(refxvector,xvector));
    //Print(_DNormal,"L2 = %9.6f \n",L2NormVector(refxvector,xvector,refxdev));
    FreeVector(refxvector);
    fclose(DiffFile);
    }
  
  RenameImage(Recon,"ReconstructedImage");
  Recon->DeltaX=itINI.DeltaX;
  Recon->DeltaY=itINI.DeltaY;
  Recon->Xmin=itINI.Xmin;
  Recon->Ymin=itINI.Ymin;
  FreeVector(L2Norm);

  return Recon;
}

/****************************************************************************
[NAME]
SLOW\_ART

[SYNOPSIS]
Image *SLOW_ART(Vector *xvector,
                Vector *bvector)

[DESCRIPTION]
This function will iterate towards a solution for the sparse system of 
equations \mb{b}=\mb{A x}, where \mb{b} is the sinogram (Radon domain)
and \mb{x} is the reconstructed image to be found. The function uses a 
slow version of ART (Algebraric Reconstruction Techniques) where the 
transformation matrix is calculated on the fly.

[USAGE]
{\tt Image=ART(TestMatrix, TestSinogram);}

Reconstructs the sinogram {\tt TestSinogram}, returns it as an image.

[REVISION]
Jan. 95, JJJ and PT
July 06, J.Schulz, Modification of ReadRefImage and the condition to
                   SaveIterions
****************************************************************************/
Image *SLOW_ART(Vector *xvector, Vector *bvector)
{
  int ARows,ACols,currentrow,currentiteration,TotalIterations,AntPrint,UseRefImage;
  float denom,lambda,brk;
  //refxdev=0.0;
  float *tempXv,*tempBv;
  char DiffFileName[200];
  FILE *DiffFile=NULL;
  Vector *refxvector=NULL,*AVector;
  Image *Recon,*RefImage=NULL;

  ARows=bvector->N;
  ACols=xvector->N;

  Print(_DNormal,"Using ART (slow) to solve %i equations with %i unknowns\n",
	ARows,ACols);

  UseRefImage=(strlen(itINI.RefFileName)!=0);
  if (UseRefImage!=0) {
    RefImage=ReadRefImage(itINI.RefFileName);
    refxvector=ImageToVector(RefImage);
    FreeImage(RefImage);
    //refxdev=DeviationVector(refxvector);
    strcpy(DiffFileName,itINI.RefFileName);
    strcat(DiffFileName,".dif");
    DiffFile=fopen(DiffFileName,"wt");
    Print(_DNormal,"Logging differences in `%s' \n", DiffFileName);
  }

  tempXv=xvector->value;
  tempBv=bvector->value;
  //srand((int)clock());
  lambda=itINI.Alpha/itINI.Beta;

  TotalIterations=itINI.Iterations*ARows;
  AntPrint=(int)(TotalIterations/93);

  InitArrays();

  for (currentiteration=0;currentiteration<TotalIterations;currentiteration++)
  { 
    if (currentiteration%ARows==0) lambda*=itINI.Beta;
    if (currentiteration%AntPrint==0)
      Print(_DNormal,"Iterating %6.2f %% done\r",
	    (currentiteration+1)*100.0/TotalIterations); 
    if (itINI.IterationType==1)
      currentrow=currentiteration%ARows; 
    else {
      //currentrow=(int)(ARows*(float)rand()/(RAND_MAX+1.0));
      GetRNGstate();
      currentrow = (int)(ARows*runif(0, 1));
      //currentrow = (int)(ARows*runif(0, RAND_MAX)/(RAND_MAX+1.0));
      PutRNGstate();
    }

    AVector=GenerateAMatrixRow(currentrow); 
    denom=MultVectorVector(AVector,AVector);
    if (fabs(denom)>1e-9)
    {
      brk=lambda*(tempBv[currentrow]-MultVectorVector(AVector,xvector))/denom;       
      ARTUpdateAddVector2(xvector, AVector, brk);
    }      
    FreeVector(AVector);

    if ( (itINI.SaveIterations) && (currentiteration != 0) && (!(currentiteration%(ARows*(itINI.SaveIterations)))) )
      SaveIteration(xvector,(int)(currentiteration/ARows),itINI.SaveIterationsName);
    
    if (UseRefImage==1)
      if (currentiteration%AntPrint==0)
	  fprintf(DiffFile,"%f %f\n",(double)currentiteration/ARows,
		    (double)L2NormVector(refxvector,xvector)); 
	  /*fprintf(DiffFile,"%f %f\n",(double)currentiteration/ARows,
		    (double)L2NormVector(refxvector,xvector,refxdev));*/
  }
  Print(_DNormal,"                                                  \r");
  Recon=VectorToImage(xvector,itINI.XSamples,itINI.YSamples);
  if (UseRefImage==1){
    Print(_DNormal,"L2 = %9.6f \n",L2NormVector(refxvector,xvector));
    //Print(_DNormal,"L2 = %9.6f \n",L2NormVector(refxvector,xvector,refxdev));
    FreeVector(refxvector);
    fclose(DiffFile);
  }
  
  RenameImage(Recon,"ReconstructedImage");
  Recon->DeltaX=itINI.DeltaX;
  Recon->DeltaY=itINI.DeltaY;
  Recon->Xmin=itINI.Xmin;
  Recon->Ymin=itINI.Ymin;
  
  return Recon;
}
