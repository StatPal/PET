/*****************************************************************
[HEADER]

This unit contains the programs associated with iterative
reconstruction with the CG-LS algorithm.

Dec. 94, JJJ and PT
July 06, J.Schulz
*****************************************************************/

#include "it.h"

/****************************************************************************
[NAME]
CGUpdateAddVector

[SYNOPSIS]
void CGUpdateAddVector(Vector *MyV1, 
                       Vector *MyV2, 
                       Vector *MyV3, 
                       float weigth)

[DESCRIPTION] 

This function preforms the following EM specific update operation
$$\mathbf{a=b - \alpha c},$$ where $\mathbf a$={\tt MyV3}, $\mathbf
b$={\tt MyV1}, $\mathbf c$={\tt MyV2} and $\alpha$={\tt weigth}.

[USAGE]
Do not use, used internally.

[REVISION]
Jan. 95, JJJ and PT
****************************************************************************/
void CGUpdateAddVector(Vector *MyV1, Vector *MyV2, Vector *MyV3, float weigth)
{
  int n, tempN;
  float *tempV1, *tempV2, *tempV3;

  if ((MyV1->N!=MyV2->N) || (MyV2->N!=MyV3->N))
    Error("Incompatible sizes encountered (CGUpdateAddVector)");

  tempV1=MyV1->value;
  tempV2=MyV2->value;
  tempV3=MyV3->value;
  tempN=MyV1->N;
  
  for(n=0; n<tempN; n++) 
    tempV3[n]=tempV1[n]+tempV2[n]*weigth; 
}


/****************************************************************************
[NAME]
FAST\_CG

[SYNOPSIS]
Image *FAST_CG(SparseMatrix *AMatrix, 
               Vector *xvector,
               Vector *bvector,)


[DESCRIPTION]
This function will iterate towards a solution for the sparse system of 
equations \mb{b}=\mb{A x}, where \mb{b} is the sinogram (Radon domain)
and \mb{x} is the reconstructed image to be found. The function uses CG-LS
(Conjugate Gradients Least Squares).

[USAGE]
{\tt Image=FAST\_CG(TestMatrix, TestSinogram);}

Reconstructs the sinogram {\tt TestSinogram}, returns it as an image.

[REVISION]
Jan. 95, JJJ and PT\\
April 7 PT Bug with regularization. Calc of ARows ACols changed.\\
Nov 1 96 PT Bug in update formulas.
July  06, J.Schulz, Modification of ReadRefImage and the condition to
                    SaveIterions
Dec  06, J.Schulz, if-Conditions for ConstrainVector was absent
****************************************************************************/
Image *FAST_CG(SparseMatrix *AMatrix,Vector *xvector, Vector *bvector)
{
  int ARows,ACols,currentiteration,UseRefImage;
  float alpha,beta,rnorm,newrnorm,qnorm;
  //float refxdev=0.0;
  char DiffFileName[200];
  FILE *DiffFile=NULL;
  Vector *refxvector=NULL, *pvector, *rvector, *svector, *qvector,*tempvec;
  Image *Recon,*RefImage;

  ARows=AMatrix->M;
  ACols=AMatrix->N;
  Print(_DNormal,"Using CG to solve %i equations with %i unknowns\n",
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

  tempvec=InitVector(ARows);
  pvector=InitVector(ACols);
  rvector=InitVector(ACols);
  qvector=InitVector(ARows);
  svector=InitVector(ARows);

  Print(_DNormal,"Initializing transformation parameters.\n");

  MultSparseMatrixVector(AMatrix,xvector,tempvec);     /*  temp = A * x0  */
  CGUpdateAddVector(bvector,tempvec,svector,-1.0);     /*  s0 = b - temp  */
  MultSparseTMatrixVector(AMatrix, svector, rvector);  /*  r0 = AT * s0   */  
  memcpy(pvector->value,rvector->value,sizeof(float)*rvector->N);/*  p0 = r0 */
  MultSparseMatrixVector(AMatrix, pvector, qvector);    /*  q0 = A * p0  */  
  rnorm=MultVectorVector(rvector,rvector);              /*  rnorm = r0 * r0 */
 
  Print(_DNormal,"CG running\n");
  for (currentiteration=0;currentiteration<itINI.Iterations;currentiteration++) {
    Print(_DNormal,"Iterating %6.2f %% done\r",(currentiteration+1)*100.0/itINI.Iterations); 
    qnorm=MultVectorVector(qvector,qvector);          /* qnorm = q0 * q0  */
    alpha=rnorm/qnorm;                                /* alpha = rnorm/qnorm */
    CGUpdateAddVector(xvector,pvector,xvector,alpha); /* xk+1 = xk + alpha pk*/
    
    /* Constrain solution  */
    if ((itINI.ConstrainMin>=0) && (itINI.ConstrainMax>=0))
        ConstrainVector(xvector,itINI.ConstrainMin,itINI.ConstrainMax);
    
    // The succeeding line was activated and the line after rk+1=AT*sk was deactivated. That is the order of the fast case in the C-routines of P.Toft. But at the moment I'm not sure, what is the correct order. In the slow case of the C-code of Peter Toft and in his PhD thesis the line sk+1=sk-alpha qk follow the line rk+1=AT*sk. If you know the correct order, please send me an email and I will remove the warning in case of mode==CG. If possible, please tell me a reference/literature with the correct algorithm.
    CGUpdateAddVector(svector,qvector,svector,-alpha);/* sk+1 = sk-alpha qk */
    MultSparseTMatrixVector(AMatrix,svector,rvector); /* rk+1 = AT * sk  */
    //CGUpdateAddVector(svector,qvector,svector,-alpha);/* sk+1 = sk-alpha qk */
    newrnorm=MultVectorVector(rvector,rvector);       /* newrnorm=rk+1 * rk+1*/
    beta=newrnorm/rnorm;                              /* beta=newrnorm/rnorm */
    rnorm=newrnorm;                                   /* rnorm=newrnorm      */
    CGUpdateAddVector(rvector,pvector,pvector,beta);  /* pk+1=rk+1 + beta pk */
    MultSparseMatrixVector(AMatrix,pvector,qvector);  /* qk+1 = A * pk+1     */

    if ( itINI.SaveIterations && (currentiteration!=0) && (!(currentiteration%(itINI.SaveIterations))) )
      SaveIteration(xvector,currentiteration,itINI.SaveIterationsName);

    if (UseRefImage==1)
	fprintf(DiffFile,"%i %e\n",currentiteration,
	        L2NormVector(refxvector,xvector)); 
	/*fprintf(DiffFile,"%i %e\n",currentiteration,
                L2NormVector(refxvector,xvector,refxdev));*/
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


/****************************************************************************
[NAME]
SLOW\_CG

[SYNOPSIS]
Image *SLOW_CG(Vector *xvector,
               Vector *bvector,)


[DESCRIPTION]
This function will iterate towards a solution for the sparse system of 
equations \mb{b}=\mb{A x}, where \mb{b} is the sinogram (Radon domain)
and \mb{x} is the reconstructed image to be found. The function uses CG-LS
(Conjugate Gradients Least Squares).

[USAGE]
{\tt Image=CG(TestMatrix, TestSinogram);}

Reconstructs the sinogram {\tt TestSinogram}, returns it as an image.

[REVISION]
Jan. 95, JJJ and PT
July 06, J.Schulz, Modification of ReadRefImage and the condition to
                   SaveIterions
Dec  06, J.Schulz, if-Conditions for ConstrainVector was absent
March 07 J.Schulz, Modification of the order of rk+1 = AT * sk and sk+1 = sk-alpha qk
****************************************************************************/
Image *SLOW_CG(Vector *xvector, Vector *bvector)
{
  int n,m,ARows,ACols,currentiteration,UseRefImage;
  float alpha,beta,rnorm,newrnorm,qnorm;
  //float refxdev=0.0;
  char DiffFileName[200];
  FILE *DiffFile=NULL;
  Vector *refxvector=NULL, *pvector, *rvector, *svector, *qvector,*tempvec;
  Vector *AVector, *ATVector;
  Image *Recon,*RefImage;

  InitArrays();

  ARows=itINI.ThetaSamples*itINI.RhoSamples;
  ACols=itINI.XSamples*itINI.YSamples;
  Print(_DNormal,"Using CG to solve %i equations with %i unknowns\n",
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

  tempvec=InitVector(ARows);
  pvector=InitVector(ACols);
  rvector=InitVector(ACols);
  qvector=InitVector(ARows);
  svector=InitVector(ARows);

  Print(_DNormal,"Initializing (transforming).\n");
  for (m=0;m<ARows;m++) {              /*  temp = A * x0       */
    AVector=GenerateAMatrixRow(m);
    tempvec->value[m]=MultVectorVector(AVector,xvector);
    FreeVector(AVector);
  }
  CGUpdateAddVector(bvector,tempvec,svector,-1.0);     /*  s0 = b - temp   */
  for (n=0;n<ACols;n++) {                               /*  r0 = AT * s0   */  
    ATVector=GenerateAMatrixColumn(n);
    rvector->value[n]=MultVectorVector(ATVector,svector);
    FreeVector(ATVector);
  }
  memcpy(pvector->value,rvector->value,sizeof(float)*rvector->N);/*  p0 = r0 */
  for (m=0;m<ARows;m++) {                                 /*  q0 = A * p0 */
    AVector=GenerateAMatrixRow(m);
    qvector->value[m]=MultVectorVector(AVector,pvector);
    FreeVector(AVector);
  }
  rnorm=MultVectorVector(rvector,rvector);              /*  rnorm = r0 * r0 */
 
  Print(_DNormal,"CG running\n");
  for (currentiteration=0;currentiteration<itINI.Iterations;currentiteration++) {
    Print(_DNormal,"Iterating %6.2f %% done\r",
          (currentiteration+1)*100.0/itINI.Iterations); 
    qnorm=MultVectorVector(qvector,qvector);        /* qnorm = q0 * q0   */
    alpha=rnorm/qnorm;                              /* alpha = rnorm/qnorm */
    CGUpdateAddVector(xvector,pvector,xvector,alpha);/* xk+1 = xk + alpha pk */
    
    /* Constrain solution  */
    if ((itINI.ConstrainMin>=0) && (itINI.ConstrainMax>=0))
        ConstrainVector(xvector,itINI.ConstrainMin,itINI.ConstrainMax);
    
    // The succeeding line was activated and the line after rk+1=AT*sk was deactivated. That is the order of the fast case in the C-routines of P.Toft. But at the moment I'm not sure, what is the correct order. In the slow case of the C-code of Peter Toft and in his PhD thesis the line sk+1=sk-alpha qk follow the line rk+1=AT*sk. If you know the correct order, please send me an email and I will remove the warning in case of mode==CG. If possible, please tell me a reference/literature with the correct algorithm.
    CGUpdateAddVector(svector,qvector,svector,-alpha);/* sk+1 = sk - alpha qk*/
    for (n=0;n<ACols;n++) {                       /*  rk+1 = AT * sk      */  
      ATVector=GenerateAMatrixColumn(n);
      rvector->value[n]=MultVectorVector(ATVector,svector);
      FreeVector(ATVector);
    }
    //CGUpdateAddVector(svector,qvector,svector,-alpha);/* sk+1 = sk - alpha qk*/
    newrnorm=MultVectorVector(rvector,rvector);   /* newrnorm = rk+1 * rk+1 */
    beta=newrnorm/rnorm;                              /* beta=newrnorm/rnorm */
    rnorm=newrnorm;                                   /* rnorm=newrnorm      */
    CGUpdateAddVector(rvector,pvector,pvector,beta);/* pk+1 = rk+1 + beta pk */
    for (m=0;m<ARows;m++) {                           /* qk+1 = A * pk+1     */
      AVector=GenerateAMatrixRow(m);
      qvector->value[m]=MultVectorVector(AVector,pvector);
      FreeVector(AVector);
    }
    if ( (itINI.SaveIterations) && (currentiteration!=0) && (!(currentiteration%(itINI.SaveIterations))) )
      SaveIteration(xvector,currentiteration,itINI.SaveIterationsName);

    if (UseRefImage==1)
	fprintf(DiffFile,"%i %e\n",currentiteration,
                L2NormVector(refxvector,xvector)); 
	/*fprintf(DiffFile,"%i %e\n",currentiteration,
                L2NormVector(refxvector,xvector,refxdev));*/
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

