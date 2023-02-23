// WHy the hell the h files do not have header protectors!


#include "iradon_smoothed.h"
#include <R.h>
#include <Rdefines.h>


#define Inf 1e+140
#define SQ(x)   ((x)*(x))
#define MAX(a, b) ((a > b) ? a : b)
#define NEG(x) ((x < 0) ? 1 : 0)
#define SIGN(x) ((x) > 0 ? 1 : -NEG(x)) /* mimics the Splus sign function */


// a is a matrix of size m x n, vectorize it to matrix b of size mn
void vectorize(int m, int n, double **a, double *b){
    for (int j = 0; j < n; ++j){
        for (int i = 0; i < m; ++i){
            b[j * m + i] = a[i][j];
        }
    }
}

// Matricize b to a
void matricize(int m, int n, double **a, double *b){
    for (int j = 0; j < n; ++j){
        for (int i = 0; i < m; ++i){
            a[i][j] = b[j * m + i];
        }
    }
}




/* 
What is the extra in this function than BackProject??

In the previous function, the values like InvMyImage->Xmin etc was not present. 
Hence this function is necessary. 
However, if this function is directly put into R, we should use a wrapper.

So, basically it only backprojects, i.e., f -> K'f
*/
Image *BackProjection(Image *MyImage)
{
  int OldHeight,OldWidth,XSamples,YSamples;
  float Xmin,Ymin;
  Image *InvMyImage;

  // Print(_DNormal,"Backproject transforming: '%s'\n",MyImage->FileName);
  Print(_DNormal,"\nMyImage dimensions in Backprojection 1: M:%i N:%i\n",MyImage->M,MyImage->N);

  //M=N=(int)((MyImage->N-1)/(float)sqrt(2))+1;

  Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples);
  
  OldHeight=IniFile.XSamples;
  OldWidth =IniFile.YSamples;
  XSamples=1<<(int)(log(IniFile.XSamples)/log(2)+1);
  YSamples=1<<(int)(log(IniFile.YSamples)/log(2)+1);
  Xmin=IniFile.Xmin+((int)((OldHeight-XSamples-1)/2))*IniFile.DeltaX;
  Ymin=IniFile.Ymin+((int)((OldWidth-YSamples-1)/2))*IniFile.DeltaY;

  /* Allocate new image, and put transformation parameters in it*/
  InvMyImage=NewFloatImage("RecImage",XSamples,YSamples,_RealArray);  // Change
  InvMyImage->Xmin=Xmin;
  InvMyImage->Ymin=Ymin;
  InvMyImage->DeltaX=IniFile.DeltaX;
  InvMyImage->DeltaY=IniFile.DeltaY;

  BackProject(MyImage,InvMyImage);
  return(InvMyImage);
}






/***************************************************************************
 * Eigenvalues of K'K
 * Okay, here 
 * `InvMyImage->Signal[m];`
 * is not used. 
 * Basically the Signal should be multiplied by the eig.
***************************************************************************/
void filter_new(Image *MyImage, double **eig){

  int i,m,n,mm,nn;
  int CenterM,CenterN;
  float *Resx,*Resy,TempFloat,*TempPoint;

  Image *InvMyImage;
  InvMyImage = BackProjection(MyImage);


  /* Eigenvalues of the Filter(K'K) of the backprojected image */
  FFTImage(InvMyImage,_FFT);  // InvMyImage -> V' InvMyImage  // Is this necessary here??
  CenterM=InvMyImage->M/2;
  CenterN=InvMyImage->N/2;

  Resx=FloatVector(InvMyImage->M);
  Resy=FloatVector(InvMyImage->N);
  
  for(m=0;m<InvMyImage->M;m++) {
    mm=m;
    if(mm>=CenterM) mm=InvMyImage->M-m;
    Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
      (InvMyImage->DeltaX*InvMyImage->DeltaX);
    /*Print(_DDebug,"Resx %5.5f\n ",Resx[m]);*/
  }
  
  for(n=0;n<InvMyImage->N;n++) {
    nn=n;
    if(nn>=CenterN) nn=InvMyImage->N-n;
    Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
      (InvMyImage->DeltaY*InvMyImage->DeltaY);
    /*Print(_DDebug,"Resy %5.5f\n",Resy[n]);*/
  }
  

  // double **eig;
  // MAKE_2ARRAY(eig, (size_t)InvMyImage->M, (size_t)InvMyImage->N);
  Print(_DNormal,"InvMyImage dimensions: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 256 256

	FILE *ffile;
  ffile=fopen("eigens_inside_filter.dat","w");
  for(m=0;m<InvMyImage->M;m++) {
    TempFloat=Resx[m];
    // TempPoint=InvMyImage->Signal[m];

    // I think n and i are same in this loop
    for(n=0,i=0;n<InvMyImage->N;n++) {
      // eig[m][i++]=sqrt(TempFloat+Resy[n]);
      eig[m][i]=sqrt(TempFloat+Resy[n]);
      fprintf(ffile,"%f ", eig[m][i]);
      i++;
    }
  }
  fclose(ffile);

  FreeImage(InvMyImage);
  Free(Resx);
  Free(Resy);

}




/***************************************************************************
***************************************************************************/
void iradon_smoothed_C(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
{
  //Image *NewImage, *InvNewImage, *spectrum, *RefImage;
}





// void only_BackProject_C(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
// {
//   Image *NewImage;
  
//   if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
//   else DebugNiveau=_DNormal;
//   ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

//   // initialization of radon-image
//   NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
//   RDoubleToImage(NewImage, InImage, *M, *N );
//   InitImage(NewImage);
  



//   // Instead of Backfilter
//   int i,m,n,mm,nn,OldHeight,OldWidth;
//   int XSamples1,YSamples1,CenterM,CenterN;
//   float Xmin1,Ymin1,Res;
//   Image *InvMyImage;

//   // Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135  
//   OldHeight=IniFile.XSamples;
//   OldWidth =IniFile.YSamples;
//   XSamples1=1<<(int)(log(IniFile.XSamples)/log(2)+1);
//   YSamples1=1<<(int)(log(IniFile.YSamples)/log(2)+1);
//   Xmin1=IniFile.Xmin+((int)((OldHeight-XSamples1-1)/2))*IniFile.DeltaX;
//   Ymin1=IniFile.Ymin+((int)((OldWidth-YSamples1-1)/2))*IniFile.DeltaY;

//   /* Allocate new image, and put transformation parameters in it*/
//   InvMyImage=NewFloatImage("RecImage",XSamples1,YSamples1,_RealArray);  // Change
//   InvMyImage->Xmin=Xmin1;
//   InvMyImage->Ymin=Ymin1;
//   InvMyImage->DeltaX=IniFile.DeltaX;
//   InvMyImage->DeltaY=IniFile.DeltaY;

//   Image *NewImagecpy;
//   NewImagecpy = CopyImage(NewImage);  // NewImagecpy = NewImage;  // Both changes
  
//   Print(_DNormal,"\nOriginal NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
//   BackProject(NewImage,InvMyImage);
//   Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x157




//   double **eig;
//   int eigen_M = InvMyImage->M;
//   int eigen_N = InvMyImage->N;
//   MAKE_2ARRAY(eig, (size_t)eigen_M, (size_t)eigen_N);
//   printf("eig dim %d, %d\n", eigen_M, eigen_N);  // 64x128
//   printf("\nCase 1 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x135
//   filter_new(NewImagecpy, eig);
//   printf("\nCase 2 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x157
//   FREE_MATRIX(eig);



//   // Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
//   ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
//   Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
//   RealImage(InvMyImage);

//   //NormImage(InvMyImage,1.0,-MeanValue(InvMyImage));  // Change
//   //PrintStats(_DDetail,InvMyImage);
  

//   ScaleImage(InvMyImage);
//   PrintStats(_DDetail,InvMyImage);
//   ImageToFloat(OutImage, InvMyImage);
//   FreeImage(InvMyImage);
  
//   FreeImage(NewImage);
//   FreeImage(NewImagecpy);
// }



void BackProjection_C(double *InImage, double *OutImage, double *backfilter, double *eig_out, int *Xdim_modified, int *Ydim_modified, 
                      char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
{
  Image *NewImage;
  
  if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
  else DebugNiveau=_DNormal;
  ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

  // initialization of radon-image
  NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
  RDoubleToImage(NewImage, InImage, *M, *N );
  InitImage(NewImage);
  


  // Instead of Backfilter
  int i,m,n,OldHeight,OldWidth,XSamples1,YSamples1;
  float Xmin1,Ymin1;
  Image *InvMyImage, *NewImagecpy;

  // Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
  // Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples); //61x105
  
  OldHeight=IniFile.XSamples;
  OldWidth =IniFile.YSamples;
  XSamples1=1<<(int)(log(IniFile.XSamples)/log(2)+1);
  YSamples1=1<<(int)(log(IniFile.YSamples)/log(2)+1);
  Xmin1=IniFile.Xmin+((int)((OldHeight-XSamples1-1)/2))*IniFile.DeltaX;
  Ymin1=IniFile.Ymin+((int)((OldWidth-YSamples1-1)/2))*IniFile.DeltaY;

  *Xdim_modified = XSamples1;
  *Ydim_modified = YSamples1;

  /* Allocate new image, and put transformation parameters in it*/
  InvMyImage=NewFloatImage("RecImage",XSamples1,YSamples1,_RealArray);  // Change
  InvMyImage->Xmin=Xmin1;
  InvMyImage->Ymin=Ymin1;
  InvMyImage->DeltaX=IniFile.DeltaX;
  InvMyImage->DeltaY=IniFile.DeltaY;


  NewImagecpy = CopyImage(NewImage);  // NewImagecpy = NewImage;  // Both changes  
  Print(_DNormal,"\nOriginal NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
  BackProject(NewImage,InvMyImage);
  Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x157
  
  ImageToFloat(backfilter, InvMyImage);  // NEW output to R
  // Do it by hand




  double **eig;
  int eigen_M = InvMyImage->M, eigen_N = InvMyImage->N;
  MAKE_2ARRAY(eig, (size_t)eigen_M, (size_t)eigen_N);
  printf("eig dim %d, %d\n", eigen_M, eigen_N);  // 64x128
  // printf("\nCase 1 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x135
  filter_new(NewImagecpy, eig);
  // printf("\nCase 2 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x157
 
  
  /* Filter the backprojected image */
  FFTImage(InvMyImage,_FFT);

  for(m=0;m<InvMyImage->M;m++) {
    for(n=0,i=0;n<InvMyImage->N;n++) {
      InvMyImage->Signal[m][i++]*=eig[m][n];
      InvMyImage->Signal[m][i++]*=eig[m][n];
      // A small difference is due to being double instead of float
    }
  }
  vectorize(eigen_M, eigen_N, eig, eig_out);  // output to R
  FREE_MATRIX(eig);


  FFTImage(InvMyImage,_IFFT);
  // Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
  ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
  Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
  RealImage(InvMyImage);


  ScaleImage(InvMyImage);
  PrintStats(_DDetail,InvMyImage);
  ImageToFloat(OutImage, InvMyImage);
  FreeImage(InvMyImage);
  
  FreeImage(NewImage);
  FreeImage(NewImagecpy);

  Print(_DNormal,"return to R.          \n");
}








// Test 
void BackProjection_C_test(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
{
  Image *NewImage;
  
  if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
  else DebugNiveau=_DNormal;
  ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

  // initialization of radon-image
  NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
  RDoubleToImage(NewImage, InImage, *M, *N );
  InitImage(NewImage);
  







  // Instead of Backfilter
  int i,m,n,mm,nn,OldHeight,OldWidth;
  int XSamples1,YSamples1,CenterM,CenterN;
  float Xmin1,Ymin1,Res,*Resx,*Resy,TempFloat,*TempPoint;
  Image *InvMyImage;

  // Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135

  // M1=N1=(int)((NewImage->N-1)/(float)sqrt(2))+1;
  // Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples); //61x105
  
  OldHeight=IniFile.XSamples;
  OldWidth =IniFile.YSamples;
  XSamples1=1<<(int)(log(IniFile.XSamples)/log(2)+1);
  YSamples1=1<<(int)(log(IniFile.YSamples)/log(2)+1);
  Xmin1=IniFile.Xmin+((int)((OldHeight-XSamples1-1)/2))*IniFile.DeltaX;
  Ymin1=IniFile.Ymin+((int)((OldWidth-YSamples1-1)/2))*IniFile.DeltaY;


  /* Allocate new image, and put transformation parameters in it*/
  InvMyImage=NewFloatImage("RecImage",XSamples1,YSamples1,_RealArray);  // Change
  InvMyImage->Xmin=Xmin1;
  InvMyImage->Ymin=Ymin1;
  InvMyImage->DeltaX=IniFile.DeltaX;
  InvMyImage->DeltaY=IniFile.DeltaY;



  Image *NewImagecpy;
  NewImagecpy = CopyImage(NewImage);  // NewImagecpy = NewImage;  // Both changes
  
  Print(_DNormal,"\nOriginal NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
  BackProject(NewImage,InvMyImage);
  Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x157




  double **eig;
  int eigen_M = InvMyImage->M;
  int eigen_N = InvMyImage->N;
  MAKE_2ARRAY(eig, (size_t)eigen_M, (size_t)eigen_N);
  printf("eig dim %d, %d\n", eigen_M, eigen_N);  // 64x128
  printf("\nCase 1 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x135
  filter_new(NewImagecpy, eig);
  printf("\nCase 2 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x157
  // FREE_MATRIX(eig);


  
  /* Filter the backprojected image */
  FFTImage(InvMyImage,_FFT);

  // CenterM=InvMyImage->M/2;
  // CenterN=InvMyImage->N/2;
  // Resx=FloatVector(InvMyImage->M);
  // Resy=FloatVector(InvMyImage->N);
  
  // for(m=0;m<InvMyImage->M;m++) {
  //   mm=m;
  //   if(mm>=CenterM) mm=InvMyImage->M-m;
  //   Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
  //     (InvMyImage->DeltaX*InvMyImage->DeltaX);
  // }
  
  // for(n=0;n<InvMyImage->N;n++) {
  //   nn=n;
  //   if(nn>=CenterN) nn=InvMyImage->N-n;
  //   Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
  //     (InvMyImage->DeltaY*InvMyImage->DeltaY);
  // }
  
  // Print(_DNormal,"InvMyImage dimensions before multiplying by Eigen: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
  // for(m=0;m<InvMyImage->M;m++) {
  //   TempFloat=Resx[m];
  //   TempPoint=InvMyImage->Signal[m];
  //   for(n=0,i=0;n<InvMyImage->N;n++) {
  //     Res=sqrt(TempFloat+Resy[n]);
  //     /*Print(_DDebug,"Res %5.5f\n",Res);*/
  //     TempPoint[i++]*=Res;
  //     TempPoint[i++]*=Res;


  //     // if(fabs(eig[m][n] - Res) > 1e-7) printf("eig %f, Res %f \t", eig[m][n], Res);
  //   }
  // }


  // The difference is due to being double instead of float
  for(m=0;m<InvMyImage->M;m++) {
    // TempPoint=InvMyImage->Signal[m];
    for(n=0,i=0;n<InvMyImage->N;n++) {
      InvMyImage->Signal[m][i++]*=eig[m][n];
      InvMyImage->Signal[m][i++]*=eig[m][n];
    }
  }
  FREE_MATRIX(eig);



  FFTImage(InvMyImage,_IFFT);
  // Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
  ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
  Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
  RealImage(InvMyImage);

  // //NormImage(InvMyImage,1.0,-MeanValue(InvMyImage));  // Change
  // //PrintStats(_DDetail,InvMyImage);
  // Free(Resx);
  // Free(Resy);
  





  ScaleImage(InvMyImage);
  PrintStats(_DDetail,InvMyImage);
  ImageToFloat(OutImage, InvMyImage);
  FreeImage(InvMyImage);
  
  FreeImage(NewImage);
  FreeImage(NewImagecpy);  
  Print(_DNormal,"return to R.          \n");

}









