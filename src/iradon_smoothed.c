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






// void Forward(Image *MyImage,Image *MySino)
// {
// 	int t,r,m,n,mmin,mmax,nmin,nmax;
// 	int T,R,M,N;
// 	float sum,x_min,rhooffset,Delta_theta,Delta_rho,costheta,sintheta;
// 	float rho_min,theta,alpha,beta,nfloat,mfloat;
// 	float Delta_x,betap,eps;
// 	float **signal,**sino;

// 	eps=1e-4;
// 	T=MySino->M;
// 	R=MySino->N;
// 	M=N=MyImage->M;

// 	rho_min=MySino->Ymin;
// 	Delta_theta=M_PI/T;
// 	Delta_rho=MySino->DeltaY;

// 	x_min=MyImage->Xmin;
// 	Delta_x=MyImage->DeltaX;
// 	signal=MyImage->Signal;
// 	sino=MySino->Signal;

// 	printf("Sino X=%i DX=%f Xmin=%f\n",MySino->M,MySino->DeltaX,MySino->Xmin);
// 	printf("Sino Y=%i DY=%f Ymin=%f\n",MySino->N,MySino->DeltaY,MySino->Ymin);
// 	printf("Image X=%i DX=%f Xmin=%f\n",MyImage->M,MyImage->DeltaX,MyImage->Xmin);
// 	printf("Image Y=%i DY=%f Ymin=%f\n",MyImage->N,MyImage->DeltaY,MyImage->Ymin);

// 	for(t=0; t<T;t++)
// 	{
// 		Print(_DNormal, "Calculating: %4d of %d \r",t,T-1);
// 		theta=t*Delta_theta;
// 		sintheta=sin(theta);
// 		costheta=cos(theta);
// 		rhooffset=x_min*(costheta+sintheta);
// 		if (sintheta>sqrt(0.5))
// 		{
// 			alpha=-costheta/sintheta;
// 			for(r=0; r<R; r++ ) 
// 			{
// 				beta=(r*Delta_rho+rho_min-rhooffset)/(Delta_x*sintheta);
// 				betap=beta+0.5;
// 				sum=0.0;
// 				if (alpha>1e-6)
// 				{
// 					mmin=(int)ceil(-(beta+0.5-eps)/alpha);
// 					mmax=1+(int)floor((N-0.5-beta-eps)/alpha);
// 				}
// 				else
// 					if (alpha<-1e-6)
// 					{
// 						mmin=(int)ceil((N-0.5-beta-eps)/alpha);
// 						mmax=1+(int)floor(-(beta+0.5-eps)/alpha);
// 					}
// 					else
// 					{
// 						mmin=0;
// 						mmax=M;
// 					}
// 				if (mmin<0) mmin=0;
// 				if (mmax>M) mmax=M;
// 				nfloat=betap+mmin*alpha;
// 				for (m=mmin;m<mmax;m++)
// 				{
// 					sum+=signal[m][(int)nfloat];
// 					nfloat+=alpha;
// 				}
// 				sino[t][r]=sum*Delta_x/fabs(sintheta);
// 			}
// 		}
// 		else
// 		{
// 			alpha=-sintheta/costheta;
// 			for(r=0; r<R; r++ ) 
// 			{
// 				beta=(r*Delta_rho+rho_min-rhooffset)/(Delta_x*costheta);
// 				betap=beta+0.5;
// 				sum=0.0;
// 				if (alpha>1e-6)
// 				{
// 					nmin=(int)ceil(-(beta+0.5-eps)/alpha);
// 					nmax=1+(int)floor((M-0.5-beta-eps)/alpha);
// 				}
// 				else 
// 					if (alpha<-1e-6)
// 					{
// 						nmin=(int)ceil((M-0.5-beta-eps)/alpha);
// 						nmax=1+(int)floor(-(beta+0.5-eps)/alpha);
// 					}
// 					else
// 					{
// 						nmin=0;
// 						nmax=N;
// 					}
// 				if (nmin<0) nmin=0;
// 				if (nmax>N) nmax=N;
// 				mfloat=betap+nmin*alpha;
// 				for (n=nmin;n<nmax;n++)
// 				{
// 					sum+=signal[(int)mfloat][n];
// 					mfloat+=alpha;
// 				}
// 				sino[t][r]=sum*Delta_x/fabs(costheta);
// 			}
// 		}
// 	}
// 	Print(_DNormal,"\n Finished         \n");
// }





// void Forward_C_orig_dim(double *InImage, double *OutImage, 
//                       char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, 
//                       double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
// {
//   Image *NewImage;
//   ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

//   // initialization of radon-image
//   NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
//   RDoubleToImage(NewImage, InImage, *M, *N );
//   InitImage(NewImage);
  

//   // Instead of Backfilter
//   int i,m,n,OldHeight,OldWidth,XSamples1,YSamples1;
//   // Image *InvMyImage, *NewImagecpy;
  
//   OldHeight=IniFile.XSamples;
//   OldWidth =IniFile.YSamples;
//   XSamples1=IniFile.XSamples;
//   YSamples1=IniFile.YSamples;


//   Image *MySino;

//   /* Allocate new image, and put transformation parameters in it*/
//   MySino=NewFloatImage("RecImage",XSamples1,YSamples1,_RealArray);  // Change
//   MySino->Xmin=IniFile.Xmin;
//   MySino->Ymin=IniFile.Ymin;
//   MySino->DeltaX=IniFile.DeltaX;
//   MySino->DeltaY=IniFile.DeltaY;


//   // NewImagecpy = CopyImage(NewImage);  // NewImagecpy = NewImage;  // Both changes  
//   // Print(_DNormal,"\nOriginal NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
//   // BackProject(NewImage,InvMyImage);
//   // Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x157
  
//   Forward(NewImage, MySino);

//   ImageToFloat(OutImage, MySino);  // NEW output to R
//   // Do it by hand

//   FreeImage(MySino);

//   Print(_DNormal,"return to R.          \n");
// }


















// void Forward(Image *MyImage,Image *MySino)
void Forward_C_orig_dim(double *sino_out, double *signal_in, double *setpar)
{
	int t,r,m,n,mmin,mmax,nmin,nmax;
	int T,R,M,N;
	float sum,x_min,rhooffset,Delta_theta,Delta_rho,costheta,sintheta;
	float rho_min,theta,alpha,beta,nfloat,mfloat;
	float Delta_x,betap,eps;
	float **signal,**sino;

	eps=1e-4;
	// T=MySino->M;
	// R=MySino->N;
	// M=N=MyImage->M;

	// rho_min=MySino->Ymin;
	// Delta_theta=M_PI/T;
	// Delta_rho=MySino->DeltaY;

	// x_min=MyImage->Xmin;
	// Delta_x=MyImage->DeltaX;

	// printf("Sino X=%i DX=%f Xmin=%f\n",MySino->M,MySino->DeltaX,MySino->Xmin);
	// printf("Sino Y=%i DY=%f Ymin=%f\n",MySino->N,MySino->DeltaY,MySino->Ymin);
	// printf("Image X=%i DX=%f Xmin=%f\n",MyImage->M,MyImage->DeltaX,MyImage->Xmin);
	// printf("Image Y=%i DY=%f Ymin=%f\n",MyImage->N,MyImage->DeltaY,MyImage->Ymin);


  // printf("Sino X=%i DX=%f Xmin=%f\n", T, MySino->DeltaX, MySino->Xmin);
	// printf("Sino Y=%i DY=%f Ymin=%f\n", R, Delta_rho, rho_min);
	// printf("Image X=%i DX=%f Xmin=%f\n", M, Delta_x,MyImage->Xmin);
	// printf("Image Y=%i DY=%f Ymin=%f\n", N,MyImage->DeltaY,MyImage->Ymin);

	
  M            = setpar[0];
  N = M;
  x_min        = setpar[1];
  Delta_x      = setpar[2];
  T            = setpar[3];
  Delta_theta  = setpar[4];
  R            = setpar[5];
  rho_min      = setpar[6];
  Delta_rho    = setpar[7];

  float theta_min;
  theta_min    = setpar[8];     // This was not used in Peter Toft's code. 
  
  // printf("\n\n\n\n\n\n\n");
  // printf("Sino X=%i DX=N Xmin=%f(N)\n", T, theta_min);
	// printf("Sino Y=%i DY=%f Ymin=%f\n", R, Delta_rho, rho_min);
	// printf("Image X=%i DX=%f Xmin=N\n", M, Delta_x);
	// printf("Image Y=%i DY=N Ymin=N\n", N);





	// signal=MyImage->Signal;
	// sino=MySino->Signal;

  // printf("Debug 0\n");

  MAKE_2ARRAY(signal, M, N);
  MAKE_2ARRAY(sino, T, R);


  // printf("Debug 1\n");
  
  // matricize(M, N, signal, signal_in);
  for (int j = 0; j < N; ++j){
    for (int i = 0; i < M; ++i){
        signal[i][j] = signal_in[j * M + i];
    }
  }

  // printf("Debug 2\n");


	for(t=0; t<T;t++)
	{
		Print(_DNormal, "Calculating: %4d of %d \r",t,T-1);

    // printf("\n\nCalculating t: %4d of %d \n", t, T-1);

		theta=t*Delta_theta;
		sintheta=sin(theta);
		costheta=cos(theta);
		rhooffset=x_min*(costheta+sintheta);
		if (sintheta>sqrt(0.5))
		{
			alpha=-costheta/sintheta;
			for(r=0; r<R; r++ ) 
			{

        printf("Calculating r: %4d of %d \r",r, R-1);

				beta=(r*Delta_rho+rho_min-rhooffset)/(Delta_x*sintheta);
				betap=beta+0.5;
				sum=0.0;
				if (alpha>1e-6)
				{
					mmin=(int)ceil(-(beta+0.5-eps)/alpha);
					mmax=1+(int)floor((N-0.5-beta-eps)/alpha);
				}
				else
					if (alpha<-1e-6)
					{
						mmin=(int)ceil((N-0.5-beta-eps)/alpha);
						mmax=1+(int)floor(-(beta+0.5-eps)/alpha);
					}
					else
					{
						mmin=0;
						mmax=M;
					}
				if (mmin<0) mmin=0;
				if (mmax>M) mmax=M;
				nfloat=betap+mmin*alpha;
				for (m=mmin;m<mmax;m++)
				{
					sum+=signal[m][(int)nfloat];
					nfloat+=alpha;
				}
				sino[t][r]=sum*Delta_x/fabs(sintheta);
			}
		}
		else
		{
			alpha=-sintheta/costheta;
      // printf("\n Toft t: %d, alpha: %f, beta: %f\n", t, alpha, (0*Delta_rho+rho_min-rhooffset)/(Delta_x*costheta));  // Same

			for(r=0; r<R; r++ ) 
			{
        // printf("Calculating r (else): %4d of %d \r",r, R-1);

				beta=(r*Delta_rho+rho_min-rhooffset)/(Delta_x*costheta);
				betap=beta+0.5;	// Not in that code
				sum=0.0;
				if (alpha>1e-6)
				{
					nmin=(int)ceil(-(beta+0.5-eps)/alpha);	// Almost same
					nmax=1+(int)floor((M-0.5-beta-eps)/alpha);	// Almost same
				}
				else 
					if (alpha<-1e-6)
					{
						nmin=(int)ceil((M-0.5-beta-eps)/alpha);
						nmax=1+(int)floor(-(beta+0.5-eps)/alpha);
					}
					else
					{
						nmin=0;
						nmax=N;
					}
				if (nmin<0) nmin=0;
				if (nmax>N) nmax=N;


        // printf("\n What!! t = %d, r = %d, beta: %f, N: %d, nmin, nmax: %d, %d", t, r, beta, N, nmin, nmax);
        

				mfloat=betap+nmin*alpha;  // Negative value
        // printf("\n Debug inside 1");

				for (n=nmin;n<nmax;n++)
				{

          // if(t==0) printf("\n Debug inside 1.5: %f, %d, %f", betap, nmin, alpha);
          // printf("\n Debug inside 1.5: %d, %d", (int)mfloat, n);
          // printf("\n Debug inside 1.5: %f, %d, %f", betap, nmin, alpha);  // -27.000000, 0, 0
					sum+=signal[(int)mfloat][n];
					mfloat+=alpha;
				}
        // printf("\n Debug inside 2");
				sino[t][r]=sum*Delta_x/fabs(costheta);
			}
		}
	}

  // printf("Debug 3\n");

  // vectorize(T, R, sino, sino_out);
  for (int j = 0; j < R; ++j){
      for (int i = 0; i < T; ++i){
          sino_out[j * T + i] = sino[i][j];
      }
  }

  // printf("Debug 4\n");


  FREE_MATRIX(signal);
  FREE_MATRIX(sino);

	Print(_DNormal,"\n Finished         \n");
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
  // Print(_DNormal,"\nMyImage dimensions in Backprojection 1: M:%i N:%i\n",MyImage->M,MyImage->N);

  //M=N=(int)((MyImage->N-1)/(float)sqrt(2))+1;

  // Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples);
  
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



Image *BackProjection_same_dim(Image *MyImage)
{
  int OldHeight,OldWidth,XSamples,YSamples;
  Image *InvMyImage;
  
  OldHeight=IniFile.XSamples;
  OldWidth =IniFile.YSamples;
  XSamples=IniFile.XSamples;
  YSamples=IniFile.YSamples;

  /* Allocate new image, and put transformation parameters in it*/
  InvMyImage=NewFloatImage("RecImage",XSamples,YSamples,_RealArray);  // Change
  InvMyImage->Xmin=IniFile.Xmin;
  InvMyImage->Ymin=IniFile.Ymin;
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
  // Print(_DNormal,"InvMyImage dimensions: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 256 256

	// FILE *ffile;
  // ffile=fopen("eigens_inside_filter.dat","w");
  for(m=0;m<InvMyImage->M;m++) {
    TempFloat=Resx[m];
    // // TempPoint=InvMyImage->Signal[m];

    // // I think n and i are same in this loop
    for(n=0,i=0;n<InvMyImage->N;n++) {
      // eig[m][i++]=sqrt(TempFloat+Resy[n]);
      eig[m][i]=sqrt(TempFloat+Resy[n]);
      // fprintf(ffile,"%f ", eig[m][i]);
      i++;
    }
  }
  // fclose(ffile);

  FreeImage(InvMyImage);
  Free(Resx);
  Free(Resy);

}




void filter_new_same_dim(Image *MyImage, double **eig){

  int i,m,n,mm,nn;
  int CenterM,CenterN;
  float *Resx,*Resy,TempFloat,*TempPoint;

  Image *InvMyImage;
  InvMyImage = BackProjection_same_dim(MyImage);

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
  }
  
  for(n=0;n<InvMyImage->N;n++) {
    nn=n;
    if(nn>=CenterN) nn=InvMyImage->N-n;
    Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
      (InvMyImage->DeltaY*InvMyImage->DeltaY);
  }

  for(m=0;m<InvMyImage->M;m++) {
    TempFloat=Resx[m];
    // // I think n and i are same in this loop
    for(n=0,i=0;n<InvMyImage->N;n++) {
      eig[m][i]=sqrt(TempFloat+Resy[n]);
      i++;
    }
  }
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





void BackProject_C_orig_dim(double *InImage, double *OutImage, double *backfilter, double *eig_out, int *Xdim_modified, int *Ydim_modified, 
                      char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
{
  Image *NewImage;
  ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

  // initialization of radon-image
  NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
  RDoubleToImage(NewImage, InImage, *M, *N );
  InitImage(NewImage);
  

  // Instead of Backfilter
  int i,m,n,OldHeight,OldWidth,XSamples1,YSamples1;
  Image *InvMyImage, *NewImagecpy;
  
  OldHeight=IniFile.XSamples;
  OldWidth =IniFile.YSamples;
  XSamples1=IniFile.XSamples;
  YSamples1=IniFile.YSamples;

  /* Allocate new image, and put transformation parameters in it*/
  InvMyImage=NewFloatImage("RecImage",XSamples1,YSamples1,_RealArray);  // Change
  InvMyImage->Xmin=IniFile.Xmin;
  InvMyImage->Ymin=IniFile.Ymin;
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
  printf("\nCase 1 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x135
  filter_new_same_dim(NewImagecpy, eig);
  printf("\nCase 2 NewImagecpy-> M %d, NewImagecpy-> N %d\n", NewImagecpy->M, NewImagecpy->N);  // 320x157
 
  
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
  Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
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



void BackProjection_C_shrinked(double *InImage, double *OutImage, double *backfilter, double *eig_out, int *Xdim_modified, int *Ydim_modified, 
                      char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
{
  Image *NewImage;
  ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

  // initialization of radon-image
  NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
  RDoubleToImage(NewImage, InImage, *M, *N );
  InitImage(NewImage);
  
  // Instead of Backfilter
  int i,m,n,OldHeight,OldWidth,XSamples1,YSamples1;
  float Xmin1,Ymin1;
  Image *InvMyImage, *NewImagecpy, *InvMyImagecpy;

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
  // Print(_DNormal,"\nOriginal NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
  BackProject(NewImage,InvMyImage);
  // Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x157
  
  NewImagecpy = CopyImage(NewImage);
  InvMyImagecpy = CopyImage(InvMyImage);
  
  // Print(_DNormal,"Original InvMyImage dimensions before shrinkimage: M:%i N:%i\n",InvMyImagecpy->M,InvMyImagecpy->N);  // 65x105
  ShrinkImage(InvMyImagecpy,OldHeight,OldWidth,_MiddleMiddle);
  // Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImagecpy->M,InvMyImagecpy->N);  // 65x105

  // Oh!! I see where the error is, after this Shrinking, the next dim are changing.
  // Copy this and store and then free
  ImageToFloat(backfilter, InvMyImagecpy);  // NEW output to R
  FreeImage(InvMyImagecpy);


  double **eig;
  int eigen_M = InvMyImage->M, eigen_N = InvMyImage->N;
  MAKE_2ARRAY(eig, (size_t)eigen_M, (size_t)eigen_N);
  // printf("eig dim %d, %d\n", eigen_M, eigen_N);  // 64x128
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



  // There is a multiple thing going on
  double *eig_tmp;
  MAKE_1ARRAY(eig_tmp, eigen_M*eigen_N);
  vectorize(eigen_M, eigen_N, eig, eig_tmp);  // vectorize eig to eig_tmp

  // To use Shrinkimage on eigen, we can turn it into a image first
  Image *eigen_Image;
  eigen_Image=NewFloatImage(IniFile.InFile, eigen_M, eigen_M, _RealArray);// initialization of radon-image
  RDoubleToImage(eigen_Image, eig_tmp, eigen_M, eigen_N );  // Put eig_tmp to eigen_Image
  InitImage(eigen_Image);

  // Print(_DNormal,"Original eigen_Image dimensions before eigen Shrink: M:%i N:%i\n",eigen_Image->M,eigen_Image->N);  // 64x64
  ShrinkImage(eigen_Image, OldHeight, OldWidth, _MiddleMiddle);
  // Print(_DNormal,"Original eigen_Image dimensions after  eigen Shrink: M:%i N:%i\n",eigen_Image->M,eigen_Image->N);  // 65x105

  ImageToFloat(eig_out, eigen_Image);  // NEW output to R
  FREE_MATRIX(eig);
  FreeImage(eigen_Image);
  // NOT CHECKED - CHECK






  FFTImage(InvMyImage,_IFFT);
  // Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
  ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
  // Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
  RealImage(InvMyImage);


  ScaleImage(InvMyImage);
  PrintStats(_DDetail,InvMyImage);
  ImageToFloat(OutImage, InvMyImage);
  FreeImage(InvMyImage);
  
  FreeImage(NewImage);
  FreeImage(NewImagecpy);

  // Print(_DNormal,"return to R.          \n");
}



void BackProjection_C(double *InImage, double *OutImage, double *backfilter, double *eig_out, int *Xdim_modified, int *Ydim_modified, 
                      char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
{
  Image *NewImage;
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
  // Print(_DNormal,"\nOriginal NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
  BackProject(NewImage,InvMyImage);
  // Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x157
  
  ImageToFloat(backfilter, InvMyImage);  // NEW output to R
  // Do it by hand




  double **eig;
  int eigen_M = InvMyImage->M, eigen_N = InvMyImage->N;
  MAKE_2ARRAY(eig, (size_t)eigen_M, (size_t)eigen_N);
  // printf("eig dim %d, %d\n", eigen_M, eigen_N);  // 64x128
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
  // Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
  RealImage(InvMyImage);


  ScaleImage(InvMyImage);
  PrintStats(_DDetail,InvMyImage);
  ImageToFloat(OutImage, InvMyImage);
  FreeImage(InvMyImage);
  
  FreeImage(NewImage);
  FreeImage(NewImagecpy);

  // Print(_DNormal,"return to R.          \n");
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









