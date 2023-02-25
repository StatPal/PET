// // WHy the hell the h files do not have header protectors!


// // #include <stdio.h>
// // #include <string.h>
// // #include <stdlib.h>
// // #include <math.h>
// // #include <R.h>
// // #include "imgtools.h"
// // #include "misc.h"
// // #include "calc.h"   // ??
// // #include "eval.h"

// #include "iradon_smoothed.h"
// #include <R.h>
// #include <Rdefines.h>


// #define Inf 1e+140
// #define SQ(x)   ((x)*(x))
// #define MAX(a, b) ((a > b) ? a : b)
// #define NEG(x) ((x < 0) ? 1 : 0)
// #define SIGN(x) ((x) > 0 ? 1 : -NEG(x)) /* mimics the Splus sign function */





// // a is a matrix of size m x n, vectorize it to matrix b of size mn
// void vectorize(int m, int n, double **a, double *b){
//     for (int j = 0; j < n; ++j){
//         for (int i = 0; i < m; ++i){
//             b[j * m + i] = a[i][j];
//         }
//     }
// }

// // Matricize b to a
// void matricize(int m, int n, double **a, double *b){
//     for (int j = 0; j < n; ++j){
//         for (int i = 0; i < m; ++i){
//             a[i][j] = b[j * m + i];
//         }
//     }
// }





// void BackProject_check(Image *Sinogram, Image *InvMyImage)
// {
//   int m,n,N,t,XSamples,YSamples,RhoInt;
//   Image *XCosTable, *YSinTable;
//   float Xmin,Ymin,DeltaX,DeltaY,Rho,DeltaRho,RhoMax,TempSin,TempCos,sum;
//   float *datam,*datan,*datat;

//   Print(_DNormal,"Backprojecting sinogram...\n");
//   PrintStats(_DDetail,Sinogram);
//   //PrintStats(_DDetail,InvMyImage);

//   Xmin=InvMyImage->Xmin;
//   Ymin=InvMyImage->Ymin;
//   DeltaX=InvMyImage->DeltaX; 
//   DeltaY=InvMyImage->DeltaY;  
//   XSamples=InvMyImage->M;    
//   YSamples=InvMyImage->N;

//   Print(_DDetail,
//         "Transformed Image: min=(%.2f,%.2f), dx=%.2f dy=%.2f M=%d N=%d\n",
// 	Xmin,Ymin,DeltaX,DeltaY,XSamples,YSamples);

//   /* Check if rho exceeded bounds in integration, if so stretch sinogram */
//   RhoMax = sqrt(sq(max((fabs(Xmin)),(fabs(Xmin+(XSamples-1)*DeltaX))))+
// 		sq(max((fabs(Ymin)),(fabs(Ymin+(YSamples-1)*DeltaY)))))+1;
  
//   Print(_DNormal," * MyImage dimensions in Backprojec_check 2.5\n: M:%i N:%i\n",Sinogram->M,Sinogram->N);
//   if (RhoMax>(-Sinogram->Ymin)) {
//     Print(_DDetail,"Required rho (%.2f) is bigger than maximum rho in sinogram (%.2f) \n"
// 	  ,RhoMax,-Sinogram->Ymin);
//     N = (int)(RhoMax/Sinogram->DeltaY)+1;
//     StretchImage(Sinogram,Sinogram->M,N*2+1,_MiddleMiddle);  // Here the image is stretched
//     //PrintStats(_DDetail,Sinogram);
//   }
//   Print(_DNormal," * MyImage dimensions in Backprojec_check 3\n: M:%i N:%i\n",Sinogram->M,Sinogram->N);


//   /* Allocate tabels and initialize */
//   YSinTable=NewFloatImage("YSin", YSamples, Sinogram->M,_RealArray);  
//   XCosTable=NewFloatImage("XCos", XSamples, Sinogram->M,_RealArray);  
//   Print(_DDetail,"Initializing tables... \n");
//   for(t=0; t<Sinogram->M; t++) {
//     TempCos=cos(PI*t/Sinogram->M)/Sinogram->DeltaY; 
//     TempSin=sin(PI*t/Sinogram->M)/Sinogram->DeltaY; 
//     for(m=0; m<XSamples; m++)
//       XCosTable->Signal[m][t]=((float)m*DeltaX+Xmin)*TempCos;
//     for(n=0; n<YSamples; n++)
//       YSinTable->Signal[n][t]=((float)n*DeltaY+Ymin)*TempSin-Sinogram->Ymin/Sinogram->DeltaY;
//   }

//   /* using linear interpolation to reduce noise in output image */
//   Print(_DDetail,"Integrating... \n");
//   for (m=0; m<XSamples; m++) {
//     //printf("(%.2f pct. done) \r", (float)(m+1)/XSamples*100);
//     if (m%5==0) Print(_DNormal,"(%.2f pct. done) \r", (float)(m+1)/XSamples*100);  
//     datam=XCosTable->Signal[m];
//     for (n=0; n<YSamples; n++) {
//       datan=YSinTable->Signal[n];
//       sum=0.0;
//       for (t=0; t<Sinogram->M; t++) {
//         RhoInt=(int)(Rho=(datam[t]+datan[t])); 
//         /* if (Rho>=Sinogram->N || Rho<0) printf("* %d %d %d %f \n", m,n,t,Rho); */
//         datat=&Sinogram->Signal[t][RhoInt];
//         //sum+=datat[0]*(1-(DeltaRho=(Rho-RhoInt)))+datat[1]*DeltaRho;
// 	  DeltaRho=(Rho-RhoInt);
// 	  sum+=datat[0]*(1-DeltaRho)+datat[1]*DeltaRho; 
//       }
//     InvMyImage->Signal[m][n]=sum*Sinogram->DeltaX;
//     }
//   }

//   FreeImage(XCosTable);
//   FreeImage(YSinTable);
// }




// /* 
// What is the extra in this function than BackProject??

// In the previous function, the values like InvMyImage->Xmin etc was not present. 
// Hence this function is necessary. 
// However, if this function is directly put into R, we should use a wrapper.

// So, basically it only backprojects, i.e., f -> K'f
// */
// Image *BackProjection(Image *MyImage)
// {
//   int i,m,n,mm,nn,M,N,OldHeight,OldWidth;
//   int XSamples,YSamples,CenterM,CenterN;
//   float Xmin,Ymin,Res,*Resx,*Resy,TempFloat,*TempPoint;
//   Image *InvMyImage;

//   Print(_DNormal,"Backproject transforming: '%s'\n",MyImage->FileName);
//   Print(_DNormal,"\n MyImage dimensions in Backprojection 1: M:%i N:%i\n",MyImage->M,MyImage->N);

//   M=N=(int)((MyImage->N-1)/(float)sqrt(2))+1;

//   Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples);
  
//   OldHeight=IniFile.XSamples;
//   OldWidth =IniFile.YSamples;
//   XSamples=1<<(int)(log(IniFile.XSamples)/log(2)+1);
//   YSamples=1<<(int)(log(IniFile.YSamples)/log(2)+1);
//   Xmin=IniFile.Xmin+((int)((OldHeight-XSamples-1)/2))*IniFile.DeltaX;
//   Ymin=IniFile.Ymin+((int)((OldWidth-YSamples-1)/2))*IniFile.DeltaY;

//   /* Allocate new image, and put transformation parameters in it*/
//   InvMyImage=NewFloatImage("RecImage",XSamples,YSamples,_RealArray);  // Change
//   InvMyImage->Xmin=Xmin;
//   InvMyImage->Ymin=Ymin;
//   InvMyImage->DeltaX=IniFile.DeltaX;
//   InvMyImage->DeltaY=IniFile.DeltaY;

//   BackProject_check(MyImage,InvMyImage);
//   return(InvMyImage);
// }






// /***************************************************************************
//  * Eigenvalues of K'K
// ***************************************************************************/
// void filter_new(Image *MyImage, double **eig){

//   int i,m,n,mm,nn,M,N,OldHeight,OldWidth;
//   int XSamples,YSamples,CenterM,CenterN;
//   float Xmin,Ymin,Res,*Resx,*Resy,TempFloat,*TempPoint;

//   Image *InvMyImage;
//   InvMyImage = BackProjection(MyImage);


//   /* Eigenvalues of the Filter(K'K) of the backprojected image */
//   FFTImage(InvMyImage,_FFT);  // InvMyImage -> V' InvMyImage  // Is this necessary here??
//   CenterM=InvMyImage->M/2;
//   CenterN=InvMyImage->N/2;

//   Resx=FloatVector(InvMyImage->M);
//   Resy=FloatVector(InvMyImage->N);
  
//   for(m=0;m<InvMyImage->M;m++) {
//     mm=m;
//     if(mm>=CenterM) mm=InvMyImage->M-m;
//     Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
//       (InvMyImage->DeltaX*InvMyImage->DeltaX);
//     /*Print(_DDebug,"Resx %5.5f\n ",Resx[m]);*/
//   }
  
//   for(n=0;n<InvMyImage->N;n++) {
//     nn=n;
//     if(nn>=CenterN) nn=InvMyImage->N-n;
//     Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
//       (InvMyImage->DeltaY*InvMyImage->DeltaY);
//     /*Print(_DDebug,"Resy %5.5f\n",Resy[n]);*/
//   }
  

//   // double **eig;
//   // MAKE_2ARRAY(eig, (size_t)InvMyImage->M, (size_t)InvMyImage->N);
//   Print(_DNormal,"InvMyImage dimensions: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 256 256

// 	FILE *ffile=fopen("eigens.dat","w");
//   for(m=0;m<InvMyImage->M;m++) {
//     TempFloat=Resx[m];
//     TempPoint=InvMyImage->Signal[m];
//     for(n=0,i=0;n<InvMyImage->N;n++) {
//       // eig[m][i++]=sqrt(TempFloat+Resy[n]);
//       eig[m][i]=sqrt(TempFloat+Resy[n]);
//       fprintf(ffile,"%f ", eig[m][i]);
//       i++;
//     }
//   }
//   fclose(ffile);


// }


// #include "Refft.h"

// double Uty_PET(Image *MyImage, double **z1)
// {

//   int i,m,n,mm,nn,M,N,OldHeight,OldWidth;
//   int XSamples,YSamples,CenterM,CenterN;
//   float Xmin,Ymin,Res,*Resx,*Resy,TempFloat,*TempPoint;

//   // y'y
//   double sum1 = 0.;
//   int j;
//   Print(_DNormal,"MyImage dimensions Uty_PET: M:%i N:%i\n",MyImage->M,MyImage->N);  // 256 256

//   for(i=0; i<MyImage->M; i++) {
// 		for(j=0; j<MyImage->N; j++) sum1+=SQ(MyImage->Signal[i][j]);
// 	}

//   FILE *ffile;
//   ffile=fopen("y.dat","w");   // This is different
//   for(m=0;m<MyImage->M;m++) {
//     for(n=0;n<MyImage->N;n++) {
//       fprintf(ffile,"%f ", MyImage->Signal[m][n]);
//     }
//   }
//   fclose(ffile);


  
//   Image *InvMyImage;
//   InvMyImage = BackProjection(MyImage);  // InvMyImage = K'y
//   Print(_DNormal,"InvMyImage dimensions Uty_PET: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 256 256
  

//   ffile=fopen("Kty.dat","w");
//   for(m=0;m<InvMyImage->M;m++) {
//     for(n=0;n<InvMyImage->N;n++) {
//       fprintf(ffile,"%f ", InvMyImage->Signal[m][n]);
//     }
//   }
//   fclose(ffile);




//   double **eig;
//   MAKE_2ARRAY(eig, (size_t)InvMyImage->M, (size_t)InvMyImage->N);
//   filter_new(MyImage, eig);  // eig is supposed to be eigenvalues of  K'K - but dim does not match
//   Print(_DNormal,"Eigens done\n\n");
//   printf("\nEigenvalues: n: %d eig_00: %f, eig_01: %f, eig_10: %f, eig_11: %f\n", InvMyImage->M, eig[0][0], eig[0][1], eig[1][0], eig[1][1]);
//   printf("\nMyImage-> M %d, MyImage-> N %d\n", MyImage->M, MyImage->N);
  
//   // OldHeight=IniFile.XSamples;
//   // OldWidth =IniFile.YSamples;
//   // ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle);    // CHECK THIS - Maybe BUG
//   // RealImage(InvMyImage);
//   // Print(_DNormal,"InvMyImage dimensions Uty_PET after shrink: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128 128



//   // Missed FFT - BUG SUbrata
//   // FFTImage(InvMyImage,_FFT);  // InvMyImage -> V' InvMyImage = V' K'y  // Check which FFT is to be used
//   // No, do the Maitra FFT here

//   double **VtKty;
//   MAKE_MATRIX(VtKty, InvMyImage->M, InvMyImage->N); 
//   // Copy InvMyImage to VtKty  
//   for(i=0; i<InvMyImage->M; i++) {
//     for (j=0; j<InvMyImage->N; j++){
//       VtKty[i][j] = InvMyImage->Signal[i][j];
//     }
//   }
//   Print(_DNormal,"Kty done\n\n");
//   Refftf2(InvMyImage->M, InvMyImage->N, VtKty);
//   Print(_DNormal,"VtKty done\n\n");






//   // THIS AREA IS TO BE CHECKED
//   // z1 is InvMyImage->Signal = D^{-1/2} V' K' y;
//   // Okay, z1 is scalar or vector?? y seems to be of dim n, not nxn!!

//   // DEBUG STARTS
//   // memcpy(z1, InvMyImage->Signal, sizeof(float)*MyImage->M*MyImage->N);
//   for(i=0; i<InvMyImage->M; i++) {
//     for (j=0; j<InvMyImage->N; j++){
//       // z1[i][j] = InvMyImage->Signal[i][j];
//       // Should have rearrangements???

//       z1[i][j] = VtKty[i][j] * eig[i][j];
//     }
//   }
//   double sum2 = 0.;
//   for(i=0;i<n;i++) {
// 		for(j=0;j<n;j++) sum2+=SQ(z1[i][j]);
// 	}


//   // FFTImage(InvMyImage,_IFFT);
//   // ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
//   // RealImage(InvMyImage);
//   // NormImage(InvMyImage,1.0,-MeanValue(InvMyImage));  // Change
//   //PrintStats(_DDetail,InvMyImage);
//   Free(Resx);
//   Free(Resy);
//   FREE_MATRIX(eig);
//   FREE_MATRIX(VtKty);
  
//   // return(InvMyImage);
//   return (sum1 - sum2);
// }









// typedef struct gcvslsdata {
// 	int npix;
//   int mpix;
// 	int nangle;
// 	int ndist;
// 	double **Z1;
// 	double Z2;
// } GCVSLSData;


// double SLSPRESS_2d(int ma, int md, int m, int n, double **U1tysq, double sumU2tysquared, double fwhm)
// {
// 	/*This calculates the GCV-estimate of the Predicted Residual Sums of Squares 
// 	  (PRESS) for the Smoothed Least Squares Estimate with Gaussian smoothing.*/
// 	int i, j;
// 	double trfwhm=0., sum=0., **eigfwhm;

// 	MAKE_MATRIX(eigfwhm, m, n);
// 	fftGauss_2d(m, n, fwhm, eigfwhm);	// eigfwhm = Omega_h 

// 	for(i=0;i<m;i++) {
// 		for(j=0;j<n;j++) {
// 			trfwhm+=eigfwhm[i][j]/(md*ma-n*n); 				// This scaling is done for c(h)
// 			sum+=SQ((1.-eigfwhm[i][j])*U1tysq[i][j]);		// Z1 part
// 		}
// 	}
// 	FREE_MATRIX(eigfwhm);
// 	sum+=sumU2tysquared*(1+2.*trfwhm+trfwhm*trfwhm);		// Z2 part
// 	return sum;
// }

// double GCVPRESS_SLS_2d(double X, void *Data)
// {
// 	if (X <=0) return Inf;
// 	else {
// 		GCVSLSData *tmp=Data;
// 		return SLSPRESS_2d(tmp->nangle, tmp->ndist, tmp->mpix, tmp->npix, tmp->Z1, tmp->Z2, X);
// 	}
// }


// void Optimize_SLS_2d(Image *MyImage, double *FWHM, int m, int n)  // I think m and n can be calculated!!
// {
// 	double fnl, fmin;
// 	const double abstol = 1e-12, lower = 1e-6, upper = 20;

// 	GCVSLSData *mydata;
// 	mydata = malloc( sizeof *mydata );
// 	mydata->mpix = m;
//   mydata->npix = n;
// 	mydata->nangle = MyImage->M;
// 	mydata->ndist = MyImage->N;
//   // MAKE_MATRIX(mydata->Z1, MyImage->M, MyImage->N);
//   // Print(_DNormal,"mydata->Z1 dimensions when created: M:%i N:%i\n",MyImage->M,MyImage->N);

//   MAKE_MATRIX(mydata->Z1, 2*m, 2*n);
//   Print(_DNormal,"mydata->Z1 dimensions when created: M:%i N:%i\n", m, n);

//   // Print(_DNormal,"Debug 1\n");
//   mydata->Z2 = Uty_PET(MyImage, mydata->Z1);
//   // Print(_DNormal,"Debug 2\n");

//   // Print(_DNormal,"Z2:%f\n", mydata->Z2);
//   // Print(_DNormal,"Debug 3\n");

// 	printf("lower val: %f\n", GCVPRESS_SLS_2d(lower, mydata));
// 	printf("upper val: %f\n", GCVPRESS_SLS_2d(upper, mydata));
// 	printf("1 val: %f\n\n", GCVPRESS_SLS_2d(1.0, mydata));

// 	fnl = Brent_fmin_new(lower, upper, GCVPRESS_SLS_2d, mydata, &fmin, abstol);
//   printf("min = %g at %g\n", fmin, fnl);
	
// 	(*FWHM) = fnl;
// 	FREE_MATRIX(mydata->Z1);
// 	free(mydata);
// }





// /***************************************************************************
// [NAME]
// Main Subroutine

// [DESCRIPTION] 

// This is the main subroutine function for the direct reconstruction
// methods. The routine have to call from R through iradon. 
// See above for all possibilities. 

// It reads the specified parameters getting from R
// and executes one of the following functions:

// \begin{itemize}
// \item CentralSlice.
// \item FilteredBack.
// \item BackFiltering.
// \item Convert.
// \item Trace.
// \item ImageInfo.
// \item Test.
// \end{itemize}

// [REVISION]
// Jan 2023, Subrata Pal, 
// ***************************************************************************/
// void iradon_smoothed_C(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
// {
//   Image *NewImage, *InvNewImage, *spectrum, *RefImage;
//   Print(_DNormal,"M:%i N:%i\n", *M, *N);  // 181, 183

  
//   if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
//   else DebugNiveau=_DNormal;
   
// //   Print(_DNormal,"\n");
// //   Print(_DNormal,"Start of inverse Radontransformation\n");
// //   Print(_DNormal,"------------------------------------\n");
// //   Print(_DNormal,"          iradon (ver 2.0)          \n");
// //   Print(_DNormal,"      Made by Jesper J. Jensen      \n");
// //   Print(_DNormal,"   Peter Philipsen and Peter Toft   \n");
// //   Print(_DNormal,"  Implemented in R by Joern Schulz  \n");
// //   Print(_DNormal,"------------------------------------\n");

//   ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 
      
//   // ==================================================================
//   // initialization of radon-image
//   NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
//   RDoubleToImage(NewImage, InImage, *M, *N );
//   InitImage(NewImage);
//   Print(_DNormal,"Initial Sinogram dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 181, 369


//   FILE *ffile;
//   ffile=fopen("y_init.dat","w");
//   int m1, n1;
//   for(m1=0;m1<NewImage->M;m1++) {
//     for(n1=0;n1<NewImage->N;n1++) {
//       fprintf(ffile,"%f ", NewImage->Signal[m1][n1]);
//     }
//   }
//   fclose(ffile);


//   //Performes an reconstruction using Filtered Backprojection.
//   // Print(_DNormal,"Using New Backprojection\n");
//   // InvNewImage=BackFilter(NewImage);   // y -> K'y
//   // No, this is the old function - NOT THIS
//   InvNewImage=BackProjection(NewImage);   // y -> K'y
//   // This changes NewImage too


//   Print(_DNormal,"NewImage dimensions after BackProjection: M:%i N:%i\n",NewImage->M,NewImage->N);  // 181, 369
//   ffile=fopen("y_init_after_backfilter.dat","w");  // NewImage Changed
//   for(m1=0;m1<NewImage->M;m1++) {
//     for(n1=0;n1<NewImage->N;n1++) {
//       fprintf(ffile,"%f ", NewImage->Signal[m1][n1]);
//     }
//   }
//   fclose(ffile);

  



//   // Print(_DNormal,"Backprojected(??????) image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples);  // 128 128
//   // These are OldHeight and NewHeight respectively
//   Print(_DNormal,"Backprojected(??????) image dim.: M:%i N:%i\n", InvNewImage->M, InvNewImage->N);  // 128 128



  

//   int i,j;
//   int m = InvNewImage->M;
//   int n = InvNewImage->N;

//   double FWHM;  // double piu=PI*(mma-1.0)/mma;
//   // Optimize_SLS_2d(NewImage, &FWHM, m, n);





//   double **fil, **f, fwhm;

//   MAKE_MATRIX(fil, m, n);
//   MAKE_MATRIX(f, m, n);

  
//   // if (fwhm == 0) {
//   //   for(i=0; i<m; i++) {
//   //     for (j=0; j<n; j++) f[i][j] = InvNewImage->Signal[i][j];
//   //   }
//   // }
//   // else {
//   //   Gausfilter_2d(m, n, fwhm, fil); // create Gaussian filter into fil 
//   //   double **g;
//   //   MAKE_MATRIX(g, m, n);
//   //   for(i=0; i<m; i++) {
//   //     for (j=0; j<n; j++) g[i][j] = InvNewImage->Signal[i][j];
//   //   }
//   //   convolve_2d(m, n, fil, g, f); /* create S_h(K'K)^{-1}K'y into f*/
//   // }
//   vectorize(m, n, f, OutImage);
//   FREE_MATRIX(fil);
//   FREE_MATRIX(f);

//   // ScaleImage(InvNewImage);
//   PrintStats(_DDetail,InvNewImage);
//   ImageToFloat(OutImage, InvNewImage);
//   FreeImage(InvNewImage);



//   FreeImage(NewImage);
//   Print(_DNormal,"return to R.          \n");

// }

// // void BackProjection_C(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
// // {
// //   Image *MyImage, *InvNewImage, *spectrum, *RefImage;

// //   if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
// //   else DebugNiveau=_DNormal;   
// //   ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 
      
// //   // ==================================================================
// //   // initialization of radon-image
// //   MyImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
// //   RDoubleToImage(MyImage, InImage, *M, *N );
// //   InitImage(MyImage);
// //   Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",MyImage->M,MyImage->N); // ThetaSamples, RhoSamples


// //   // // Old fn
// //   int i,m,n,mm,nn,M1,N1,OldHeight,OldWidth;
// //   int XSamples1,YSamples1,CenterM,CenterN;
// //   float Xmin1,Ymin1,Res,*Resx,*Resy,TempFloat,*TempPoint;
// //   Image *InvMyImage;

// //   Print(_DNormal,"Backproject transforming: '%s'\n",MyImage->FileName);
// //   Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",MyImage->M,MyImage->N);

// //   M1=N1=(int)((MyImage->N-1)/(float)sqrt(2))+1;
// //   Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples);
  
// //   OldHeight=IniFile.XSamples;
// //   OldWidth =IniFile.YSamples;
// //   XSamples1=1<<(int)(log(IniFile.XSamples)/log(2)+1);
// //   YSamples1=1<<(int)(log(IniFile.YSamples)/log(2)+1);
// //   Xmin1=IniFile.Xmin+((int)((OldHeight-XSamples1-1)/2))*IniFile.DeltaX;
// //   Ymin1=IniFile.Ymin+((int)((OldWidth-YSamples1-1)/2))*IniFile.DeltaY;

// //   /* Allocate new image, and put transformation parameters in it*/
// //   InvMyImage=NewFloatImage("RecImage",XSamples1,YSamples1,_RealArray);  // Change
// //   InvMyImage->Xmin=Xmin1;
// //   InvMyImage->Ymin=Ymin1;
// //   InvMyImage->DeltaX=IniFile.DeltaX;
// //   InvMyImage->DeltaY=IniFile.DeltaY;
// //   BackProject(MyImage,InvMyImage);

  
// //   FFTImage(InvMyImage,_FFT);
// //   CenterM=InvMyImage->M/2;
// //   CenterN=InvMyImage->N/2;
// //   Resx=FloatVector(InvMyImage->M);
// //   Resy=FloatVector(InvMyImage->N);
  
// //   for(m=0;m<InvMyImage->M;m++) {
// //     mm=m;
// //     if(mm>=CenterM) mm=InvMyImage->M-m;
// //     Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
// //       (InvMyImage->DeltaX*InvMyImage->DeltaX);
// //   }
  
// //   for(n=0;n<InvMyImage->N;n++) {
// //     nn=n;
// //     if(nn>=CenterN) nn=InvMyImage->N-n;
// //     Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
// //       (InvMyImage->DeltaY*InvMyImage->DeltaY);
// //   }
// //   for(m=0;m<InvMyImage->M;m++) {
// //     TempFloat=Resx[m];
// //     TempPoint=InvMyImage->Signal[m];
// //     for(n=0,i=0;n<InvMyImage->N;n++) {
// //       Res=sqrt(TempFloat+Resy[n]);
// //       TempPoint[i++]*=Res;
// //       TempPoint[i++]*=Res;
// //     }
// //   }

// //   FFTImage(InvMyImage,_IFFT);
// //   ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
// //   RealImage(InvMyImage);

// //   Free(Resx);
// //   Free(Resy);
  




// //   // vectorize(InvMyImage->M, InvMyImage->N, InvMyImage->Signal, OutImage);



// //   // WHAT??
// //   // Backfilter gives back something
// //   // ScaleImage(InvNewImage);
// //   ScaleImage(InvMyImage);
// //   PrintStats(_DDetail,InvNewImage);
// //   ImageToFloat(OutImage, InvNewImage);
// //   FreeImage(InvNewImage);
// // }





// // TEST function
// void BF_C(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
// {
//   Image *NewImage, *InvNewImage, *spectrum, *RefImage;
  
//   if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
//   else DebugNiveau=_DNormal;
//   ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

//   // initialization of radon-image
//   NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
//   RDoubleToImage(NewImage, InImage, *M, *N );
//   InitImage(NewImage);
  







//   // Instead of Backfilter
//   int i,m,n,mm,nn,M1,N1,OldHeight,OldWidth;
//   int XSamples1,YSamples1,CenterM,CenterN;
//   float Xmin1,Ymin1,Res,*Resx,*Resy,TempFloat,*TempPoint;
//   Image *InvMyImage;

//   Print(_DNormal,"Filter after Backproject transforming: '%s'\n",NewImage->FileName);
//   Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135

//   M1=N1=(int)((NewImage->N-1)/(float)sqrt(2))+1;
//   Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples); //65x105
  
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

//   Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
//   BackProject(NewImage,InvMyImage);
//   Print(_DNormal,"Original InvMyImage dimensions after Backproject: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128





  
//   /* Filter the backprojected image */
//   FFTImage(InvMyImage,_FFT);
//   Print(_DNormal,"Original InvMyImage dimensions after FFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128

//   CenterM=InvMyImage->M/2;
//   CenterN=InvMyImage->N/2;
//   Resx=FloatVector(InvMyImage->M);
//   Resy=FloatVector(InvMyImage->N);
  
//   for(m=0;m<InvMyImage->M;m++) {
//     mm=m;
//     if(mm>=CenterM) mm=InvMyImage->M-m;
//     Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
//       (InvMyImage->DeltaX*InvMyImage->DeltaX);
//   }
  
//   for(n=0;n<InvMyImage->N;n++) {
//     nn=n;
//     if(nn>=CenterN) nn=InvMyImage->N-n;
//     Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
//       (InvMyImage->DeltaY*InvMyImage->DeltaY);
//   }
  
//   Print(_DNormal,"Original InvMyImage dimensions before filter: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128

//   for(m=0;m<InvMyImage->M;m++) {
//     TempFloat=Resx[m];
//     TempPoint=InvMyImage->Signal[m];
//     for(n=0,i=0;n<InvMyImage->N;n++) {
//       Res=sqrt(TempFloat+Resy[n]);
//       /*Print(_DDebug,"Res %5.5f\n",Res);*/
//       TempPoint[i++]*=Res;
//       TempPoint[i++]*=Res;
//     }
//   }

//   FFTImage(InvMyImage,_IFFT);
//   Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
//   ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
//   Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
//   RealImage(InvMyImage);
//   Print(_DNormal,"Original InvMyImage dimensions after RealImage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105

//   //NormImage(InvMyImage,1.0,-MeanValue(InvMyImage));  // Change
//   //PrintStats(_DDetail,InvMyImage);
//   Free(Resx);
//   Free(Resy);
  









//   Print(_DNormal,"Using Filtering after Backprojection\n");
//   // InvNewImage=BackFilter(NewImage);
//   // ScaleImage(InvNewImage);
//   // PrintStats(_DDetail,InvNewImage);
//   // ImageToFloat(OutImage, InvNewImage);
//   // FreeImage(InvNewImage);

//   ScaleImage(InvMyImage);
//   PrintStats(_DDetail,InvMyImage);
//   ImageToFloat(OutImage, InvMyImage);
//   FreeImage(InvMyImage);
  

  
//   FreeImage(NewImage);
//   Print(_DNormal,"return to R.          \n");

// }






// void BackProjection_C(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
// {
//   Image *NewImage, *InvNewImage, *spectrum, *RefImage;
  
//   if (strstr(*DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
//   else DebugNiveau=_DNormal;
//   ReadIradonArgs("RadonData",*mode, *DebugLevel, InterPol, *FilterTyp, Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples); 

//   // initialization of radon-image
//   NewImage=NewFloatImage(IniFile.InFile, *M, *N,_RealArray);
//   RDoubleToImage(NewImage, InImage, *M, *N );
//   InitImage(NewImage);
  







//   // Instead of Backfilter
//   int i,m,n,mm,nn,M1,N1,OldHeight,OldWidth;
//   int XSamples1,YSamples1,CenterM,CenterN;
//   float Xmin1,Ymin1,Res,*Resx,*Resy,TempFloat,*TempPoint;
//   Image *InvMyImage;

//   Print(_DNormal,"Filter after Backproject transforming: '%s'\n",NewImage->FileName);
//   Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135

//   M1=N1=(int)((NewImage->N-1)/(float)sqrt(2))+1;
//   Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples); //65x105
  
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

//   Print(_DNormal,"Original NewImage dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135
//   BackProject(NewImage,InvMyImage);
//   Print(_DNormal,"Original InvMyImage dimensions after Backproject: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128





  
//   /* Filter the backprojected image */
//   FFTImage(InvMyImage,_FFT);
//   Print(_DNormal,"Original InvMyImage dimensions after FFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128

//   CenterM=InvMyImage->M/2;
//   CenterN=InvMyImage->N/2;
//   Resx=FloatVector(InvMyImage->M);
//   Resy=FloatVector(InvMyImage->N);
  
//   for(m=0;m<InvMyImage->M;m++) {
//     mm=m;
//     if(mm>=CenterM) mm=InvMyImage->M-m;
//     Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
//       (InvMyImage->DeltaX*InvMyImage->DeltaX);
//   }
  
//   for(n=0;n<InvMyImage->N;n++) {
//     nn=n;
//     if(nn>=CenterN) nn=InvMyImage->N-n;
//     Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
//       (InvMyImage->DeltaY*InvMyImage->DeltaY);
//   }
  
//   Print(_DNormal,"Original InvMyImage dimensions before filter: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128

//   for(m=0;m<InvMyImage->M;m++) {
//     TempFloat=Resx[m];
//     TempPoint=InvMyImage->Signal[m];
//     for(n=0,i=0;n<InvMyImage->N;n++) {
//       Res=sqrt(TempFloat+Resy[n]);
//       /*Print(_DDebug,"Res %5.5f\n",Res);*/
//       TempPoint[i++]*=Res;
//       TempPoint[i++]*=Res;
//     }
//   }

//   FFTImage(InvMyImage,_IFFT);
//   Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
//   ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
//   Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
//   RealImage(InvMyImage);
//   Print(_DNormal,"Original InvMyImage dimensions after RealImage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105

//   //NormImage(InvMyImage,1.0,-MeanValue(InvMyImage));  // Change
//   //PrintStats(_DDetail,InvMyImage);
//   Free(Resx);
//   Free(Resy);
  









//   Print(_DNormal,"Using Filtering after Backprojection\n");
//   // InvNewImage=BackFilter(NewImage);
//   // ScaleImage(InvNewImage);
//   // PrintStats(_DDetail,InvNewImage);
//   // ImageToFloat(OutImage, InvNewImage);
//   // FreeImage(InvNewImage);

//   ScaleImage(InvMyImage);
//   PrintStats(_DDetail,InvMyImage);
//   ImageToFloat(OutImage, InvMyImage);
//   FreeImage(InvMyImage);
  

  
//   FreeImage(NewImage);
//   Print(_DNormal,"return to R.          \n");

// }














// // double Uty_PET_new(Image *MyImage, double **z1)
// // {

// //   int i,m,n,mm,nn,M,N,OldHeight,OldWidth;
// //   int XSamples,YSamples,CenterM,CenterN;
// //   float Xmin,Ymin,Res,*Resx,*Resy,TempFloat,*TempPoint;

// //   double sum1 = 0.; int j;
// //   for(i=0; i<MyImage->M; i++) {
// // 		for(j=0; j<MyImage->N; j++) sum1+=SQ(MyImage->Signal[i][j]);
// // 	} // y'y

// //   Image *InvMyImage;
// //   InvMyImage = BackProjection(MyImage);  // InvMyImage = K'y

// //   // Missed FFT - BUG SUbrata
// //   FFTImage(InvMyImage,_FFT);  // InvMyImage -> V' InvMyImage = V' K'y  // Check which FFT is to be used

// //   CenterM=InvMyImage->M/2;
// //   CenterN=InvMyImage->N/2;

// //   Resx=FloatVector(InvMyImage->M);
// //   Resy=FloatVector(InvMyImage->N);
  
// //   for(m=0;m<InvMyImage->M;m++) {
// //     mm=m;
// //     if(mm>=CenterM) mm=InvMyImage->M-m;
// //     Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
// //       (InvMyImage->DeltaX*InvMyImage->DeltaX);
// //     /*Print(_DDebug,"Resx %5.5f\n ",Resx[m]);*/
// //   }
  
// //   for(n=0;n<InvMyImage->N;n++) {
// //     nn=n;
// //     if(nn>=CenterN) nn=InvMyImage->N-n;
// //     Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
// //       (InvMyImage->DeltaY*InvMyImage->DeltaY);
// //     /*Print(_DDebug,"Resy %5.5f\n",Resy[n]);*/
// //   }
  
// //   Print(_DNormal,"Original InvMyImage dimensions: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);

// //   for(m=0;m<InvMyImage->M;m++) {
// //     TempFloat=Resx[m];
// //     TempPoint=InvMyImage->Signal[m];
// //     for(n=0,i=0;n<InvMyImage->N;n++) {
// //       Res=sqrt(sqrt(TempFloat+Resy[n]));  // This is the Sqrt of previous Res

// //       /*Print(_DDebug,"Res %5.5f\n",Res);*/
// //       TempPoint[i++]*=Res;
// //       TempPoint[i++]*=Res;
// //     }
// //   }
// //   Print(_DNormal,"Original InvMyImage dimensions after filter: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);


// //   // THIS AREA IS TO BE CHECKED

// //   // z1 is InvMyImage->Signal = D^{-1/2} V' K' y;
// //   // Okay, z1 is scalar or vector?? y seems to be of dim n, not nxn!!
// //   // memcpy(z1, InvMyImage->Signal, sizeof(float)*MyImage->M*MyImage->N);



// //   // DEBUG STARTS


// //   Print(_DNormal,"MyImage->M %d\n ",MyImage->M);
// //   Print(_DNormal,"MyImage->N %d\n ",MyImage->N);

// //   Print(_DNormal,"InvMyImage->M %d\n ",InvMyImage->M);
// //   Print(_DNormal,"InvMyImage->N %d\n ",InvMyImage->N);

  
// //   // 256 x 256 -> 128 x 128
// //   for(i=0; i<InvMyImage->M; i++) {
// //     for (j=0; j<InvMyImage->N; j++){ 
// //       // z1[i][j] = InvMyImage->Signal[i][j]; 
// //       // Here would be a rearrangement according to the different FFT


// //     }
// //   }
// //   double sum2 = 0.;
// //   for(i=0;i<n;i++) {
// // 		for(j=0;j<n;j++) sum2+=SQ(z1[i][j]);
// // 	}

// //   Free(Resx);
// //   Free(Resy);

// //   return (sum1 - sum2);
// // }


// // f is of size 2n, g is of size n
// // f -> g
// // void rearrange_1(int n, double *f, double *g){
  
// //   g[0]=f[0];
// // 	if ((n % 2)==0)
// // 		f[n/2]=creal(zf[n/2]);    
	
// // 	for(i=1;i<n;i++) {
// // 		f[i]=creal(f[i]+f[2*n-i])/sqrt(2.);
// // 		f[n-i]=cimag(zf[i]-zf[n-i])/sqrt(2.); 
// // 	}

// // }



// // Test 
// void BackProjection_C_test(double *InImage, double *OutImage, char **mode, int *InterPol , char **FilterTyp, char **DebugLevel, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *M, int *N, int *XSamples, int *YSamples)
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
//   float Xmin1,Ymin1,Res,*Resx,*Resy,TempFloat,*TempPoint;
//   Image *InvMyImage;

//   // Print(_DNormal,"Filter after Backproject transforming: '%s'\n",NewImage->FileName);
//   Print(_DNormal,"Sinogram dimensions: M:%i N:%i\n",NewImage->M,NewImage->N);  // 320x135

//   // M1=N1=(int)((NewImage->N-1)/(float)sqrt(2))+1;
//   Print(_DNormal,"Backprojected image dim.: M:%i N:%i\n",IniFile.XSamples,IniFile.YSamples); //61x105
  
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

//   // // Printing eig coming from filter_new
//   // FILE *ffile;
//   // ffile=fopen("eig_in_backprojection.dat","w");
//   // for(i=0; i<InvMyImage->M; ++i){
//   //   for(int j=0; j<InvMyImage->N; ++j){
//   //     fprintf(ffile,"%f ", eig[i][j]);
//   //   }
//   // }
//   // fclose(ffile);
//   // FREE_MATRIX(eig);




  
//   /* Filter the backprojected image */
//   FFTImage(InvMyImage,_FFT);

//   CenterM=InvMyImage->M/2;
//   CenterN=InvMyImage->N/2;
//   Resx=FloatVector(InvMyImage->M);
//   Resy=FloatVector(InvMyImage->N);
  
//   for(m=0;m<InvMyImage->M;m++) {
//     mm=m;
//     if(mm>=CenterM) mm=InvMyImage->M-m;
//     Resx[m]=(((float)mm)/InvMyImage->M)*(((float)mm)/InvMyImage->M)/
//       (InvMyImage->DeltaX*InvMyImage->DeltaX);
//   }
  
//   for(n=0;n<InvMyImage->N;n++) {
//     nn=n;
//     if(nn>=CenterN) nn=InvMyImage->N-n;
//     Resy[n]=(((float)nn)/InvMyImage->N)*(((float)nn)/InvMyImage->N)/
//       (InvMyImage->DeltaY*InvMyImage->DeltaY);
//   }
  
//   for(m=0;m<InvMyImage->M;m++) {
//     TempFloat=Resx[m];
//     TempPoint=InvMyImage->Signal[m];
//     for(n=0,i=0;n<InvMyImage->N;n++) {
//       Res=sqrt(TempFloat+Resy[n]);
//       /*Print(_DDebug,"Res %5.5f\n",Res);*/
//       TempPoint[i++]*=Res;
//       TempPoint[i++]*=Res;
//     }
//   }

//   FFTImage(InvMyImage,_IFFT);
//   Print(_DNormal,"Original InvMyImage dimensions after iFFT: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 128x128
//   ShrinkImage(InvMyImage,OldHeight,OldWidth,_MiddleMiddle); 
//   Print(_DNormal,"Original InvMyImage dimensions after shrinkimage: M:%i N:%i\n",InvMyImage->M,InvMyImage->N);  // 65x105
//   RealImage(InvMyImage);

//   //NormImage(InvMyImage,1.0,-MeanValue(InvMyImage));  // Change
//   //PrintStats(_DDetail,InvMyImage);
//   Free(Resx);
//   Free(Resy);
  









//   Print(_DNormal,"Using Filtering after Backprojection\n");
//   // InvNewImage=BackFilter(NewImage);
//   // ScaleImage(InvNewImage);
//   // PrintStats(_DDetail,InvNewImage);
//   // ImageToFloat(OutImage, InvNewImage);
//   // FreeImage(InvNewImage);

//   ScaleImage(InvMyImage);
//   PrintStats(_DDetail,InvMyImage);
//   ImageToFloat(OutImage, InvMyImage);
//   FreeImage(InvMyImage);
  

  
//   FreeImage(NewImage);
//   FreeImage(NewImagecpy);
//   Print(_DNormal,"return to R.          \n");

// }

