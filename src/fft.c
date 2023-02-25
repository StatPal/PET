#include<stdio.h>
#include<math.h>
#include"array.h"
#include<complex.h>
#include<fftw3.h>

/* calculates the FFT of a complex periodic 1-, 2- and 3-dimensional arrays, 
   using FFTW, described at www.fftw.org. Note that the FFT's are scaled 
   which means that applying the FFT and then its inverse provides the 
   sequence back. This is a variation from the output of FFTW. 
   
   Written by: Ranjan Maitra, Ames IA 50014. 
   Dated: 2005/08/26 

   All rights reserved.
*/


void fft(int n,double complex *f,int sign)
{
  /* calculates the forward and backward Fast Fourier Transform of a complex 
     1-dimensional sequence (using FFTW, the Fastest Fourier Transform of the 
     West). 
     The third argument sign takes one of FFTW_FORWARD and FFTW_BACKWARD.
     Note that the input sequence is replaced by its forward or inverse FFT on
     return.  */
  
  fftw_complex *ff,*out;
  fftw_plan p;
  int i;

  if (!abs(sign)) {
    printf("ERROR: Please specify one of FFTW_FORWARD or FFTW_BACKWARD for forward/inverse transform\n");
    exit(1);   
  }

  ff = fftw_malloc(n * sizeof *ff);
  out = fftw_malloc(n * sizeof *out);

  p = fftw_plan_dft_1d(n, ff, out, sign, FFTW_MEASURE);
  for (i=0;i<n;i++) ff[i]=f[i];

  fftw_execute(p); /* repeat as needed */
  for (i=0;i<n;i++) f[i]=out[i]/sqrt((double)n);
  fftw_destroy_plan(p);
  fftw_free(ff);
  fftw_free(out);
  
}

void fft2(int m, int n, double complex **f, int sign)
{
  /* calculates the forward and backward Fast Fourier Transform of a complex 
     2-dimensional m x n array (using FFTW, the Fastest Fourier Transform of 
     the West). 
     The fourth argument sign takes one of FFTW_FORWARD and FFTW_BACKWARD.
     Note that the input array is replaced by its two-dimensional forward or 
     inverse FFT on return.  */
  
	fftw_complex *ff,*out;
	fftw_plan p;
	int i,j;
	
	if (!abs(sign)) {
		printf("ERROR: Please specify one of FFTW_FORWARD or FFTW_BACKWARD for forward/inverse transform\n");
		exit(1);
	}
	
	ff = fftw_malloc(m * n * sizeof *ff);
	out= fftw_malloc(m * n * sizeof *out);
	
	p = fftw_plan_dft_2d(m,n, ff, out, sign, FFTW_MEASURE);
	
	for (i=0;i<m;i++) {
		for(j=0;j<n;j++) ff[i*n+j]=f[i][j];
	}
	
	fftw_execute(p); /* repeat as needed */
	for (i=0;i<m;i++) {
		for(j=0;j<n;j++) f[i][j]=out[i*n+j]/sqrt((double)m*n);
	}
	fftw_destroy_plan(p);
	fftw_free(ff);
	fftw_free(out);
	
}

void fft3(int m, int n, int p, double complex ***f, int sign)
{
	/* calculates the forward and backward Fast Fourier Transform of a complex 
	   3-dimensional m x n x p array (using FFTW, the Fastest Fourier Transform of 
	   the West). 
	   The fourth argument sign takes one of FFTW_FORWARD and FFTW_BACKWARD.
	   Note that the input array is replaced by its three-dimensional forward or 
	   inverse FFT on return.  */
	
	fftw_complex *ff, *out; 
	fftw_plan pp;
	int i, j, k;
	
	if (!abs(sign)) {
		printf("ERROR: Please specify one of FFTW_FORWARD or FFTW_BACKWARD for forward/inverse transform\n");
		exit(1);
	}
	
	ff = fftw_malloc(m * n * p * sizeof *ff);
	out= fftw_malloc(m * n * p * sizeof *out);
	
	pp = fftw_plan_dft_3d(m, n, p, ff, out, sign, FFTW_MEASURE);
	
	for (i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			for (k = 0; k < p; k++) ff[i * n * p + j * p + k] = f[i][j][k];
		}
	}
	fftw_execute(pp); /* repeat as needed */
	for (i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			for (k = 0; k < p; k++) f[i][j][k] = out[i * n * p + j * p + k]/sqrt((double)m * n * p);
		}
	}
	fftw_destroy_plan(pp);
	fftw_free(ff);
	fftw_free(out);
}

