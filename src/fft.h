#include<stdio.h>
#include<math.h>
#include<complex.h>
#include"array.h"
#include<fftw3.h>

void fft(int n,double complex *f,int sign);
void fft2(int m,int n,double complex **f,int sign);
void fft3(int m,int n, int p, double complex ***f,int sign);
