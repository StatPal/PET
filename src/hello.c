#include "hello.h"

void test_vec_C(double *a, double *b, int *na, double *val){

    // double tmp;
    for(int i = 0; i < *na; i++){
        a[i] += *val;
        // tmp = a[i];
        // printf("%i %f   ", i, a[i]);
    }
    *b += *val;
}

#include <math.h>

void Gausfilter_test(int *n1, double *fwhm, double *f)
{
	/* one-dimensional Gaussian kernel, given in terms of FWHM pixels */
	int i;
    int n = *n1;
	double sum=0,spar=*fwhm/2.3548;

	if (spar < 1e-16) {
		f[0] = 1.;
		for (i=1; i<n; i++) f[i]=0.;
	}
	else {
		for(i=0;i<=n/2;i++) f[i] = exp(-0.5*i*i/(spar*spar));
		for(i=(n/2+1);i<n;i++) f[i]=f[(n-i)];
		sum=0;
		for(i=0;i<n;i++) sum+=f[i];
		for(i=0;i<n;i++) f[i]/=sum;
	}
}
