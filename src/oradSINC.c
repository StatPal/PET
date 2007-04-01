/* This routine is from Peter Toft and was modified from Joern Schulz. 
   
   Peter Toft, "The Radon Transform - Theory and Implementation", Ph.D. thesis,
   Department of Mathematical Modelling, Technical University of Denmark, 
   June 1996, 326 pages. See also: http://pto.linux.dk/PhD/
   
   July 2006, Joern Schulz, Humboldt-Universität zu Berlin.
    
   NN = Nearest Neighbour

   The calling syntax in R is:

       rdata <- .C("radonSINC", 
                   y = matrix(0, nrow=setpar[4,1], ncol=setpar[6,1]), 
                   as.double(odata), 
                   as.double(setpar))$y

*/

/******************** Includes Header-Files *********************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include "misc.h"

/********************* The c-programm-code **********************/

/* Subroutine for Mainroutine ordaon */
double sinc(double x)
{

  if (x<0) x=-x;
  if (x<1e-10) 
    return 1.0;
  else
    return sin(x)/x;
}

/* Main Routine */
void radonSINC(double *v, double *u, double *setpar)
{
  int m,n, M;
  int r, R;
  int t, T;
  // int mmin,mmax,nmin,nmax,adr;
  // double dx_isintheta,dx_icostheta,x,y
  double sum,x_min,Delta_x,Delta_rho,rho_min,Delta_theta;
  double costheta,sintheta,theta;
  double rho,rhooffset;
  double theta_min;
  double idx,mintrig,psi;
  double *xar,scal;

  M            = setpar[0];
  x_min        = setpar[1];
  Delta_x      = setpar[2];
  T            = setpar[3];
  Delta_theta  = setpar[4];
  R            = setpar[5];
  rho_min      = setpar[6];
  Delta_rho    = setpar[7];
  theta_min    = setpar[8];

  scal=M_PI/Delta_x;
  xar=(double *)R_alloc(M, sizeof(double));

  for ( m = 0 ; m < M ; m++ ) 
    xar[m]=m*Delta_x+x_min;

  idx=1.0/Delta_x;
  for(t=0; t<T;t++)
  {
    Print(_DNormal,"t=%i (%i)\n",t,T);
    theta=t*Delta_theta+theta_min;
    sintheta=sin(theta);
    costheta=cos(theta);
    rhooffset=x_min*(costheta+sintheta);
    if (sintheta>sqrt(0.5))
      mintrig=1.0/sintheta;
    else
      mintrig=1.0/fabs(costheta);

    for(r=0; r<R; r++ )
    {
      rho=r*Delta_rho+rho_min;
      sum=0.0;
      for (m=0;m<M;m++)
      {
        rhooffset=rho-xar[m]*costheta;
        for (n=0;n<M;n++)
        {
          psi=scal*(rhooffset-xar[n]*sintheta);
          sum+=sinc(mintrig*psi)*u[m+n*M];
	  /* sum+=psi*u[m+n*M];*/
        }
      }
      v[t+T*r]+=sum*mintrig*Delta_x;
    }
  }
}
