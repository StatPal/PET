/* This routine is from Peter Toft and was modified from Joern Schulz. 
   
   Peter Toft, "The Radon Transform - Theory and Implementation", Ph.D. thesis,
   Department of Mathematical Modelling, Technical University of Denmark, 
   June 1996, 326 pages. See also: http://pto.linux.dk/PhD/
   
   Jörn Schulz, July 2006, Humboldt-Universität zu Berlin.
    
   NN = Nearest Neighbour

   The calling syntax in R is:

   hdata <- .C("Hough1", 
                y = matrix(0, nrow=setpar[4,1], ncol=setpar[6,1]), 
                as.double(odata), 
                as.double(setpar))$y

*/

/******************** Includes Header-Files *********************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

/*static void oHoughO1(double *v,double *u,double *setpar);
static void oHoughO2(double *v,double *u,double *setpar);
static void oHoughO3(double *v,double *u,double *setpar);
static void oHoughO4(double *v,double *u,double *setpar);*/


/******************** Main Routine C-Routine ********************/
void Hough4(double *v,double *u,double *setpar)
{

  int m,n, M;
  int R;
  int t, T,moffset,noffset;
  double x_min,Delta_x,Delta_rho,rho_min,Delta_theta;
  double *costheta,*sintheta,theta,x,y;
  double theta_min;
  double sigval;
  double *xc,*ys;
  double cst;
  
  M            = setpar[0];
  x_min        = setpar[1];
  Delta_x      = setpar[2];
  T            = setpar[3];
  Delta_theta  = setpar[4];
  R            = setpar[5];
  rho_min      = setpar[6];
  Delta_rho    = setpar[7];
  theta_min    = setpar[8];

  for (t=0;t<R*T;t++)
    v[t] = 0.0;

  if (R<=Delta_x/Delta_rho*sqrt(2.0)*(M-1))
  {
    Rprintf("Number of samples in the rho-direction R is to small.\n");
  }
  else
  {
    costheta=(double *)Calloc(T, double);
    sintheta=(double *)Calloc(T, double);
    xc=(double *)Calloc(M*T, double);
    ys=(double *)Calloc(M*T, double);
    
    for (t=0;t<T;t++)
    {
      theta=t*Delta_theta;
      costheta[t]=cos(theta);
      sintheta[t]=sin(theta);
    }
    
    for (m=0;m<M;m++)
    {
      x=(x_min+m*Delta_x)/Delta_rho; /* x(m) = y(m) */    
      cst=rho_min/Delta_rho+0.5;
      for (t=0;t<T;t++)
      {
        xc[m*T+t]=x*costheta[t];
        ys[m*T+t]=x*sintheta[t]-cst;
      }
    }
    for (m=0;m<M;m++)
    {
      x=x_min+m*Delta_x;
      moffset=m*T;
      for (n=0;n<M;n++)
      {
        sigval=u[m+n*M];
        if (fabs(sigval)>1E-30)
        {
          y=x_min+n*Delta_x;
          noffset=n*T;
          for (t=0;t<T;t++)
            v[t+T*(int)(xc[moffset+t]+ys[noffset+t])]+=sigval;
        }
      }
    }
    Free(costheta);
    Free(sintheta);
    Free(xc);
    Free(ys);
  }
}

void Hough3(double *v,double *u,double *setpar)
{

  int m,n, M;
  int R;
  int t, T;
  double x_min,Delta_x,Delta_rho,rho_min,Delta_theta;
  double *costheta,*sintheta,theta,x,y;
  double theta_min;
  double sigval,offset;

  M            = setpar[0];
  x_min        = setpar[1];
  Delta_x      = setpar[2];
  T            = setpar[3];
  Delta_theta  = setpar[4];
  R            = setpar[5];
  rho_min      = setpar[6];
  Delta_rho    = setpar[7];
  theta_min    = setpar[8];

  for (t=0;t<R*T;t++)
    v[t] = 0.0;

  if (R<=Delta_x/Delta_rho*sqrt(2.0)*(M-1))
    Rprintf("Number of samples in the rho-direction R is to small\n");
  else
  {
    costheta=(double *)Calloc(T, double);
    sintheta=(double *)Calloc(T, double);
    for (t=0;t<T;t++)
    {
      theta=t*Delta_theta;
      costheta[t]=cos(theta);
      sintheta[t]=sin(theta);
    }
    
    offset=0.5-rho_min/Delta_rho;
    for (m=0;m<M;m++)
    {
      x=(x_min+m*Delta_x)/Delta_rho;
      for (n=0;n<M;n++)
      {
        y=(x_min+n*Delta_x)/Delta_rho;
        sigval=u[m+n*M];
        if (fabs(sigval)>1E-30)
          for (t=0;t<T;t++)
            v[t+T*(int)(offset+x*costheta[t]+y*sintheta[t])]+=sigval;
      }
    }
    Free(sintheta);
    Free(costheta);
  }
}

void Hough2(double *v,double *u,double *setpar)
{

  int m,n, M;
  int r, R;
  int t, T;
  double x_min,Delta_x,Delta_rho,rho_min,Delta_theta;
  double *costheta,*sintheta,theta,x,y;
  double theta_min;
  double sigval,off;

  M            = setpar[0];
  x_min        = setpar[1];
  Delta_x      = setpar[2];
  T            = setpar[3];
  Delta_theta  = setpar[4];
  R            = setpar[5];
  rho_min      = setpar[6];
  Delta_rho    = setpar[7];
  theta_min    = setpar[8];

  costheta=(double *)Calloc(T, double);
  sintheta=(double *)Calloc(T, double);

  for (t=0;t<R*T;t++)
    v[t] = 0.0;

  for (t=0;t<T;t++)
  {
    theta=t*Delta_theta;
    costheta[t]=cos(theta);
    sintheta[t]=sin(theta);
  }

  off=0.5-rho_min/Delta_rho;
  for (m=0;m<M;m++)
  {
    x=(x_min+m*Delta_x)/Delta_rho;
    for (n=0;n<M;n++)
    {
      y=(x_min+n*Delta_x)/Delta_rho;
      sigval=u[m+n*M];
      if ((sigval>1E-30) || (sigval<-1E-30))
        for (t=0;t<T;t++)
        {
          r=(int)floor(off+x*costheta[t]+y*sintheta[t]);
          if ((r>=0) && (r<R)) v[t+r*T]+=sigval;
        }
    }
  }
  Free(sintheta);
  Free(costheta);
}

void Hough1(double *v,double *u,double *setpar)
{

  int m,n, M;
  int r, R;
  int t, T;
  double x_min,Delta_x,Delta_rho,rho_min,Delta_theta;
  double theta,x,y;
  double theta_min;
  double sigval;

  M            = setpar[0];
  x_min        = setpar[1];
  Delta_x      = setpar[2];
  T            = setpar[3];
  Delta_theta  = setpar[4];
  R            = setpar[5];
  rho_min      = setpar[6];
  Delta_rho    = setpar[7];
  theta_min    = setpar[8];

  for (t=0;t<R*T;t++)
    v[t] = 0.0;

  for (m=0;m<M;m++)
  {
    x=x_min+m*Delta_x;
    for (n=0;n<M;n++)
    {
      y=x_min+n*Delta_x;
      sigval=u[m+n*M];
      if (fabs(sigval)>1E-30)
        for (t=0;t<T;t++)
        {
          theta=t*Delta_theta;
          r=(int)floor(0.5+(x*cos(theta)+y*sin(theta)-rho_min)/Delta_rho);
          if ((r>=0) && (r<R)) v[t+r*T]+=sigval;
        }
    }
  }
}




