/* This routine is from Peter Toft and was modified from Joern Schulz. 
   
   Peter Toft, "The Radon Transform - Theory and Implementation", Ph.D. thesis,
   Department of Mathematical Modelling, Technical University of Denmark, 
   June 1996, 326 pages. See also: http://pto.linux.dk/PhD/
   
   July 2006, Joern Schulz, Humboldt-Universität zu Berlin.
    
   NN = Nearest Neighbour

   The calling syntax in R is:

       rdata <- .C("radonLI", 
                   y = matrix(0, nrow=setpar[4,1], ncol=setpar[6,1]), 
                   as.double(odata), 
                   as.double(setpar))$y

*/

/******************** Includes Header-Files *********************/

#include <stdio.h>
#include <math.h>


/******************** Main Routine C-Routine ********************/
void radonLI(double *v,double *u,double *setpar)
{
  int m,n, M;
  int r, R;
  int t, T;
  int mmin,mmax,nmin,nmax,adr;
  double sum,x_min,Delta_x,Delta_rho,rho_min,Delta_theta;
  double costheta,sintheta,theta;
  //double rho,x,y;
  double dx_isintheta,dx_icostheta,rhooffset;
  double alpha,beta,theta_min;
  double eps=1e-9,mfloat,nfloat,reldx;
  double idx;
  M            = setpar[0];
  x_min        = setpar[1];
  Delta_x      = setpar[2];
  T            = setpar[3];
  Delta_theta  = setpar[4];
  R            = setpar[5];
  rho_min      = setpar[6];
  Delta_rho    = setpar[7];
  theta_min    = setpar[8];

  idx=1.0/Delta_x;
  for(t=0; t<T;t++)   // T=ma?
  {
    theta=t*Delta_theta+theta_min;
    sintheta=sin(theta);
    costheta=cos(theta);
    rhooffset=x_min*(costheta+sintheta);  // Common factor
    if (sintheta>sqrt(0.5))               // Project onto x axis
    {
      dx_isintheta=Delta_x/sintheta;
      alpha=-costheta/sintheta;           // Digital slope is set
      for(r=0; r<R; r++ )                 // For all values of rho
      {
        beta=(r*Delta_rho+rho_min-rhooffset)/(Delta_x*sintheta);    // Offset is set

        // mmin and mmax are set using eqn 1.19
        if (alpha>eps)
        {
          mmin=(int)ceil((eps-beta)/alpha);
          mmax=1+(int)floor((M-beta-1-eps)/alpha);
        }
        else
          if (alpha<-eps)
          {
            mmin=(int)ceil((M-beta-1-eps)/alpha);
            mmax=1+(int)floor((eps-beta)/alpha);
          }
          else
            if ((beta>1) && (beta<M-2))
            {
              mmin=0;
              mmax=M;
            }
            else
            {
              mmin=0;
              mmax=-1;
            }
        if (mmin<0) mmin=0;
        if (mmax>M) mmax=M;
        // mmin and mmax set


        sum=0.0;
        for (m=mmin;m<mmax;m++)   // For all legal values of x
        {
          nfloat=beta+m*alpha;
          n=(int)nfloat;
          reldx=(nfloat-n)*idx;
          adr=m+M*n;
          sum+=u[adr]*(1-reldx)+u[adr+M]*reldx;  // Increment sum
        }
        v[t+T*r]=sum*dx_isintheta;                // Update matrix element
      }
    }
    else                                          // Project onto y-axis
    {
      dx_icostheta=Delta_x/fabs(costheta);        // DIFFERENT
      alpha=-sintheta/costheta;                   // DIFFERENT
      for(r=0; r<R; r++ ) 
      {
        beta=(r*Delta_rho+rho_min-rhooffset)/(Delta_x*costheta);    // DIFFERENT
        if (alpha>eps)
        {
          nmin=(int)ceil((eps-beta)/alpha);
          nmax=1+(int)floor((M-beta-1-eps)/alpha);
        }
        else 
          if (alpha<-eps)
          {
            nmin=(int)ceil((M-beta-1-eps)/alpha);
            nmax=1+(int)floor((eps-beta)/alpha);
          }
          else
            if ((beta>1) && (beta<M-2))
            {
              nmin=0;
              nmax=M;
            }
            else
            {
              nmin=0;
              nmax=-1;
            }
        if (nmin<0) nmin=0;
        if (nmax>M) nmax=M;

        sum=0.0;
        for (n=nmin;n<nmax;n++)
        {
          mfloat=beta+n*alpha;
          m=(int)mfloat;
          reldx=(mfloat-m)*idx;
          adr=m+M*n;
          sum+=u[adr]*(1-reldx)+u[adr+1]*reldx;
        }
        v[t+T*r]=sum*dx_icostheta;                // DIFFERENT
      }
    }
  }
}
