/* itgetpar.c */

#include "it.h"

char Algorithmar[3][50]={"Algebraric Reconstruction",
                           "Expection Maximization",
                           "Conjugate Gradients"};

char RadonKernelar[6][50]={"Sinc interpolation",
                             "Binary Nearest neighboar",
                             "Ray driven nearest neighboar",
                             "Ray driven linear interpolation",
                             "Pixel driven (rays)",
                             "Pixel driven (for all points)"};

char TrueFalsear[2][6]={"False",
                          "True"};

extern itINItype itINI;

/*****************************************************************************
[NAME]
ReadItArgs

[SYNOPSIS]
void ReadIt

[DESCRIPTION]

This function copy all parameter from sending from R to  the global structure 
{\tt IniFile}.


[REVISION]
April 96 JJJ and PT
July 2006, J.S., completly modification of the routine for the using with R
*****************************************************************************/
void ReadItArgs(char *inFileName, char *Mode, int *UseFast, char *RadonKernel, char *IterationsType, int *Iterations, int *SaveIterations, char *SaveIterationsName, double *LowestALevel, double *ConstrainMin, double *ConstrainMax, double *Alpha, double *Beta, double *Regularization, int *KernelFileSave, char *KernelFileName, char *RefFileName, int *ThetaSamples, double *ThetaMin, double *DeltaTheta, int *RhoSamples, double *RhoMin, double *DeltaRho, double *Xmin, double *Ymin, double *DeltaX, double *DeltaY, int *XSamples, int *YSamples, int *OverSamp, char *DebugLevel)
{
  char Temp[100];
  char *strptr;
  //Image *tempImage, *refImage=NULL;
  
  // ******************************
  // DebugLevel
  strcpy(Temp, DebugLevel); 
  if (strstr(Temp,"Normal"))             DebugNiveau=_DNormal;
  else if (strstr(Temp,"Detail"))        DebugNiveau=_DDetail;
  else if (strstr(Temp,"HardCore"))      DebugNiveau=_DHardCore;
  else Error("Unknown DebugLevel: `%s'",Temp);
  
  // ******************************
  // names
  strcpy(itINI.InFileName, inFileName);
  Print(_DDetail,"`inFileName' entry is set to `%s' \n", itINI.InFileName);
  strcpy(itINI.RefFileName, RefFileName);
  Print(_DDetail,"`RefFileName' entry is set to `%s'\n", 
        itINI.RefFileName);
  strcpy(itINI.OutFileName,"OutputData");
  
  // ******************************
  // Input sampling parameters
  Print(_DDetail,"`ThetaSamples' entry is set to %i\n",(itINI.ThetaSamples=*ThetaSamples));
  Print(_DDetail,"`ThetaMin' entry is set to %f\n",(itINI.ThetaMin=(float) *ThetaMin));
  Print(_DDetail,"`DeltaTheta' entry is set to %f\n",
                 (itINI.DeltaTheta=(float) *DeltaTheta));
  Print(_DDetail,"`RhoSamples' entry is set to %i\n",(itINI.RhoSamples=*RhoSamples));
  Print(_DDetail,"`RhoMin' entry is set to %f\n",(itINI.RhoMin=(float) *RhoMin));
  Print(_DDetail,"`DeltaRho' entry is set to %f\n",(itINI.DeltaRho=(float) *DeltaRho));
  
  // ******************************
  // Output sampling parameters
  Print(_DDetail,"`XSamples' entry is set to %i\n",(itINI.XSamples=*XSamples));
  Print(_DDetail,"`YSamples' entry is set to %i\n",(itINI.YSamples=*YSamples));
  Print(_DDetail,"`Xmin' entry is set to %f\n",(itINI.Xmin=(float) *Xmin));
  Print(_DDetail,"`Ymin' entry is set to %f\n",(itINI.Ymin=(float) *Ymin));
  Print(_DDetail,"`DeltaX' entry is set to %f\n",(itINI.DeltaX=(float) *DeltaX));
  Print(_DDetail,"`DeltaY' entry is set to %f\n",(itINI.DeltaY=(float) *DeltaX));
    
  // ******************************
  // LowestALevel ConstrainMin ConstrainMax (all optional)
  Print(_DDetail,"`LowestALevel' entry is set to %f\n",
                 (itINI.LowestALevel=(float)  *LowestALevel));
  Print(_DDetail,"`ConstainMin' entry is set to %f\n",
                 (itINI.ConstrainMin=(float) *ConstrainMin));
  Print(_DDetail,"`ConstainMax' entry is set to %f\n",
                 (itINI.ConstrainMax=(float) *ConstrainMax));

  // ******************************
  // Iterations SaveIterations
  Print(_DDetail,"`Iterations' entry is set to %i\n",itINI.Iterations=*Iterations);
  Print(_DDetail,"`SaveIterations' entry is set to %i\n", itINI.SaveIterations=*SaveIterations);
  if (itINI.SaveIterations) 
  {
    strcpy(itINI.SaveIterationsName, SaveIterationsName);
    if (itINI.Iterations/itINI.SaveIterations>MaxSaves) {
       Print(_DNormal,"To large number of saves requested, aborting SaveIterations\n");
                      itINI.SaveIterations=0;
    }
    else
       Print(_DDetail,"Saving temporary result every %i iterations\n",
                      itINI.SaveIterations);
  }
  else Print(_DDetail,"No temporary saving \n");
  
  
  // ******************************
  // Alpha Beta Regularization
  Print(_DDetail,"`Alpha' entry is set to %f\n",(itINI.Alpha=(float) *Alpha));
  Print(_DDetail,"`Beta' entry is set to %f\n",(itINI.Beta=(float) *Beta));
  Print(_DDetail,"`Regularization' entry is set to %f\n",
                 (itINI.Regularization=(float) *Regularization));

  // ******************************
  // KernelFileSave KernelFileName
  Print(_DDetail,"`KernelFileSave' entry is set to %i\n",(itINI.KernelFileSave= *KernelFileSave));
  
  //strcpy(itINI.KernelFileName, "XTest.sif");
  strcpy(itINI.KernelFileName, KernelFileName);
  if (itINI.KernelFileSave){
      if (!(strptr=strstr(itINI.KernelFileName,".sif")))      
          Error("Error in IniFile: `%s'", itINI.KernelFileName);
      itINI.KernelFileName[strptr-itINI.KernelFileName]='\0';
     Print(_DDetail,"`KernelFileName' entry is set to `%s'\n",itINI.KernelFileName);
  }
  
        
  // ******************************
  // IterationType Mode
  if (strcmp(IterationsType,"cyclic")==0)        itINI.IterationType = 1;
  else if (strcmp(IterationsType,"random")==0)   itINI.IterationType = 0;
  else Error("Unknown IterationsType: `%s'",IterationsType); 
  Print(_DDetail,"`Numbering' entry is set to %i\n",itINI.IterationType);
  
  if (strcmp(Mode,"EM")==0)                      itINI.Algorithm=_EM;
  else if (strcmp(Mode,"ART")==0)                itINI.Algorithm=_ART;
  else if (strcmp(Mode,"CG")==0)                 itINI.Algorithm=_CG;
  else Error("Unknown Mode: `%s'",Mode);   
  Print(_DDetail,"`Algorithm' entry is set to %s\n",Algorithmar[itINI.Algorithm]);

  // ******************************
  // UseFast RadonKernel 
  Print(_DDetail,"`UseFast' entry is set to %i\n",(itINI.IsFast=*UseFast));

  if (strcmp(RadonKernel,"SINC")==0)             itINI.RadonKernel=_SINC;
  else if (strcmp(RadonKernel,"NN")==0)          itINI.RadonKernel=_NN;
  else if (strcmp(RadonKernel,"RNN")==0)         itINI.RadonKernel=_RNN;
  else if (strcmp(RadonKernel,"RL")==0)          itINI.RadonKernel=_RL;
  else if (strcmp(RadonKernel,"P1")==0)          itINI.RadonKernel=_P1;
  else if (strcmp(RadonKernel,"P2")==0)          itINI.RadonKernel=_P2;
  else Error("Unknown RadonKernel: `%s'",RadonKernel);   
  Print(_DDetail,"`RadonKernel' entry is set to %s\n",
                 RadonKernelar[itINI.RadonKernel]);
  
  // ******************************
  // OverSamp
  Print(_DDetail,"`OverSamp' entry is set to %i\n", (itINI.OverSamp=*OverSamp));
  Print(_DDetail,"\n");
  
}
