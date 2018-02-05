/******************************************************************************
[HEADER] 

This library contains different common utilities used by most of the
other functions.

******************************************************************************/

#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include "misc.h"
#include "imgtools.h"
#include "calc.h"
#include <R.h>

/*****************************************************************************
[NAME]
ReadIradonArgs

[SYNOPSIS]
void ReadIradonArgs

[DESCRIPTION]

This function copy all parameter from sending from R to  the global structure 
{\tt IniFile}.


[REVISION]
Sep. 94, JJJ and PAP\\
April 12, 96 PT typo
July 2006, J.S., completly modification of the routine for the using with R
*****************************************************************************/
void ReadIradonArgs(char *inFile,char *mode, char *DebugLevell, int *InterPol, char *FilterTyp, double *Xminn, double *Yminn, double *DeltaXn, double *DeltaYn, int *M, int *N)
{
  

  strcpy(IniFile.DebugLevel, DebugLevell);     
	 if (strstr(IniFile.DebugLevel,"Normal"))   DebugNiveau=_DNormal;
  else if (strstr(IniFile.DebugLevel,"Detail"))    DebugNiveau=_DDetail;
  else if (strstr(IniFile.DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
  else Error("Unknown DebugLevel: '%s'",IniFile.DebugLevel);
    
  strcpy(IniFile.InFile, inFile);
  strcpy(IniFile.Function, mode);
 
  //Interpolation level. Used by Filtered Backprojection.
  IniFile.InterPol=*InterPol;
  IniFile.FilterCutoff=1.0;
  //IniFile.FilterType=_Ramp;
  if (strcmp(FilterTyp,"Ramp")==0) IniFile.FilterType=_Ramp;
  else if (strcmp(FilterTyp,"Hanning")==0) IniFile.FilterType=_Hanning;
  else if (strcmp(FilterTyp,"Hamming")==0) IniFile.FilterType=_Hamming;
  else Error("Unknown FilterTyp: '%s'", FilterTyp);
  
  IniFile.SliceNumber=1;
   
  IniFile.Xmin=*Xminn;
  IniFile.Ymin=*Yminn;
  IniFile.XSamples=*M;
  IniFile.YSamples=*N;
  IniFile.DeltaX=*DeltaXn;
  IniFile.DeltaY=*DeltaYn;
  
  
}



/*******************************************************************************
[NAME]
GetDateTime

[SYNOPSIS]
void GetDateTime(char *str,
                 int DateTimeFormat);
 
[DESCRIPTION]

This function gets the current time and date, and returns it in the
pointer {\tt str}. The return string are formatted accordingly to {\tt
DateTime} Format. Three formats are avaliable as of now, {\tt
\_LongDate}, {\tt \_Time} and \tc{RealTime}. {\tt \_Longdate} returns the date and
time in a string as `Day, dd. month, hh:mm:ss', {\tt \_Time} will
return `hh:mm:ss' and \tc{ \_RealTime} returns the number of seconds
from last midnight.

[USAGE]
{\tt GetDateTime(str,\_Time);}

Returns the current time in {\tt str} formatted as ``12:34:56''.

[REVISION]
Oct. 94, JJJ
Mar 02/07 J.Schulz Memory allocation by R (Calloc)
Feb 2018  J.Schulz Replace sprintf by Rprintf
*******************************************************************************/
void GetDateTime(char *str, int DateTimeFormat)
{
  time_t *timep;
  struct tm *times;

  if (!(timep=(time_t *)Calloc(1, time_t)))
    Error(" Memory allocation error (GetDateTime)");
/*  if (!(times=(struct tm *)malloc(sizeof(struct tm))))
    Error(" Memory allocation error (GetDateTime)");*/
  time(timep);
  times=localtime(timep);
  if (DateTimeFormat==_LongDate) 
    strftime(str, 55, "%A, %d. %B, %H.%M:%S", times);
  if (DateTimeFormat==_Time) 
    strftime(str, 20,"%H.%M:%S", times);
  if (DateTimeFormat==_RealTime)
    sprintf(str,"%d",times->tm_sec+times->tm_min*60+times->tm_hour*3600);
    //Rprintf(str, times->tm_sec+times->tm_min*60+times->tm_hour*3600, " \n") 
  Free_PT(timep);
}



/*******************************************************************************
[NAME]
Print

[SYNOPSIS]
void Print(int Niveau,
           char *fmt,
           ...);

[DESCRIPTION]

This function handles all the writing to the screen. 
{\tt fmt} and {\tt ...} are the same arguments as to {\tt
print}. The variable {\tt \_Niveau} combined with the global variable
{\tt DebugNiveau} chooses what should be written. {\tt DebugNiveau}
takes precedence over {\tt Niveau}. The following choises exist

\begin{tabbing}
{\tt \_DHardCore} \== Nothing written.\\
{\tt \_DNormal}   \>= Some information on screen.\\
{\tt \_DDetail}    \>= Full information both on screen.
\end{tabbing}


[REVISION]
Oct. 94, JJJ
March 06, JS
*******************************************************************************/
void Print(int Niveau, char *fmt, ...)
{ 
  //char LogString[255];
  va_list ap;
  va_start(ap,fmt);

  //vsprintf(LogString, fmt, ap);  
  if (((DebugNiveau&(_DDetail)) && (Niveau&(_DDetail|_DNormal))) 
     || ((DebugNiveau&(_DNormal)) && (Niveau&(_DNormal)))) {
    Rvprintf(fmt,ap);
    //vprintf(fmt,ap);
    //fflush(stdout);
  }
  va_end(ap);
} 


/*******************************************************************************
[NAME]
Error

[SYNOPSIS]
void Error(char *fmt,
           ...);

[DESCRIPTION]

This function is called when a fatal error occurs. It prints the error
message on the screen, closes the log and exits. {\tt fmt} and {\tt
...} are the same arguments as to {\tt print}.

[USAGE]
{\tt Error("Memory allocation problems");}

Prints the above message, closes the log and exits.


[REVISION]
Oct. 94, JJJ and PAP
*******************************************************************************/
void Error(char *fmt, ...)
{
  char LogString[255];
  va_list ap;

  va_start(ap,fmt);
  vsprintf(LogString, fmt, ap); 
  va_end(ap);
  
  error(LogString);
  //printf(LogString);
  //printf("\n");
  //exit(1);
}


/**********************************************************
[NAME]
MultReStore

[SYNOPSIS]
void MultReStore(float *p1,float *p2)

[DESCRIPTION]
The function will multiply the two complex numbers $A$
and $B$ and store the result at the location of $A$. Here
$A=p1[0]~+~i~p1[1]$ and $B=p2[0]~+~i~p2[1]$.

[USAGE]

{\t tMultReStore(arr1,arr2);}

Preforms the multiplication {\tt arr1=arr1*arr2}.

[REVISION]
Nov. 94, JJJ and PT
**********************************************************/
void MultReStore(float *p1,float *p2)
{
  multtemp=p1[0]*p2[0]-p1[1]*p2[1];
  p1[1]=p1[0]*p2[1]+p1[1]*p2[0];
  p1[0]=multtemp;
}


/**********************************************************
[NAME]
MultNew

[SYNOPSIS]
void MultNew(float *p1,float *p2,float *p3)

[DESCRIPTION] The function will multiply the two complex numbers $A$
and $B$, so that $C=AB$. The result is stored at the location of
$C$. Here $A=p1[0]~+~i~p1[1]$, $B=p2[0]~+~i~p2[1]$ and
$C=p3[0]~+~i~p3[1]$.

[USAGE]
{\tt MultNew(arr1,arr2,arr3);}

Preforms the multiplication {\tt arr3=arr1*arr2}.

[REVISION]
Nov. 94, JJJ and PT 
**********************************************************/
void MultNew(float *p1,float *p2,float *p3)
{
  p3[0]=p1[0]*p2[0]-p1[1]*p2[1];
  p3[1]=p1[0]*p2[1]+p1[1]*p2[0];
}


/********************************************************************************
[NAME]
FloatVector

[SYNOPSIS]
float *FloatVector(int Size);

[DESCRIPTION]

Allocates and returns a vector of floats, with \tc{Size} number of
elements. The vector is initialized to 0. Checking is made to ensure
no memory problems.

[USEAGE]
{\tt Test=FloatVector(99);}

Allocates \tc{Test} as a float vector 99 elements long.

[REVISION]
Oct. 94, JJJ
Mar 02/07 J.Schulz Memory allocation by R (calloc --> Calloc)
*****************************************************************************/
float *FloatVector(int Size)
{
  float *data;

  if (!(data=(float*)Calloc(Size, float)))
    Error("Memory allocation problems (FloatVector). %i elements.",Size);
  return data;
}


/*****************************************************************************
[NAME]
IntVector

[SYNOPSIS]
float *IntVector(int Size);

[DESCRIPTION]

Allocates and returns a vector of integers, with \tc{Size} number of
elements. The vector is initialized to 0. Adequate checking is done.

[USEAGE]
{\tt Test=IntVector(99);}

Allocates \tc{Test} as an int vector 99 elements long.

[REVISION]
Oct. 94, JJJ
Mar 02/07 J.Schulz Memory allocation by R (calloc --> Calloc)
*****************************************************************************/
int *IntVector(int Size)
{
  int *data;

  if (!(data=(int*)Calloc(Size, int)))
    Error("Memory allocation problems (IntVector). %i elements",Size);
  return data;
}


/*****************************************************************************
[NAME]
convert2

[SYNOPSIS]
void convert2(char *ind);

[DESCRIPTION]

The function swaps the two bytes at the location {\tt ind}. Used for
converting short int's between HP/PC.

[USEAGE]
{\tt convert2((char *)\&X);}

Converts {\tt X}, where {\tt X} is a short int.

[REVISION]
Oct. 94, PAP and PT
*****************************************************************************/
void convert2(char *ind)
{
  char tmp;
 
  tmp=ind[0];
  ind[0]=ind[1];
  ind[1]=tmp;
}


/******************************************************************************
[NAME]
convert4

[SYNOPSIS]
void convert4(char *ind);

[DESCRIPTION]

The function swaps the four bytes at the location {\tt ind} around the
center. Used for converting short floats between HP/PC.

[USEAGE]
{\tt convert4((char *)\&X);}

Converts {\tt X}, where {\tt X} is a float.

[REVISION]
Oct. 94, PAP and PT
*****************************************************************************/
void convert4(char *ind)
{
  char tmp;

  tmp=ind[0];
  ind[0]=ind[3];
  ind[3]=tmp;
  tmp=ind[1];
  ind[1]=ind[2];
  ind[2]=tmp;
}


/*****************************************************************************
[NAME]
convert8

[SYNOPSIS]

void convert8(char *ind);

[DESCRIPTION]

The function swaps the eight bytes at the location {\tt ind} around
the center. Used for converting doubles between HP/PC.

[USEAGE]
{\tt convert8((char *)\&X);}

Converts {\tt X}, where {\tt X} is a double.

[REVISION]
Oct. 94, PAP and PT
*****************************************************************************/
void convert8(char *ind)
{
  char tmp;

  tmp=ind[0];
  ind[0]=ind[7];
  ind[7]=tmp;
  tmp=ind[1];
  ind[1]=ind[6];
  ind[6]=tmp;
  tmp=ind[2];
  ind[3]=ind[5];
  ind[5]=tmp;
  tmp=ind[3];
  ind[3]=ind[4];
  ind[4]=tmp;
}


/*****************************************************************************
[NAME]
convert\_UNIX\_PC

[SYNOPSIS]
void convert_UNIX_PC(char *ind, int lgd);

[DESCRIPTION]

The function swaps the {\tt lgd} bytes at the location {\tt ind}
around the center. Used for converting arbitrary length numbers between
HP/PC.

[USEAGE]
{\tt convert\_UNIX\_PC((char *)\&X, 8);}

Converts {\tt X}, where {\tt X} is 8 bytes long eg. a double.

[REVISION]
Oct. 94, PAP and PT
****************************************************************************/
void convert_UNIX_PC(char *ind, int lgd)
{
  int i,lgdh,lgdm1;
  char tmp;

  if (lgd>1)
  {
    lgdh=lgd/2;
    lgdm1=lgd-1;
    for (i=0;i<lgdh;i++)
    {
      tmp=ind[i];
      ind[i]=ind[lgdm1-i];
      ind[lgdm1-i]=tmp;
    }
  }
}


/*****************************************************************************
[NAME]
strequal

[SYNOPSIS]
int strequal(char *str1, char *str2);

[DESCRIPTION]
This function will compare two strings by length and see whether one is
contained in the other.

[USEAGE]
{\tt ok=strequal("Radon","Radon transform");}
This will return a zero. Drop transform and it will return a 1.

[REVISION]
April 9, 96 PToft
****************************************************************************/
int strequal(char *str1, char *str2)
{
  if ((strlen(str1)==strlen(str2)) && (strstr(str1,str2)==str1))
    return 1;
  else
    return 0;
}
