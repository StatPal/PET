/* Misc.h */

typedef struct {
  char  InFile[100];
  char  Function[50];
  char  DebugLevel[50];
  float Param[10];
  float Xmin;
  float Ymin;
  float DeltaX;
  float DeltaY;
  int   XSamples;
  int   YSamples;
  int   InterPol;
  int   FilterType;
  float  FilterCutoff;
  int   SliceNumber;
} INI;

#define _DHardCore  1
#define _DNormal    2
#define _DDetail    4

#define _LongDate 100
#define _Time     101
#define _RealTime 102

#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)<(b) ? (a) : (b))
#define sq(a) ((a)*(a))

extern int DebugNiveau;
extern INI IniFile;
extern int DebugNiveau;
extern float multtemp;

extern void ReadIradonArgs(char *, char *, char *, int *, char *, double *, double *, double *, double *, int *, int *);
extern void  GetDateTime(char *str, int DateTimeFormat);

extern void Print(int, char*, ...);
extern void Error(char*, ...);

extern void MultNew(float *p1,float *p2,float *p3);
extern void MultReStore(float *p1,float *p2);
extern float *FloatVector(int Size);
extern int *IntVector(int Size);

extern void convert2(char *);
extern void convert4(char *);
extern void convert8(char *);
extern void convert_UNIX_PC(char *ind, int lgd);
int strequal(char *, char *);










