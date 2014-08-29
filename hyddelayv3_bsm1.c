/*
 * hyddelayv3_bsm2 is a C-file S-function for first order reaction of flow and conc.  
 * In this version the loads are first calculated and the first
 * order reaction is used for the load and flow. After this the concentrations are
 * recalculated based on the delayed flow and load. Note that the state
 * vector internally here represents mass. State no 14 (TSS) is a dummy state to 
 * maintain consistence with the normal vector size used in reactors.
 * For temperature T(out) = T(in) independently of parameter TEMPMODEL.
 *  
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *
 *
 * The file has been modified in order to include two step nitrification and four step denitrification
 * according to the principles stated in Hiatt et al., 2008
 *
 * Copyright: Xavier Flores-Alsina, modelEAU, Universite Laval, Quebec, Canada
 *                                  IEA, Lund University, Lund, Sweden 
 *
 */
 

#define S_FUNCTION_NAME hyddelayv3_bsm1

#include "simstruc.h"

#define XINIT   ssGetArg(S,0)
#define PAR     ssGetArg(S,1)
#define T       ssGetArg(S,2)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 26);  /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 26);  /* number of inputs                      */
    ssSetNumOutputs(       S, 26);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 0);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 3);   /* number of input arguments             */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements   */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements*/
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements*/
}

/*
 * mdlInitializeSampleTimes - initialize the sample times array
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}


/*
 * mdlInitializeConditions - initialize the states
 */
static void mdlInitializeConditions(double *x0, SimStruct *S)
{
int i;

for (i = 0; i < 26; i++) {
   x0[i] = mxGetPr(XINIT)[i];
}
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
  double X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS;
  int i,j ;

  X_I2TSS = mxGetPr(PAR)[19];
  X_S2TSS = mxGetPr(PAR)[20];
  X_BH2TSS = mxGetPr(PAR)[21];
  X_BA2TSS = mxGetPr(PAR)[22];
  X_P2TSS = mxGetPr(PAR)[23];

  for (i = 0; i < 13; i++) {
      y[i] = x[i]/x[14];
  }

  y[13] = (X_I2TSS*x[2]+X_S2TSS*x[3]+X_BH2TSS*x[4]+X_BA2TSS*x[5]+X_P2TSS*x[6])/x[14];
  y[14] = x[14]; /* Flow rate */
  
  y[15] = x[15]; /* Temp */
  
  /* dummy states */
  
  
  for (j = 16; j < 26; j++) {
   
      if (x[j]/x[14] > 0.0)
          y[j] =  x[j]/x[14];
      else
          y[j] =  0;
      
       }
       
  
//   y[16] = x[16]/x[14];
//   y[17] = x[17]/x[14];
//   y[18] = x[18]/x[14];
//   y[19] = x[19]/x[14];
//   y[20] = x[20]/x[14];
//   
//   y[21] = x[21]/x[14];
//   y[22] = x[22]/x[14];
//   y[23] = x[23]/x[14];
//   y[24] = x[24]/x[14];
//   y[25] = x[25]/x[14];
}

/*
 * mdlUpdate - perform action at major integration time step
 */

static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}

/*
 * mdlDerivatives - compute the derivatives
 */
static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{
int i;
double timeconst;

timeconst = mxGetPr(T)[0];
if (timeconst > 0.000001) {
  dx[0] = (u[0]*u[14]-x[0])/timeconst;
  dx[1] = (u[1]*u[14]-x[1])/timeconst;
  dx[2] = (u[2]*u[14]-x[2])/timeconst;
  dx[3] = (u[3]*u[14]-x[3])/timeconst;
  dx[4] = (u[4]*u[14]-x[4])/timeconst;
  dx[5] = (u[5]*u[14]-x[5])/timeconst;
  dx[6] = (u[6]*u[14]-x[6])/timeconst;
  dx[7] = (u[7]*u[14]-x[7])/timeconst;
  dx[8] = (u[8]*u[14]-x[8])/timeconst;
  dx[9] = (u[9]*u[14]-x[9])/timeconst;
  dx[10] = (u[10]*u[14]-x[10])/timeconst;
  dx[11] = (u[11]*u[14]-x[11])/timeconst;
  dx[12] = (u[12]*u[14]-x[12])/timeconst;
  
  dx[13] = (u[13]*u[14]-x[13])/timeconst; /* TSS */
  
  dx[14] = (u[14]-x[14])/timeconst;       /* Flow rate */
  
  dx[15] = (u[15]-x[15])/timeconst;       /* Temp */
  
  /* dummy states */
  dx[16] = (u[16]*u[14]-x[16])/timeconst;
  dx[17] = (u[17]*u[14]-x[17])/timeconst;
  dx[18] = (u[18]*u[14]-x[18])/timeconst;
  dx[19] = (u[19]*u[14]-x[19])/timeconst;
  dx[20] = (u[20]*u[14]-x[20])/timeconst;
  
  dx[21] = (u[21]*u[14]-x[21])/timeconst;
  dx[22] = (u[22]*u[14]-x[22])/timeconst;
  dx[23] = (u[23]*u[14]-x[23])/timeconst;
  dx[24] = (u[24]*u[14]-x[24])/timeconst;
  dx[25] = (u[25]*u[14]-x[25])/timeconst;
  }

else {
  dx[0] = 0;
  dx[1] = 0;
  dx[2] = 0;
  dx[3] = 0;
  dx[4] = 0;
  dx[5] = 0;
  dx[6] = 0;
  dx[7] = 0;
  dx[8] = 0;
  dx[9] = 0;
  dx[10] = 0;
  dx[11] = 0;
  dx[12] = 0;
  dx[13] = 0; 
  dx[14] = 0; 
  dx[15] = 0; 
  dx[16] = 0; 
  dx[17] = 0; 
  dx[18] = 0; 
  dx[19] = 0; 
  dx[20] = 0; 
  
  dx[21] = 0; 
  dx[22] = 0; 
  dx[23] = 0; 
  dx[24] = 0; 
  dx[25] = 0; 
  
  x[0] = u[0]*u[14];
  x[1] = u[1]*u[14];
  x[2] = u[2]*u[14];
  x[3] = u[3]*u[14];
  x[4] = u[4]*u[14];
  x[5] = u[5]*u[14];
  x[6] = u[6]*u[14];
  x[7] = u[7]*u[14];
  x[8] = u[8]*u[14];
  x[9] = u[9]*u[14];
  x[10] = u[10]*u[14];
  x[11] = u[11]*u[14];
  x[12] = u[12]*u[14];
  x[13] = u[13]*u[14];
  x[14] = u[14];
  x[15] = u[15];
  x[16] = u[16]*u[14];
  x[17] = u[17]*u[14];
  x[18] = u[18]*u[14];
  x[19] = u[19]*u[14];
  x[20] = u[20]*u[14];
  
  x[21] = u[21]*u[14];
  x[22] = u[22]*u[14];
  x[23] = u[23]*u[14];
  x[24] = u[24]*u[14];
  x[25] = u[25]*u[14];
}
}


/*
 * mdlTerminate - called when the simulation is terminated.
 */
static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

