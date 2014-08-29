/*
 * flowsplitter_bsm2.c calculates the output concentrations when splitting two flow  
 * streams. Output temperature is always identical in both flows,
 * i.e. parameter TEMPMODEL is not used. If any output flow rate is
 * (less or) equal to zero then all outputs are set to zero.
 * TYPE = 0: specific fraction of influent flow rate to first output (value between 0 and 1)
 * TYPE = 1: specific flow rate (m3/d) to first output
 * TYPE = 2: flow rate (m3/d) above specific limit to first output
 * If more flow than available is trying to be diverted then the flow is 
 * automatically limited to the maximum available flow rate.
 * Input(22)=u[21] is the requested bypass flow rate (control input).
 * 
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *
 * The file has been modified in order to include two step nitrification and four step denitrification
 * according to the principles stated in Hiatt and Grady 2008
 *
 * Copyright: Xavier Flores-Alsina, modelEAU, Universite Laval, Quebec, Canada
 *                                  IEA, Lund University, Lund, Sweden
 */

#define S_FUNCTION_NAME flowsplitter_bsm1

#include "simstruc.h"

#define TYPE   ssGetArg(S,0)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 27);  /* number of inputs                      */
    ssSetNumOutputs(       S, 52);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 1);   /* number of input arguments             */
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
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
    
  int type, i;

  type = mxGetPr(TYPE)[0];
  
  for (i = 0; i < 26; i++) {
      y[i] = u[i];
      y[i+26] = u[i];
  }
  
  if (type == 0) {          /* divert specific fraction of influent flow */
    y[14]=u[14]*u[26];
  }
  else if (type == 1) {     /* divert specific flow rate of influent flow */
    if (u[26] <= u[14]) {
      y[14]=u[26];
      }
    else {
      y[14]=u[14];
      }
    }
  else if (type == 2) {     /* divert specific flow rate ABOVE limit value */
    if (u[26] <= u[14]) {
      y[14]=u[14]-u[26];
      }
    else {
      y[14]=0.0;
      }
    }

  y[40]=u[14]-y[14];
      
  if (y[14] <= 0.00001) {
      for (i = 0; i < 26; i++) {
      y[i] = 0.0;
      }
  }

  if (y[40] <= 0.00001) {
      for (i = 0; i < 26; i++) {
      y[i+26] = 0.0;
      }
  }
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

