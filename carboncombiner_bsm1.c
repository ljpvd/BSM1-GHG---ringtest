/*
 * carboncombiner_bsm2.c calculates the concentration (Ss) and flow rate when adding an external  
 * carbon source flow rate (with a constant concentration) to the general vector.
 *  
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *
 * The file has been modified in order to include two step nitrification and four step denitrification
 * according to the principles stated in Hiatt and Grady 2008
 *
 * Copyright: Xavier Flores-Alsina, modelEAU, Universite Laval, Quebec, Canada
 *                                  IEA, Lund University, Lund, Sweden 
 */


#define S_FUNCTION_NAME carboncombiner_bsm1

#include "simstruc.h"

#define CARBONSOURCECONC   ssGetArg(S,0)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 27);   /* number of inputs                      */
    ssSetNumOutputs(       S, 26);   /* number of outputs                     */
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
  double csourceconc;
  
  csourceconc = mxGetPr(CARBONSOURCECONC)[0];
    
  if ((u[0] > 0) || (u[15] > 0)) {
    y[0]=(u[1]*u[15])/(u[0]+u[15]); // SI
    y[1]=(u[2]*u[15] + csourceconc*u[0])/(u[0]+u[15]); // SS
    y[2]=(u[3]*u[15])/(u[0]+u[15]); // XI
    y[3]=(u[4]*u[15])/(u[0]+u[15]); // XS
    y[4]=(u[5]*u[15])/(u[0]+u[15]); // XBH
    y[5]=(u[6]*u[15])/(u[0]+u[15]); // XBA1
    y[6]=(u[7]*u[15])/(u[0]+u[15]); // XBA2
    y[7]=(u[8]*u[15])/(u[0]+u[15]); // XU
    y[8]=(u[9]*u[15])/(u[0]+u[15]); // SO
    y[9]=(u[10]*u[15])/(u[0]+u[15]); // SNO3
    y[10]=(u[11]*u[15])/(u[0]+u[15]); // SNH
    y[11]=(u[12]*u[15])/(u[0]+u[15]); // SND
    y[12]=(u[13]*u[15])/(u[0]+u[15]); // XND
    y[13]=(u[14]*u[15])/(u[0]+u[15]); // SALK
    
    y[14]=(u[0]+u[15]); // Q
    y[15]=u[16]; // T
    
    y[16]=(u[17]*u[15])/(u[0]+u[15]); // SNO2
    y[17]=(u[18]*u[15])/(u[0]+u[15]); // SNO
    y[18]=(u[19]*u[15])/(u[0]+u[15]); // SN20
    y[19]=(u[20]*u[15])/(u[0]+u[15]); // SN2
    y[20]=(u[21]*u[15])/(u[0]+u[15]); // XBA2
    
    y[21]=(u[22]*u[15])/(u[0]+u[15]); // D1
    y[22]=(u[23]*u[15])/(u[0]+u[15]); // D2
    y[23]=(u[24]*u[15])/(u[0]+u[15]); // D3
    y[24]=(u[25]*u[15])/(u[0]+u[15]); // D4
    y[25]=(u[26]*u[15])/(u[0]+u[15]); // D5
  }
  else {
    y[0]=0;
    y[1]=0;
    y[2]=0;
    y[3]=0;
    y[4]=0;
    y[5]=0;
    y[6]=0;
    y[7]=0;
    y[8]=0;
    y[9]=0;
    y[10]=0;
    y[11]=0;
    y[12]=0;
    
    y[13]=0;
    y[14]=0;
    
    y[15]=0;
    y[16]=0;
    y[17]=0;
    y[18]=0;
    y[19]=0;
    
    y[20]=0;
    y[21]=0;
    y[22]=0;
    y[23]=0;
    y[24]=0;
    y[25]=0;
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
