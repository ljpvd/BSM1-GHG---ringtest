/*
 * combiner_bsm2.c calculates the output concentrations when adding two flow  
 * streams together. Output temperature always based on 'heat content' of the
 * influent flows, i.e. parameter TEMPMODEL is not used. If all input flow rates are
 * less or equal to zero then all outputs are zero.
 *
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *
 * The file has been modified in order to include two step nitrification and four step denitrification
 * according to the principles stated in Hiatt and Grady 2008
 *
 * Copyright: Xavier Flores-Alsina, modelEAU, Universite Laval, Quebec, Canada
 *                                  IEA, Lund University, Lund, Sweden
 */

#define S_FUNCTION_NAME combiner_bsm1


#include "simstruc.h"

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 52);  /* number of inputs                      */
    ssSetNumOutputs(       S, 26);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 0);   /* number of input arguments             */
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
    
  if ((u[14] > 0) || (u[40] > 0)) {
    y[0]=(u[0]*u[14] + u[26]*u[40])/(u[14]+u[40]); // SI
    y[1]=(u[1]*u[14] + u[27]*u[40])/(u[14]+u[40]); // SS
    y[2]=(u[2]*u[14] + u[28]*u[40])/(u[14]+u[40]); // XI
    y[3]=(u[3]*u[14] + u[29]*u[40])/(u[14]+u[40]); // XS
    y[4]=(u[4]*u[14] + u[30]*u[40])/(u[14]+u[40]); // XBH
    y[5]=(u[5]*u[14] + u[31]*u[40])/(u[14]+u[40]); // XBA1
    y[6]=(u[6]*u[14] + u[32]*u[40])/(u[14]+u[40]); // XU
    y[7]=(u[7]*u[14] + u[33]*u[40])/(u[14]+u[40]); // SO
    y[8]=(u[8]*u[14] + u[34]*u[40])/(u[14]+u[40]); // SNO3
    y[9]=(u[9]*u[14] + u[35]*u[40])/(u[14]+u[40]); // SNH
    y[10]=(u[10]*u[14] + u[36]*u[40])/(u[14]+u[40]); // SND
    y[11]=(u[11]*u[14] + u[37]*u[40])/(u[14]+u[40]); // XND
    y[12]=(u[12]*u[14] + u[38]*u[40])/(u[14]+u[40]); // SALK
    y[13]=(u[13]*u[14] + u[39]*u[40])/(u[14]+u[40]); // XTSS
    
    y[14]=u[14]+u[40];                               /* Flow rate */
    
    y[15]=(u[15]*u[14] + u[41]*u[40])/(u[14]+u[40]);    /* Temp */
         
    y[16]=(u[16]*u[14] + u[42]*u[40])/(u[14]+u[40]); // SNO2
    y[17]=(u[17]*u[14] + u[43]*u[40])/(u[14]+u[40]); // SNO
    y[18]=(u[18]*u[14] + u[44]*u[40])/(u[14]+u[40]); // SN2O
    y[19]=(u[19]*u[14] + u[45]*u[40])/(u[14]+u[40]); // SN2
    y[20]=(u[20]*u[14] + u[46]*u[40])/(u[14]+u[40]); // XBA2
    
    /* Dummy states */
    
    y[21]=(u[21]*u[14] + u[47]*u[40])/(u[14]+u[40]); // D1
    y[22]=(u[22]*u[14] + u[48]*u[40])/(u[14]+u[40]); // D2
    y[23]=(u[23]*u[14] + u[49]*u[40])/(u[14]+u[40]); // D3
    y[24]=(u[24]*u[14] + u[50]*u[40])/(u[14]+u[40]); // D4
    y[25]=(u[25]*u[14] + u[51]*u[40])/(u[14]+u[40]); // D5
    
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

