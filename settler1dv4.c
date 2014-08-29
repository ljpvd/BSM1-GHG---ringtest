/*
 * SETTLER1D is a C-file S-function for defining a 10 layer settler model.  
 * can simulate 0, 1 or 10 layers for the solubles by using MODELTYPE
 *
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 */

#define S_FUNCTION_NAME settler1dv4

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetArg(S,0)
#define PAR	ssGetArg(S,1)
#define DIM	ssGetArg(S,2)
#define LAYER	ssGetArg(S,3)
#define MODELTYPE	ssGetArg(S,4)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 80);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 17);   /* number of inputs                      */
    ssSetNumOutputs(       S, 113);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 5);   /* number of input arguments             */
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

for (i = 0; i < 80; i++) {
   x0[i] = mxGetPr(XINIT)[i];
}

}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
  double gamma, gamma_eff, modeltype;
  int i;

  gamma = x[9]/u[13];
  gamma_eff = x[0]/u[13];

  modeltype = mxGetPr(MODELTYPE)[0];

  if (modeltype < 0.5) {
     /* underflow */
     y[0]=x[19];
     y[1]=x[29];
     y[2]=u[2]*gamma;
     y[3]=u[3]*gamma;
     y[4]=u[4]*gamma;
     y[5]=u[5]*gamma;
     y[6]=u[6]*gamma;
     y[7]=x[39]; /* use oxygen in return sludge flow */
     y[8]=x[49];
     y[9]=x[59];
     y[10]=x[69];
     y[11]=u[11]*gamma;
     y[12]=x[79];
     y[13]=x[9];
     y[14]=u[15];  /* Q_r */
     y[15]=u[16];  /* Q_w */
  
     /* effluent */
     y[16]=x[10];
     y[17]=x[20];
     y[18]=u[2]*gamma_eff;
     y[19]=u[3]*gamma_eff;
     y[20]=u[4]*gamma_eff;
     y[21]=u[5]*gamma_eff;
     y[22]=u[6]*gamma_eff;
     y[23]=x[30]; /* use oxygen in effluent flow */
     y[24]=x[40];
     y[25]=x[50];
     y[26]=x[60];
     y[27]=u[11]*gamma_eff;
     y[28]=x[70];
     y[29]=x[0];
     y[30]=u[14]-u[15]-u[16];  /* Q_e */

     /* internal TSS states */
     y[31]=x[0];
     y[32]=x[1];
     y[33]=x[2];
     y[34]=x[3];
     y[35]=x[4];
     y[36]=x[5];
     y[37]=x[6];
     y[38]=x[7];
     y[39]=x[8];
     y[40]=x[9];

     y[41]=gamma;
     y[42]=gamma_eff;

     for (i = 10; i < 80; i++) 
        y[i+33] = x[i];
  }

  else if ((modeltype > 0.5) && (modeltype < 1.5)) {
     /* underflow */
     y[0]=x[10];
     y[1]=x[20];
     y[2]=u[2]*gamma;
     y[3]=u[3]*gamma;
     y[4]=u[4]*gamma;
     y[5]=u[5]*gamma;
     y[6]=u[6]*gamma;
     y[7]=x[30]; /* use oxygen in return sludge flow */
     y[8]=x[40];
     y[9]=x[50];
     y[10]=x[60];
     y[11]=u[11]*gamma;
     y[12]=x[70];
     y[13]=x[9];
     y[14]=u[15];  /* Q_r */
     y[15]=u[16];  /* Q_w */
  
     /* effluent */
     y[16]=x[10];
     y[17]=x[20];
     y[18]=u[2]*gamma_eff;
     y[19]=u[3]*gamma_eff;
     y[20]=u[4]*gamma_eff;
     y[21]=u[5]*gamma_eff;
     y[22]=u[6]*gamma_eff;
     y[23]=x[30]; /* use oxygen in effluent flow */
     y[24]=x[40];
     y[25]=x[50];
     y[26]=x[60];
     y[27]=u[11]*gamma_eff;
     y[28]=x[70];
     y[29]=x[0];
     y[30]=u[14]-u[15]-u[16];  /* Q_e */

     /* internal TSS states */
     y[31]=x[0];
     y[32]=x[1];
     y[33]=x[2];
     y[34]=x[3];
     y[35]=x[4];
     y[36]=x[5];
     y[37]=x[6];
     y[38]=x[7];
     y[39]=x[8];
     y[40]=x[9];

     y[41]=gamma;
     y[42]=gamma_eff;

     for (i = 10; i < 20; i++)
        y[i+33] = x[10];
     for (i = 20; i < 30; i++)
        y[i+33] = x[20];
     for (i = 30; i < 40; i++)
        y[i+33] = x[30];
     for (i = 40; i < 50; i++)
        y[i+33] = x[40];
     for (i = 50; i < 60; i++)
        y[i+33] = x[50];
     for (i = 60; i < 70; i++)
        y[i+33] = x[60];
     for (i = 70; i < 80; i++)
        y[i+33] = x[70];
  }

  else if (modeltype > 1.5) {
     /* underflow */
     y[0]=u[0];
     y[1]=u[1];
     y[2]=u[2]*gamma;
     y[3]=u[3]*gamma;
     y[4]=u[4]*gamma;
     y[5]=u[5]*gamma;
     y[6]=u[6]*gamma;
     y[7]=u[7]; /* use oxygen in return sludge flow */
     y[8]=u[8];
     y[9]=u[9];
     y[10]=u[10];
     y[11]=u[11]*gamma;
     y[12]=u[12];
     y[13]=x[9];
     y[14]=u[15];  /* Q_r */
     y[15]=u[16];  /* Q_w */
  
     /* effluent */
     y[16]=u[0];
     y[17]=u[1];
     y[18]=u[2]*gamma_eff;
     y[19]=u[3]*gamma_eff;
     y[20]=u[4]*gamma_eff;
     y[21]=u[5]*gamma_eff;
     y[22]=u[6]*gamma_eff;
     y[23]=u[7]; /* use oxygen in effluent flow */
     y[24]=u[8];
     y[25]=u[9];
     y[26]=u[10];
     y[27]=u[11]*gamma_eff;
     y[28]=u[12];
     y[29]=x[0];
     y[30]=u[14]-u[15]-u[16];  /* Q_e */

     /* internal TSS states */
     y[31]=x[0];
     y[32]=x[1];
     y[33]=x[2];
     y[34]=x[3];
     y[35]=x[4];
     y[36]=x[5];
     y[37]=x[6];
     y[38]=x[7];
     y[39]=x[8];
     y[40]=x[9];

     y[41]=gamma;
     y[42]=gamma_eff;

     for (i = 10; i < 20; i++)
        y[i+33] = u[0];
     for (i = 20; i < 30; i++)
        y[i+33] = u[1];
     for (i = 30; i < 40; i++)
        y[i+33] = u[7];
     for (i = 40; i < 50; i++)
        y[i+33] = u[8];
     for (i = 50; i < 60; i++)
        y[i+33] = u[9];
     for (i = 60; i < 70; i++)
        y[i+33] = u[10];
     for (i = 70; i < 80; i++)
        y[i+33] = u[12];
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

double v0_max, v0, r_h, r_p, f_ns, X_t, area, h, feedlayer, volume, modeltype;
double Q_f, Q_e, Q_u, v_up, v_dn, v_in, eps;
int i;
double vs[10];
double Js[11];
double Jstemp[10];
double Jflow[11];

v0_max = mxGetPr(PAR)[0];
v0 = mxGetPr(PAR)[1];
r_h = mxGetPr(PAR)[2];
r_p = mxGetPr(PAR)[3];
f_ns = mxGetPr(PAR)[4];
X_t = mxGetPr(PAR)[5];
area = mxGetPr(DIM)[0];
h = mxGetPr(DIM)[1]/mxGetPr(LAYER)[1];
feedlayer = mxGetPr(LAYER)[0];
modeltype = mxGetPr(MODELTYPE)[0];
volume = area*mxGetPr(DIM)[1];

eps = 0.01;
v_in = u[14]/area;
Q_f = u[14];
Q_u = u[15] + u[16];
Q_e = u[14] - Q_u;
v_up = Q_e/area;
v_dn = Q_u/area;

for (i = 0; i < 10; i++) {
   vs[i] = v0*(exp(-r_h*(x[i]-f_ns*u[13]))-exp(-r_p*(x[i]-f_ns*u[13])));
   if (vs[i] > v0_max)     
      vs[i] = v0_max;
   else if (vs[i] < 0)
      vs[i] = 0;
}

for (i = 0; i < 10; i++) {
   Jstemp[i] = vs[i]*x[i];
}

// bulk movement
for (i = 0; i < 11; i++) {
   if (i < (feedlayer-eps))     
      Jflow[i] = v_up*x[i];
   else
      Jflow[i] = v_dn*x[i-1];
}


Js[0] = 0;
Js[10] = 0;
for (i = 0; i < 9; i++) {
   if ((i < (feedlayer-1-eps)) && (x[i+1] <= X_t))
      Js[i+1] = Jstemp[i];
   else if (Jstemp[i] < Jstemp[i+1])     
      Js[i+1] = Jstemp[i];
   else
      Js[i+1] = Jstemp[i+1];
}


for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i] = (-Jflow[i]+Jflow[i+1]+Js[i]-Js[i+1])/h;
   else if (i > (feedlayer-eps))
      dx[i] = (Jflow[i]-Jflow[i+1]+Js[i]-Js[i+1])/h;
   else
      dx[i] = (v_in*u[13]-Jflow[i]-Jflow[i+1]+Js[i]-Js[i+1])/h;
}

/* soluble component S_I */
if (modeltype < 0.5) {
   for (i = 10; i < 20; i++) {
      if (i < (feedlayer-1+10-eps))
         dx[i] = (-v_up*x[i]+v_up*x[i+1])/h;
      else if (i > (feedlayer+10-eps))
         dx[i] = (v_dn*x[i-1]-v_dn*x[i])/h;
      else
         dx[i] = (v_in*u[0]-v_up*x[i]-v_dn*x[i])/h;
   }
}
else if ((modeltype > 0.5) && (modeltype < 1.5)) {
   dx[10] = (Q_f*(u[0]-x[10]))/volume;
   for (i = 11; i < 20; i++) 
      dx[i] = 0;
}
else if (modeltype > 1.5) {
   for (i = 10; i < 20; i++) 
      dx[i] = 0;
}


/* soluble component S_S */
if (modeltype < 0.5) {
   for (i = 20; i < 30; i++) {
      if (i < (feedlayer-1+20-eps))
         dx[i] = (-v_up*x[i]+v_up*x[i+1])/h;
      else if (i > (feedlayer+20-eps))
         dx[i] = (v_dn*x[i-1]-v_dn*x[i])/h;
      else
         dx[i] = (v_in*u[1]-v_up*x[i]-v_dn*x[i])/h;
   }
}
else if ((modeltype > 0.5) && (modeltype < 1.5)) {
   dx[20] = (Q_f*(u[1]-x[20]))/volume;
   for (i = 21; i < 30; i++) 
      dx[i] = 0;
}
else if (modeltype > 1.5) {
   for (i = 20; i < 30; i++) 
      dx[i] = 0;
}

/* soluble component S_O */
if (modeltype < 0.5) {
   for (i = 30; i < 40; i++) {
      if (i < (feedlayer-1+30-eps))
         dx[i] = (-v_up*x[i]+v_up*x[i+1])/h;
      else if (i > (feedlayer+30-eps))
         dx[i] = (v_dn*x[i-1]-v_dn*x[i])/h;
      else
         dx[i] = (v_in*u[7]-v_up*x[i]-v_dn*x[i])/h;
   }
}
else if ((modeltype > 0.5) && (modeltype < 1.5)) {
   dx[30] = (Q_f*(u[7]-x[30]))/volume;
   for (i = 31; i < 40; i++) 
      dx[i] = 0;
}
else if (modeltype > 1.5) {
   for (i = 30; i < 40; i++) 
      dx[i] = 0;
}

/* soluble component S_NO */
if (modeltype < 0.5) {
   for (i = 40; i < 50; i++) {
      if (i < (feedlayer-1+40-eps))
         dx[i] = (-v_up*x[i]+v_up*x[i+1])/h;
      else if (i > (feedlayer+40-eps))
         dx[i] = (v_dn*x[i-1]-v_dn*x[i])/h;
      else
         dx[i] = (v_in*u[8]-v_up*x[i]-v_dn*x[i])/h;
   }
}
else if ((modeltype > 0.5) && (modeltype < 1.5)) {
   dx[40] = (Q_f*(u[8]-x[40]))/volume;
   for (i = 41; i < 50; i++) 
      dx[i] = 0;
}
else if (modeltype > 1.5) {
   for (i = 40; i < 50; i++) 
      dx[i] = 0;
}

/* soluble component S_NH */
if (modeltype < 0.5) {
   for (i = 50; i < 60; i++) {
      if (i < (feedlayer-1+50-eps))
         dx[i] = (-v_up*x[i]+v_up*x[i+1])/h;
      else if (i > (feedlayer+50-eps))
         dx[i] = (v_dn*x[i-1]-v_dn*x[i])/h;
      else
         dx[i] = (v_in*u[9]-v_up*x[i]-v_dn*x[i])/h;
   }
}
else if ((modeltype > 0.5) && (modeltype < 1.5)) {
   dx[50] = (Q_f*(u[9]-x[50]))/volume;
   for (i = 51; i < 60; i++) 
      dx[i] = 0;
}
else if (modeltype > 1.5) {
   for (i = 50; i < 60; i++) 
      dx[i] = 0;
}

/* soluble component S_ND */
if (modeltype < 0.5) {
   for (i = 60; i < 70; i++) {
      if (i < (feedlayer-1+60-eps))
         dx[i] = (-v_up*x[i]+v_up*x[i+1])/h;
      else if (i > (feedlayer+60-eps))
         dx[i] = (v_dn*x[i-1]-v_dn*x[i])/h;
      else
         dx[i] = (v_in*u[10]-v_up*x[i]-v_dn*x[i])/h;
   }
}
else if ((modeltype > 0.5) && (modeltype < 1.5)) {
   dx[60] = (Q_f*(u[10]-x[60]))/volume;
   for (i = 61; i < 70; i++) 
      dx[i] = 0;
}
else if (modeltype > 1.5) {
   for (i = 60; i < 70; i++) 
      dx[i] = 0;
}

/* soluble component S_ALK */
if (modeltype < 0.5) {
   for (i = 70; i < 80; i++) {
      if (i < (feedlayer-1+70-eps))
         dx[i] = (-v_up*x[i]+v_up*x[i+1])/h;
      else if (i > (feedlayer+70-eps))
         dx[i] = (v_dn*x[i-1]-v_dn*x[i])/h;
      else
         dx[i] = (v_in*u[12]-v_up*x[i]-v_dn*x[i])/h;
   }
}
else if ((modeltype > 0.5) && (modeltype < 1.5)) {
   dx[70] = (Q_f*(u[12]-x[70]))/volume;
   for (i = 71; i < 80; i++) 
      dx[i] = 0;
}
else if (modeltype > 1.5) {
   for (i = 70; i < 80; i++) 
      dx[i] = 0;
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
