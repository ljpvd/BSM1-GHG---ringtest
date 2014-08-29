/*
 * settler1d_reac is a C-file S-function for defining a reactive 10 layer sedimentation tank model
 * Compatibility: ASM1 model
 * The model is developed based on the settler1dv4.c model of Ulf Jeppsson (IEA, LTH, Sweden)
 *
 * 25 October 2003
 * Krist V. Gernaey, CAPEC, Dept. of Chemical Engineering, DTU, Denmark
 *
 * The file has been modified in order to include two step nitrification and four step denitrification
 * according to the principles stated in Hiatt et al., 2008
 *
 * Copyright: Xavier Flores-Alsina, modelEAU, Universite Laval, Quebec, Canada
 *                                  IEA, Lund University, Lund, Sweden 

 */

#define S_FUNCTION_NAME settler1d_bsm1

#include "simstruc.h"
#include <math.h>



#define XINIT   ssGetArg(S,0)   /* initial values                    */
#define SEDPAR	 ssGetArg(S,1)  /* parameters sedimentation model    */
#define DIM	    ssGetArg(S,2)   /* dimensions clarifier              */
#define LAYER   ssGetArg(S,3)   /* Number of layers                  */
#define TEMPMODEL  ssGetArg(S,4)
#define ACTIVATE  ssGetArg(S,5)  /* reactive settler*/



/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(  S, 250);   /* number of continuous states */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(       S, 28);   /* number of inputs */
    ssSetNumOutputs(     S, 253);   /* number of outputs: 140 states (13 ASM1 components + TSS for each layer) + effluent flow rate + return sludge flow rate + waste sludge flow rate */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 6);   /* number of input arguments             */
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

for (i = 0; i < 250; i++) {
   x0[i] = mxGetPr(XINIT)[i];
}

}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
double tempmodel;
double activate;

int i,j,k,z;

tempmodel  = mxGetPr(TEMPMODEL)[0];
activate =  mxGetPr(ACTIVATE)[0];



    /* Clarifier state variables */
	for (i=0;i<10;i++){
        
		for (j=0;j<15;j++){
		y[(i*25)+j]=x[i+j*10];
		}
        
	    if (tempmodel < 0.5) {  
        y[(i*25)+15]=u[14];
        }
        else {
        y[(i*25)+15]=x[i+14*10];
        }
        
        for (z=15;z<20;z++){
		y[(i*25)+z]=x[i+z*10];
		}
        
        if (activate < 0.5) {
        for (k=20;k<25;k++){
		y[(i*25)+k]=0 ;}
        }
        else {
        for (k=20;k<25;k++){   
        y[(i*25)+k]=x[i+k*10];
        }
        }           
            
    }


    /* Flow rates out of the clarifier */
   y[250]=u[14]-u[27]-u[26];     /* Q_e  */
   y[251]=u[26];                 /* Q_r  */
   y[252]=u[27];                 /* Q_w  */
  
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

double v0_max, v0, r_h, r_p, f_ns, X_t, area, h, feedlayer, volume;
double Q_f, Q_e, Q_u, v_up, v_dn, v_in, eps;
int i;

double vs[10];
double Js[11];
double Jstemp[10];
double Jflow[11];

double reac1[10], reac2[10], reac3[10], reac4[10], reac5[10], reac6[10], reac7[10], reac8[10];
double reac9[10], reac10[10], reac11[10], reac12[10], reac13[10], reac14[10];
double reac16[10], reac17[10], reac18[10], reac19[10], reac20[10];
double reac21[10], reac22[10], reac23[10], reac24[10], reac25[10];

double xtemp[250];

v0_max = mxGetPr(SEDPAR)[0];
v0 = mxGetPr(SEDPAR)[1];
r_h = mxGetPr(SEDPAR)[2];
r_p = mxGetPr(SEDPAR)[3];
f_ns = mxGetPr(SEDPAR)[4];
X_t = mxGetPr(SEDPAR)[5];
area = mxGetPr(DIM)[0];
h = mxGetPr(DIM)[1]/mxGetPr(LAYER)[1];
feedlayer = mxGetPr(LAYER)[0];
volume = area*mxGetPr(DIM)[1];

eps = 0.01;
v_in = u[14]/area;
Q_f = u[14];
Q_u = u[26] + u[27];  /* underflow flow rate */
Q_e = u[14] - Q_u;
v_up = Q_e/area;
v_dn = Q_u/area;


for (i = 0; i < 250; i++) {
   if (x[i] < 0.0)
     xtemp[i] = 0.0;
   else
     xtemp[i] = x[i];
}

/* Reaction rates ASM1 model */

for (i = 0; i < 10; i++) {
   
reac1[i] = 0.0;
reac2[i] = 0.0;
reac3[i] = 0.0;
reac4[i] = 0.0;
reac5[i] = 0.0;
reac6[i] = 0.0;
reac7[i] = 0.0;
reac8[i] = 0.0;
reac9[i] = 0.0;
reac10[i] = 0.0;
reac11[i] = 0.0;
reac12[i] = 0.0;
reac13[i] = 0.0;
reac14[i] = 0.0;

reac16[i] = 0.0;
reac17[i] = 0.0;
reac18[i] = 0.0;
reac19[i] = 0.0;
reac20[i] = 0.0;

reac21[i]  = 0.0;
reac22[i]  = 0.0;
reac23[i]  = 0.0;
reac24[i]  = 0.0;
reac25[i]  = 0.0;


}

/* calculation of the sedimentation velocity for each of the layers */
for (i = 0; i < 10; i++) {
   vs[i] = v0*(exp(-r_h*(xtemp[i+130]-f_ns*u[13]))-exp(-r_p*(xtemp[i+130]-f_ns*u[13])));      /* u[13] = influent SS concentration */
   if (vs[i] > v0_max)     
      vs[i] = v0_max;
   else if (vs[i] < 0.0)
      vs[i] = 0.0;
}

/* calculation of the sludge flux due to sedimentation for each layer (not taking into account X limit) */
for (i = 0; i < 10; i++) {
   Jstemp[i] = vs[i]*xtemp[i+130];
}

/* calculation of the sludge flux due to the liquid flow (upflow or downflow, depending on layer) */
for (i = 0; i < 11; i++) {
   if (i < (feedlayer-eps))     
      Jflow[i] = v_up*xtemp[i+130];
   else
      Jflow[i] = v_dn*xtemp[i-1+130];
}

/* calculation of the sludge flux due to sedimentation for each layer */
Js[0] = 0.0;
Js[10] = 0.0;
for (i = 0; i < 9; i++) {
   if ((i < (feedlayer-1-eps)) && (xtemp[i+1+130] <= X_t))
      Js[i+1] = Jstemp[i];
   else if (Jstemp[i] < Jstemp[i+1])     
      Js[i+1] = Jstemp[i];
   else
      Js[i+1] = Jstemp[i+1];
}

/* Reaction rates ASM1 model */



/* ASM1 model component balances over the layers */
/* ASM1: [Si Ss Xi Xs Xbh Xba Xp So Sno Snh Snd Xnd Salk TSS Q_in]   */

/* soluble component S_I */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i]  = (-v_up*xtemp[i] +v_up*xtemp[i+1])/h+reac1[i];                              
   else if (i > (feedlayer-eps)) 
      dx[i]  = (v_dn*xtemp[i-1] -v_dn*xtemp[i])/h+reac1[i];
   else
      dx[i]  = (v_in*u[0] -v_up*xtemp[i] -v_dn*xtemp[i])/h+reac1[i];
}

/* soluble component S_S */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+10]  = (-v_up*xtemp[i+10] +v_up*xtemp[i+1+10])/h +reac2[i];
   else if (i > (feedlayer-eps)) 
      dx[i+10]  = (v_dn*xtemp[i-1+10] -v_dn*xtemp[i+10])/h +reac2[i];
   else
      dx[i+10]  = (v_in*u[1] -v_up*xtemp[i+10] -v_dn*xtemp[i+10])/h +reac2[i];
}

/* particulate component X_I */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+20] = ((xtemp[i+20]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+20]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+20]/xtemp[i+1+130])*Jflow[i+1])/h +reac3[i];
   else if (i > (feedlayer-eps)) 
      dx[i+20] = ((xtemp[i+20]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+20]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac3[i];
   else
      dx[i+20] = ((xtemp[i+20]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+20]/xtemp[i-1+130])*Js[i]+v_in*u[2])/h +reac3[i];
}

/* particulate component X_S */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
       dx[i+30] = ((xtemp[i+30]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+30]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+30]/xtemp[i+1+130])*Jflow[i+1])/h +reac4[i];
   else if (i > (feedlayer-eps)) 
       dx[i+30] = ((xtemp[i+30]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+30]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac4[i];
   else
       dx[i+30] = ((xtemp[i+30]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+30]/xtemp[i-1+130])*Js[i]+v_in*u[3])/h +reac4[i];
}

/* particulate component X_BH */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
       dx[i+40] = ((xtemp[i+40]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+40]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+40]/xtemp[i+1+130])*Jflow[i+1])/h +reac5[i];
   else if (i > (feedlayer-eps)) 
       dx[i+40] = ((xtemp[i+40]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+40]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac5[i];
   else
       dx[i+40] = ((xtemp[i+40]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+40]/xtemp[i-1+130])*Js[i]+v_in*u[4])/h +reac5[i];
}

/* particulate component X_BA */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
       dx[i+50] = ((xtemp[i+50]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+50]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+50]/xtemp[i+1+130])*Jflow[i+1])/h +reac6[i];
   else if (i > (feedlayer-eps)) 
       dx[i+50] = ((xtemp[i+50]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+50]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac6[i];
   else
       dx[i+50] = ((xtemp[i+50]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+50]/xtemp[i-1+130])*Js[i]+v_in*u[5])/h +reac6[i];
}

/* particulate component X_P */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
       dx[i+60] = ((xtemp[i+60]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+60]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+60]/xtemp[i+1+130])*Jflow[i+1])/h +reac7[i];
   else if (i > (feedlayer-eps)) 
       dx[i+60] = ((xtemp[i+60]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+60]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac7[i];
   else
       dx[i+60] = ((xtemp[i+60]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+60]/xtemp[i-1+130])*Js[i]+v_in*u[6])/h +reac7[i];
}

/* soluble component S_O */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+70]  = (-v_up*xtemp[i+70] +v_up*xtemp[i+1+70])/h +reac8[i];
   else if (i > (feedlayer-eps)) 
      dx[i+70]  = (v_dn*xtemp[i-1+70] -v_dn*xtemp[i+70])/h +reac8[i];
   else
      dx[i+70]  = (v_in*u[7] -v_up*xtemp[i+70] -v_dn*xtemp[i+70])/h +reac8[i];
}

/* soluble component S_NO */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+80]  = (-v_up*xtemp[i+80] +v_up*xtemp[i+1+80])/h +reac9[i];
   else if (i > (feedlayer-eps)) 
      dx[i+80]  = (v_dn*xtemp[i-1+80] -v_dn*xtemp[i+80])/h +reac9[i];
   else
      dx[i+80]  = (v_in*u[8] -v_up*xtemp[i+80] -v_dn*xtemp[i+80])/h +reac9[i];
}

/* soluble component S_NH */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+90]  = (-v_up*xtemp[i+90] +v_up*xtemp[i+1+90])/h +reac10[i];
   else if (i > (feedlayer-eps)) 
      dx[i+90]  = (v_dn*xtemp[i-1+90] -v_dn*xtemp[i+90])/h +reac10[i];
   else
      dx[i+90]  = (v_in*u[9] -v_up*xtemp[i+90] -v_dn*xtemp[i+90])/h +reac10[i];
}

/* soluble component S_ND */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+100]  = (-v_up*xtemp[i+100] +v_up*xtemp[i+1+100])/h +reac11[i];
   else if (i > (feedlayer-eps)) 
      dx[i+100]  = (v_dn*xtemp[i-1+100] -v_dn*xtemp[i+100])/h +reac11[i];
   else
      dx[i+100]  = (v_in*u[10] -v_up*xtemp[i+100] -v_dn*xtemp[i+100])/h +reac11[i];
}

/* particulate component X_ND */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+110] = ((xtemp[i+110]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+110]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+110]/xtemp[i+1+130])*Jflow[i+1])/h +reac12[i];
   else if (i > (feedlayer-eps)) 
      dx[i+110] = ((xtemp[i+110]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+110]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac12[i];
   else
      dx[i+110] = ((xtemp[i+110]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+110]/xtemp[i-1+130])*Js[i]+v_in*u[11])/h +reac12[i];
}

/* soluble component S_ALK */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+120]  = (-v_up*xtemp[i+120] +v_up*xtemp[i+1+120])/h +reac13[i];
   else if (i > (feedlayer-eps)) 
      dx[i+120]  = (v_dn*xtemp[i-1+120] -v_dn*xtemp[i+120])/h +reac13[i];
   else
      dx[i+120]  = (v_in*u[12] -v_up*xtemp[i+120] -v_dn*xtemp[i+120])/h +reac13[i];
}

/* particulate component X_TSS */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+130] = ((-Jflow[i]-Js[i+1])+Js[i]+Jflow[i+1])/h +reac14[i];
   else if (i > (feedlayer-eps)) 
      dx[i+130] = ((-Jflow[i+1]-Js[i+1])+(Jflow[i]+Js[i]))/h +reac14[i];
   else
      dx[i+130] = ((-Jflow[i]-Jflow[i+1]-Js[i+1])+Js[i]+v_in*u[13])/h +reac14[i];
}

/* Sol. comp. T     ; ASM1mod */

for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+140]  = (-v_up*xtemp[i+140] +v_up*xtemp[i+1+140])/h ;
   else if (i > (feedlayer-eps)) 
      dx[i+140]  = (v_dn*xtemp[i-1+140] -v_dn*xtemp[i+140])/h ;
   else
      dx[i+140]  = (v_in*u[15] -v_up*xtemp[i+140] -v_dn*xtemp[i+140])/h ;
}


/* Sol. comp. SNO2     ; ASM1mod */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+150]  = (-v_up*xtemp[i+150] +v_up*xtemp[i+1+150])/h +reac16[i];
   else if (i > (feedlayer-eps)) 
      dx[i+150]  = (v_dn*xtemp[i-1+150] -v_dn*xtemp[i+150])/h +reac16[i];
   else
      dx[i+150]  = (v_in*u[16] -v_up*xtemp[i+150] -v_dn*xtemp[i+150])/h +reac16[i];
}

/* Sol. comp. SNO     ; ASM1mod */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+160]  = (-v_up*xtemp[i+160] +v_up*xtemp[i+1+160])/h +reac17[i];
   else if (i > (feedlayer-eps)) 
      dx[i+160]  = (v_dn*xtemp[i-1+160] -v_dn*xtemp[i+160])/h +reac17[i];
   else
      dx[i+160]  = (v_in*u[17] -v_up*xtemp[i+160] -v_dn*xtemp[i+160])/h +reac17[i];
}
   
/* Sol. comp. SN2O     ; ASM1mod */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+170]  = (-v_up*xtemp[i+170] +v_up*xtemp[i+1+170])/h +reac18[i];
   else if (i > (feedlayer-eps)) 
      dx[i+170]  = (v_dn*xtemp[i-1+170] -v_dn*xtemp[i+170])/h +reac18[i];
   else
      dx[i+170]  = (v_in*u[18] -v_up*xtemp[i+170] -v_dn*xtemp[i+170])/h +reac18[i];
}   

/* Sol. comp. SN2     ; ASM1mod */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps))
      dx[i+180]  = (-v_up*xtemp[i+180] +v_up*xtemp[i+1+180])/h +reac19[i];
   else if (i > (feedlayer-eps)) 
      dx[i+180]  = (v_dn*xtemp[i-1+180] -v_dn*xtemp[i+180])/h +reac19[i];
   else
      dx[i+180]  = (v_in*u[19] -v_up*xtemp[i+180] -v_dn*xtemp[i+180])/h +reac19[i];
}   


 /* Part. comp. XBA2    ; ASM1mode */
      dx[190] = ((x[190]/x[130])*(-Jflow[0]-Js[1])+(x[191]/x[131])*Jflow[1])/h +reac20[0];
for (i = 1; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+190] = ((xtemp[i+190]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+190]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+190]/xtemp[i+1+130])*Jflow[i+1])/h +reac20[i];
   else if (i > (feedlayer-eps)) 
      dx[i+190] = ((xtemp[i+190]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+190]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac20[i];
   else
      dx[i+190] = ((xtemp[i+190]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+190]/xtemp[i-1+130])*Js[i]+v_in*u[20])/h +reac20[i]; 
}  


/* Sol. comp. D1   ; ASM1mod */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+200]  = (-v_up*xtemp[i+200] +v_up*xtemp[i+1+200])/h + reac21[i];
   else if (i > (feedlayer-eps)) 
      dx[i+200]  = (v_dn*xtemp[i-1+200] -v_dn*xtemp[i+200])/h +  reac21[i];
   else
      dx[i+200]  = (v_in*u[21] -v_up*xtemp[i+200] -v_dn*xtemp[i+200])/h + + reac21[i];
}  
            
/* Sol. comp. D2   ; ASM1mod */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+210]  = (-v_up*xtemp[i+210] +v_up*xtemp[i+1+210])/h + reac22[i];
   else if (i > (feedlayer-eps)) 
      dx[i+210]  = (v_dn*xtemp[i-1+210] -v_dn*xtemp[i+210])/h +  reac22[i];
   else
      dx[i+210]  = (v_in*u[22] -v_up*xtemp[i+210] -v_dn*xtemp[i+210])/h +  reac22[i];
}   
      
/* Sol. comp. D3   ; ASM1mod */
for (i = 0; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+220]  = (-v_up*xtemp[i+220] +v_up*xtemp[i+1+220])/h + reac23[i];
   else if (i > (feedlayer-eps)) 
      dx[i+220]  = (v_dn*xtemp[i-1+220] -v_dn*xtemp[i+220])/h +  reac23[i];
   else
      dx[i+220]  = (v_in*u[23] -v_up*xtemp[i+220] -v_dn*xtemp[i+220])/h +  reac23[i];
}    
      
      
/* Part. comp. D4    ; ASM1mod */
      dx[230] = ((x[230]/x[130])*(-Jflow[0]-Js[1])+(x[231]/x[131])*Jflow[1])/h +reac24[0];
for (i = 1; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+230] = ((xtemp[i+230]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+230]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+230]/xtemp[i+1+130])*Jflow[i+1])/h +reac24[i];
   else if (i > (feedlayer-eps)) 
      dx[i+230] = ((xtemp[i+230]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+230]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac24[i];
   else
      dx[i+230] = ((xtemp[i+230]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+230]/xtemp[i-1+130])*Js[i]+v_in*u[23])/h +reac24[i];
}  

 /* Part. comp. D5    ; ASM1mod */
      dx[240] = ((x[240]/x[130])*(-Jflow[0]-Js[1])+(x[241]/x[131])*Jflow[1])/h +reac25[0];
for (i = 1; i < 10; i++) {
   if (i < (feedlayer-1-eps)) 
      dx[i+240] = ((xtemp[i+240]/xtemp[i+130])*(-Jflow[i]-Js[i+1])+(xtemp[i-1+240]/xtemp[i-1+130])*Js[i]+(xtemp[i+1+240]/xtemp[i+1+130])*Jflow[i+1])/h +reac25[i];
   else if (i > (feedlayer-eps)) 
      dx[i+240] = ((xtemp[i+240]/xtemp[i+130])*(-Jflow[i+1]-Js[i+1])+(xtemp[i-1+240]/xtemp[i-1+130])*(Jflow[i]+Js[i]))/h +reac25[i];
   else
      dx[i+240] = ((xtemp[i+240]/xtemp[i+130])*(-Jflow[i]-Jflow[i+1]-Js[i+1])+(xtemp[i-1+240]/xtemp[i-1+130])*Js[i]+v_in*u[24])/h +reac25[i]; 
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
