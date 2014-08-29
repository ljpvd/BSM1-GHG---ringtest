/*
 * asm1_bsm2 is a C-file S-function for IAWQ AS Model No 1 with temperature 
 * dependencies of the kinetic parameters. In addition to the ASM1 states, TSS
 * and dummy states are included. TEMPMODEL defines how temperature changes
 * in the input affects the reactor temperature. Temperature dependency for 
 * oxygen saturation concentration and KLa has also been added in accordance
 * with BSM2 documentation.
 *
 * Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
 *
 * The file has been modified in order to include two step nitrification and four step denitrification
 * according to the principles stated in Hiatt et al., 2008
 *
 * Copyright: Xavier Flores-Alsina, modelEAU, Universite Laval, Quebec, Canada
 *                                  IEA, Lund University, Lund, Sweden 
 */

#define S_FUNCTION_NAME asm1_bsm1

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetArg(S,0)
#define PAR	ssGetArg(S,1)
#define V	ssGetArg(S,2)
#define SOSAT	ssGetArg(S,3)
#define TEMPMODEL  ssGetArg(S,4)
#define ACTIVATE  ssGetArg(S,5)

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 26);  /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 27);  /* number of inputs                      */
    ssSetNumOutputs(       S, 29);  /* number of outputs                     */
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
  double tempmodel, activate;
  
  double D_N2, D_N2O , D_NO, D_O2, H_N2, H_N2O, H_NO, P_N2O_air, P_N2_air, P_NO_air;
  double Temp_Ref, theta_kla;
  double Kla_N2O, Kla_NO ,Kla_N2;
  double Flux_NO, Flux_N2O, Flux_N2;
  
  double KLa_temp, vol;
  double S_FA, KB_2_KW, K_A, S_FNA, pH;

  int i;

X_I2TSS =           mxGetPr(PAR)[57];
X_S2TSS =           mxGetPr(PAR)[58];
X_BH2TSS =          mxGetPr(PAR)[59];
X_BA2TSS =          mxGetPr(PAR)[60];
X_P2TSS =           mxGetPr(PAR)[61];

/////////////////////////////////////////////////////
D_N2	=           mxGetPr(PAR)[3];
D_N2O	=           mxGetPr(PAR)[4];
D_NO	=           mxGetPr(PAR)[5];
D_O2    =           mxGetPr(PAR)[6];
H_N2	=           mxGetPr(PAR)[10];
H_N2O	=           mxGetPr(PAR)[11];
H_NO	=           mxGetPr(PAR)[12];
P_N2O_air =     	mxGetPr(PAR)[51];
P_N2_air  =         mxGetPr(PAR)[52];
P_NO_air  =       	mxGetPr(PAR)[53];

Temp_Ref	=       mxGetPr(PAR)[68];
theta_kla	=       mxGetPr(PAR)[72];

pH	    =           mxGetPr(PAR)[50];

vol =               mxGetPr(V)[0];

////////////////////////////////////////////////////////

  tempmodel = mxGetPr(TEMPMODEL)[0];
  activate = mxGetPr(ACTIVATE)[0];
  
  for (i = 0; i < 13; i++) {                    // ASM state variables: SI, SS, XI, XS, XBH, XBA1, XU, SO, SNO3, SNH, SND, XND, SALK
      y[i] = x[i];
  }

  y[13] = X_I2TSS*x[2]+X_S2TSS*x[3]+X_BH2TSS*x[4]+X_BA2TSS*x[5]+X_P2TSS*x[6]+X_P2TSS*x[20];  // XTSS
 
  
  y[14] = u[14];                                  /* Flow */

  if (tempmodel < 0.5)                            /* Temp */ 
     y[15] = u[15];                                  
  else 
     y[15] = x[15]; 
         
  
  
  y[16] = x[16];                                   // 2 step nitrification // four step denitrification (additional) variables: SNO2, SNO, SN2O,SN2 and XBA2
  y[17] = x[17];    
  y[18] = x[18];    
  y[19] = x[19];    
  y[20] = x[20];    
  
  /* dummy states, only give outputs if ACTIVATE = 1 */
  if (activate > 0.5) {
      
  y[21] = x[21];
  y[22] = x[22];
  y[23] = x[23];
  y[24] = x[24];
  y[25] = x[25];
 
      
  }
  else if (activate < 0.5) {
  
  y[21] = 0.0;
  y[22] = 0.0;
  y[23] = 0.0;
  y[24] = 0.0;
  y[25] = 0.0;
  
  }
  
  KLa_temp = u[26]*pow(theta_kla, (u[15]-Temp_Ref)); // not x[15]?
  
  Kla_N2O = pow(D_N2O,0.5) / pow(D_O2,0.5) * KLa_temp;
  Kla_NO =  pow(D_NO,0.5) / pow(D_N2O,0.5) * Kla_N2O;
  Kla_N2 =  pow(D_N2,0.5) / pow(D_N2O,0.5) * Kla_N2O;
  
  Flux_NO =  - Kla_NO *  ((P_NO_air *  14 / H_NO) -  x[17]) * vol;
  Flux_N2O = - Kla_N2O * ((P_N2O_air * 28 / H_N2O) - x[18]) * vol;
  Flux_N2 =  - Kla_N2 *  ((P_N2_air * 28 /  H_N2) -  x[19]) * vol;
  
// calculation of free ammonia   
KB_2_KW = exp(6344.0/(273.15 + u[15]));
S_FA = (x[9]*pow(10,pH))/(KB_2_KW + pow(10,pH));

   
// calculation of free nitrous acid //
K_A = exp(-2300.0/(273.15+ u[15]));
S_FNA = (x[16] * 1.0 / (1.0 + K_A * pow(10,pH)));

  y[26] =Flux_NO; //KB_2_KW; 
  y[27] =Flux_N2O; //K_A; 
  y[28] =Flux_N2; //S_FNA ;
  
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


double proc1, proc2,proc2x1, proc2x2, proc2x3, proc2x4, proc3,proc3x1,proc3x2, proc4, proc5, proc5x1, proc5x2, proc6, proc7, proc8, proc9, proc10;
double reac1, reac2, reac3, reac4, reac5, reac6, reac7, reac8, reac9, reac10, reac11, reac12, reac13, reac16, reac17, reac18, reac19, reac20, reac21, reac22, reac23, reac24, reac25;

double vol, SO_sat, SO_sat_temp, KLa_temp;

double a_KlaN2O,	b_A1,	b_A2,	b_H,	b_KlaN2O,	D_N2,	D_N2O,	D_NO,	D_O2,	F_BOD_COD,	f_P,	F_TSS_COD,	H_N2,	H_N2O,	H_NO,	H_O2,	i_X_B,	i_X_P,	KlaN2O_anoxic,	k_a	, K_FA ,	K_FNA, k_h,	K_I10FA,	K_I10FNA,	K_I3NO,	K_I4NO,	K_I5NO,	K_I9FA,	K_I9FNA,	K_N2O,	K_NO,	K_NO2,	K_NO3,	K_OA1,	K_OA2,	K_OH,	K_OH1,	K_OH2,	K_OH3,	K_OH4,	K_OH5,	K_S1,	K_S2,	K_S3,	K_S4,	K_S5,	K_X	,mu_A1,	mu_A2,	mu_H,	n_g2,	n_g3,	n_g4,	n_g5,	n_h,	n_Y,	pH,	P_N2O_air,	P_N2_air,	P_NO_air,	P_O2_air, Y_A1, Y_A2, Y_H;
double X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS;
double b_Ratkowsky_mu_A1,	b_Ratkowsky_mu_A2,	b_Ratkowsky_mu_H,	c_Ratkowsky_mu_A1,	c_Ratkowsky_mu_A2,	c_Ratkowsky_mu_H,	Temp_Ref,	theta_b_A1,	theta_b_A2,	theta_b_H,	theta_kla,	theta_k_a,	theta_k_h,	Tmax_Ratkowsky_mu_A1,	Tmax_Ratkowsky_mu_A2,	Tmax_Ratkowsky_mu_H,	Tmin_Ratkowsky_mu_A1,	Tmin_Ratkowsky_mu_A2,	Tmin_Ratkowsky_mu_H;
double mu_H_Temp, mu_A1_Temp, mu_A2_Temp, b_H_Temp, b_A1_Temp, b_A2_Temp, k_h_Temp, k_a_Temp;
double K_SNH_aob1,	K_SNH_aob2,	K_SNO2_aob,	K_SNO_aob,	K_SO_aob1,	K_SO_aob2,	Y_aob,	Y_nob;


double S_FA, KB_2_KW, K_A, S_FNA;
double S_NOX;

double Kla_N2O, Kla_NO, Kla_N2;
double FluxN2O_gas, FluxN2_gas, FluxNO_gas;
double K_SO_AOBden1,  K_IO_AOBden1,  K_SO_AOBden2,  K_IO_AOBden2,  n_AOB,  n_Y_AOB, K_FNA_aob, K_FA_aob;

double xtemp[26];
double tempmodel;
int i;

b_A1	=           mxGetPr(PAR)[0];
b_A2 =          	mxGetPr(PAR)[1];
b_H =           	mxGetPr(PAR)[2];
D_N2	=           mxGetPr(PAR)[3];
D_N2O	=           mxGetPr(PAR)[4];
D_NO	=           mxGetPr(PAR)[5];
D_O2    =           mxGetPr(PAR)[6];
F_BOD_COD	=       mxGetPr(PAR)[7];
f_P  	=           mxGetPr(PAR)[8];
F_TSS_COD	=       mxGetPr(PAR)[9];
H_N2	=           mxGetPr(PAR)[10];
H_N2O	=           mxGetPr(PAR)[11];
H_NO	=           mxGetPr(PAR)[12];
i_X_B	=           mxGetPr(PAR)[13];
i_X_P	 =          mxGetPr(PAR)[14];
k_a              	 =mxGetPr(PAR)[15];
K_FA              	 =mxGetPr(PAR)[16];
K_FNA             	=              mxGetPr(PAR)[17];
k_h  	 =          mxGetPr(PAR)[18];
K_I10FA	=       mxGetPr(PAR)[19];
K_I10FNA	=           mxGetPr(PAR)[20];
K_I3NO	=           mxGetPr(PAR)[21];
K_I4NO	=           mxGetPr(PAR)[22];
K_I5NO	 =           mxGetPr(PAR)[23];
K_I9FA 	=           mxGetPr(PAR)[24];
K_I9FNA	=           mxGetPr(PAR)[25];
K_N2O	=           mxGetPr(PAR)[26];
K_NO	=           mxGetPr(PAR)[27];
K_NO2	=           mxGetPr(PAR)[28];
K_NO3   	=           mxGetPr(PAR)[29];
K_OA1	=           mxGetPr(PAR)[30];
K_OA2	=           mxGetPr(PAR)[31];
K_OH	=           mxGetPr(PAR)[32];
K_OH1	=           mxGetPr(PAR)[33];
K_OH2	=           mxGetPr(PAR)[34];
K_OH3	=           mxGetPr(PAR)[35];
K_OH4	=           mxGetPr(PAR)[36];
K_OH5	=           mxGetPr(PAR)[37];
K_S1	=           mxGetPr(PAR)[38];
K_S2	=           mxGetPr(PAR)[39];
K_S3	=           mxGetPr(PAR)[40];
K_S4	=           mxGetPr(PAR)[41];
K_S5	=               mxGetPr(PAR)[42];
	
K_X	=           mxGetPr(PAR)[43];
n_g2	=           mxGetPr(PAR)[44];
n_g3	=           mxGetPr(PAR)[45];
n_g4	=           mxGetPr(PAR)[46];
n_g5	    =           mxGetPr(PAR)[47];
n_h	    =           mxGetPr(PAR)[48];
n_Y	    =           mxGetPr(PAR)[49];
pH      =           mxGetPr(PAR)[50];
P_N2O_air 	=         mxGetPr(PAR)[51];
P_N2_air  	= mxGetPr(PAR)[52];
P_NO_air         	=           mxGetPr(PAR)[53];
Y_A1	=           mxGetPr(PAR)[54];
Y_A2	=               mxGetPr(PAR)[55];
Y_H	=           mxGetPr(PAR)[56]; 
	
X_I2TSS 	=           mxGetPr(PAR)[57];
X_S2TSS 	=          mxGetPr(PAR)[58];
X_BH2TSS 	=          mxGetPr(PAR)[59];
X_BA2TSS 	 =           mxGetPr(PAR)[60];
X_P2TSS	= mxGetPr(PAR)[61];
	
	
b_Ratkowsky_mu_A1=	mxGetPr(PAR)[62];
b_Ratkowsky_mu_A2=	mxGetPr(PAR)[63];
b_Ratkowsky_mu_H=	mxGetPr(PAR)[64];
c_Ratkowsky_mu_A1=	mxGetPr(PAR)[65];
c_Ratkowsky_mu_A2=	mxGetPr(PAR)[66];
c_Ratkowsky_mu_H=	mxGetPr(PAR)[67];
	
Temp_Ref	=       mxGetPr(PAR)[68];
	
theta_b_A1	=       mxGetPr(PAR)[69];
theta_b_A2	=       mxGetPr(PAR)[70];
theta_b_H	=       mxGetPr(PAR)[71];
theta_kla	=       mxGetPr(PAR)[72];
theta_k_a	=       mxGetPr(PAR)[73];
theta_k_h	=   mxGetPr(PAR)[74];
Tmax_Ratkowsky_mu_A1	=   mxGetPr(PAR)[75];
Tmax_Ratkowsky_mu_A2	    =   mxGetPr(PAR)[76];
Tmax_Ratkowsky_mu_H	=   mxGetPr(PAR)[77];
Tmin_Ratkowsky_mu_A1	=   mxGetPr(PAR)[78];
Tmin_Ratkowsky_mu_A2	    =   mxGetPr(PAR)[79];
Tmin_Ratkowsky_mu_H	=        mxGetPr(PAR)[80];
	
K_SNH_aob1 	=         mxGetPr(PAR)[81];
K_SNH_aob2	=         mxGetPr(PAR)[82];
K_SNO2_aob	= mxGetPr(PAR)[83];
K_SNO_aob	=          mxGetPr(PAR)[84];
K_SO_aob1	=          mxGetPr(PAR)[85];
K_SO_aob2	= mxGetPr(PAR)[86];
	
	
K_SO_AOBden1 	= mxGetPr(PAR)[87];
K_IO_AOBden1 	 = mxGetPr(PAR)[88];
K_SO_AOBden2	= mxGetPr(PAR)[89];
K_IO_AOBden2 	 = mxGetPr(PAR)[90];
n_AOB       	 = mxGetPr(PAR)[91];
n_Y_AOB     	 = mxGetPr(PAR)[92];
K_FNA_aob   	  = mxGetPr(PAR)[93];
K_FA_aob   	= mxGetPr(PAR)[94];


vol =               mxGetPr(V)[0];

SO_sat =            mxGetPr(SOSAT)[0];

tempmodel =         mxGetPr(TEMPMODEL)[0];

/* temperature compensation */

mu_H_Temp  =    pow((b_Ratkowsky_mu_H *  (u[15]- Tmin_Ratkowsky_mu_H) *  (1.0 - exp(c_Ratkowsky_mu_H *  (u[15]- Tmax_Ratkowsky_mu_H)))),2);
mu_A1_Temp =    pow((b_Ratkowsky_mu_A1 * (u[15]- Tmin_Ratkowsky_mu_A1) * (1.0 - exp(c_Ratkowsky_mu_A1 * (u[15]- Tmax_Ratkowsky_mu_A1)))),2);
mu_A2_Temp =    pow((b_Ratkowsky_mu_A2 * (u[15]- Tmin_Ratkowsky_mu_A2) * (1.0 - exp(c_Ratkowsky_mu_A2 * (u[15]- Tmax_Ratkowsky_mu_A2)))),2);
b_H_Temp =      b_H  * pow(theta_b_H,u[15]- Temp_Ref);
b_A1_Temp =     b_A1 * pow(theta_b_A1,u[15]- Temp_Ref);
b_A2_Temp =     b_A2 * pow(theta_b_A2,u[15]- Temp_Ref);
k_h_Temp =      k_h * pow(theta_k_h,u[15]- Temp_Ref);
k_a_Temp =      k_a * pow(theta_k_a,u[15]- Temp_Ref);

SO_sat_temp =   0.9997743214*   (8.0/  10.5*(  56.12*  6791.5*  exp(-  66.7354 + 87.4755/((  u[15]+ 273.15)/  100.0) + 24.4526*  log((u[15]+273.15)/   100.0)))); /* van't Hoff equation */
KLa_temp =      u[26]*pow(theta_kla, (u[15]-Temp_Ref));


// calculation of free ammonia   
KB_2_KW = exp(6344.0/(273.15 + u[15]));
S_FA = (x[9]*pow(10,pH))/(KB_2_KW + pow(10,pH));

   
// calculation of free nitrous acid //
K_A = exp(-2300.0/(273.15+ u[15]));
S_FNA = (x[16] * 1.0 / (1.0 + K_A * pow(10,pH)));


// calculation of the total oxidazed nitrogen
S_NOX = x[8]+x[16]+x[17]+x[18];
   
// specific KLa //
Kla_N2O = pow(D_N2O,0.5) / pow(D_O2,0.5) * KLa_temp;
Kla_NO =  pow(D_NO,0.5) / pow(D_N2O,0.5) * Kla_N2O;
Kla_N2 =  pow(D_N2,0.5) / pow(D_N2O,0.5) * Kla_N2O;

// procesess //

for (i = 0; i < 26; i++) {
   if (x[i] < 0.0)
     xtemp[i] = 0.0;
   else
     xtemp[i] = x[i];
}

if (u[26] < 0.0)
      x[7] = fabs(u[26]);


proc1 =     mu_H_Temp *        xtemp[4]* (xtemp[1] / (K_S1 + xtemp[1])) *(xtemp[7]/ (K_OH1 + xtemp[7]));
proc2x1 =   mu_H_Temp * n_g2 * xtemp[4]* (xtemp[1] / (K_S2 + xtemp[1])) *(xtemp[8] /(K_NO3 + xtemp[8])) * (K_OH2 / (K_OH2 + xtemp[7]));
proc2x2 =   mu_H_Temp * n_g3 * xtemp[4]* (xtemp[1] / (K_S3 + xtemp[1])) *(xtemp[16]/(K_NO2 + xtemp[16])) *(K_OH3 / (K_OH3 + xtemp[7])) * (K_I3NO/(K_I3NO + xtemp[17]));
proc2x3 =   mu_H_Temp * n_g4 * xtemp[4]* (xtemp[1] / (K_S4 + xtemp[1])) *(xtemp[17]/(K_NO +  xtemp[17]  + (pow(xtemp[17],2)/K_I4NO)))  * (K_OH4 /(K_OH4 + xtemp[7]  ));
proc2x4 =   mu_H_Temp * n_g5 * xtemp[4]* (xtemp[1] / (K_S5 + xtemp[1])) *(xtemp[18]/(K_N2O + xtemp[18])) *(K_OH5 / (K_OH5 + xtemp[7])) * (K_I5NO/(K_I5NO + xtemp[17]));
proc3x1 =   mu_A1_Temp * xtemp[5]  * ( S_FA / (K_FA  + S_FA +  pow(S_FA,2)  / K_I9FA  )) * (xtemp[7] / (K_OA1 + xtemp[7])) * (K_I9FNA  / (K_I9FNA + S_FNA));
proc3x2 =   mu_A2_Temp * xtemp[20] * ( S_FNA /(K_FNA + S_FNA + pow(S_FNA,2) / K_I10FNA)) * (xtemp[7] / (K_OA2 + xtemp[7])) * (K_I10FA  / (K_I10FA + S_FA ));
proc4 =     b_H_Temp * xtemp[4];
proc5x1 =   b_A1_Temp * xtemp[5];
proc5x2 =   b_A2_Temp * xtemp[20];
proc6 =     k_a_Temp * xtemp[10]*xtemp[4];
proc7 =     k_h_Temp *(xtemp[3] /xtemp[4])/(K_X +(xtemp[3]/xtemp[4]))*((xtemp[7]/(K_OH+xtemp[7]))+n_h*(K_OH/(K_OH+xtemp[7]))*(S_NOX/(K_NO3+S_NOX)))*xtemp[4];
proc8 =     proc7*(xtemp[11]/xtemp[3]);

proc9 =     n_AOB * (mu_A1_Temp * xtemp[5] * (S_FNA / (K_FNA_aob + S_FNA)) * (S_FA / (K_FA_aob + S_FA)) * (xtemp[7]/(K_SO_AOBden1 + (1-2*pow((K_SO_AOBden1/K_IO_AOBden1),0.5)) * xtemp[7] + (pow(xtemp[7],2.0)/K_IO_AOBden1))));
proc10 =    n_AOB * (mu_A1_Temp * xtemp[5] * (xtemp[17] / (K_SNO_aob +  xtemp[17])) * (S_FA / (K_FA_aob + S_FA)) * (xtemp[7]/(K_SO_AOBden2 + (1-2*pow((K_SO_AOBden2/K_IO_AOBden2),0.5)) * xtemp[7] + (pow(xtemp[7],2.0)/K_IO_AOBden2))));

// reactions //

/* SI */    reac1 = 0.0;
/* SS */    reac2 = (- 1 / Y_H)*proc1 + ( - 1 / (Y_H *n_Y))*proc2x1 + (- 1 / (Y_H * n_Y))*proc2x2 + (- 1 / (Y_H * n_Y))*proc2x3 + (- 1 / (Y_H * n_Y))*proc2x4 + proc7;;
/* XI */    reac3 = 0.0;
/* XS */    reac4 = (1 - f_P)*proc4 + (1 - f_P)*proc5x1 + (1 - f_P)*proc5x2 + (-1)*proc7;;
/* XBH */   reac5 = (1)*proc1 + (1)*proc2x1 + (1)*proc2x2 + (1)*proc2x3 + (1)*proc2x4 + (-1)*proc4;;
/* XBA1 */  reac6 =  proc3x1 - proc5x1 + proc9 + proc10  ;
/* XU */    reac7 = (f_P)*proc4 + (f_P)*proc5x1 + (f_P)*proc5x2;
/* SO2 */   reac8 =  - (1 - Y_H) / Y_H*proc1  - (3.4285714 - Y_A1) /Y_A1 *proc3x1  - (1.1428571 - Y_A2) / Y_A2*proc3x2 +(- ((2.2857143 / (Y_A1*n_Y_AOB)) - 1.0))*proc9 + (- ((2.2857143 / (Y_A1*n_Y_AOB)) - 1.0))*proc10;
/* SNO3 */  reac9 = ( - (1 -Y_H * n_Y) / (1.14385714 * Y_H * n_Y))*proc2x1 + (1 / Y_A2)*proc3x2;
/* SNH  */  reac10 = ( - i_X_B)*proc1 + ( - i_X_B)*proc2x1 + ( - i_X_B)*proc2x2 + ( - i_X_B)*proc2x3 + ( - i_X_B)*proc2x4 + ( - i_X_B - (1 / Y_A1))*proc3x1 + ( - i_X_B)*proc3x2 + (1) * proc6 + ( - (1.0 / (Y_A1*n_Y_AOB)) - i_X_B)*proc9 + ( - (1.0 / (Y_A1*n_Y_AOB)) - i_X_B)*proc10 ;
/* SND  */  reac11 = -proc6+proc8;
/* XND */   reac12 = (i_X_B - f_P * i_X_P)*proc4 + (i_X_B - f_P * i_X_P)*proc5x1 + (i_X_B - f_P * i_X_P)*proc5x2 +(-1) *proc8;
/* SALK */  reac13 = ( - i_X_B / 14)*proc1 + (- i_X_B / 14)*proc2x1 + ( - (i_X_B / 14) + (1 - Y_H * n_Y) / (14 * ((3.428571429 - 2.8571429) * Y_H * n_Y)))*proc2x2 + ( - i_X_B / 14)*proc2x3 + ( - i_X_B / 14)*proc2x4 + (( - i_X_B) / 14 - 1 / (7.0 * Y_A1))*proc3x1 + ( - i_X_B / 14.0)*proc3x2 + (1.0 / 14.0)*proc6 + (- i_X_B / 14.0)*proc9 + ( - i_X_B / 14.0) - (1.0 / (7.0 * (Y_A1*n_Y_AOB)))*proc10;


/* SNO2 */reac16 =  ((1 - Y_H * n_Y) / (1.14385714 * Y_H * n_Y))*proc2x1 + ( - (1 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x2 + (1 / Y_A1)*proc3x1 + ( - 1 / Y_A2)*proc3x2 + (- 1.0 / (Y_A1*n_Y_AOB))*proc9 + (1.0 / (Y_A1*n_Y_AOB))*proc10 ;
/* SNO  */reac17 =  ((1 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x2 + ( - (1 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x3 + ( 2.0 / (Y_A1*n_Y_AOB))*proc9 + ( - 2.0 / (Y_A1*n_Y_AOB))*proc10;;
/* SN2O */reac18 =  ((1 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x3 + ( - (1 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x4 + ( 2.0 / (Y_A1*n_Y_AOB))*proc10;
/* SN2  */reac19 =  ((1 - Y_H * n_Y) / ((3.4285714 - 2.8571429) * Y_H * n_Y))*proc2x4;
/* XBA2 */reac20 = proc3x2 -proc5x2;

reac21 = 0.0;
reac22 = 0.0;
reac23 = 0.0;
reac24 = 0.0;
reac25 = 0.0;




// mass balances //

dx[0] = 1.0/vol*(u[14]*(u[0]-x[0]))+reac1;                                                  /* SI */                                      
dx[1] = 1.0/vol*(u[14]*(u[1]-x[1]))+reac2;                                                  /* SS */     
dx[2] = 1.0/vol*(u[14]*(u[2]-x[2]))+reac3;                                                  /* XI */     
dx[3] = 1.0/vol*(u[14]*(u[3]-x[3]))+reac4;                                                  /* XS */    
dx[4] = 1.0/vol*(u[14]*(u[4]-x[4]))+reac5;                                                  /* XBH */    
dx[5] = 1.0/vol*(u[14]*(u[5]-x[5]))+reac6;                                                  /* XBA1 */    
dx[6] = 1.0/vol*(u[14]*(u[6]-x[6]))+reac7;                                                  /* XU */    
if (u[21] < 0.0)
      dx[7] = 0.0;
else
      dx[7] = 1.0/vol*(u[14]*(u[7]-x[7]))+reac8+KLa_temp*(SO_sat_temp - x[7]);

dx[8] = 1.0/vol*(u[14]*(u[8]-x[8]))+reac9;                                                  /* SNO3 */    
dx[9] = 1.0/vol*(u[14]*(u[9]-x[9]))+reac10;                                                 /* SNH */    
dx[10] = 1.0/vol*(u[14]*(u[10]-x[10]))+reac11;                                              /* SND */    
dx[11] = 1.0/vol*(u[14]*(u[11]-x[11]))+reac12;                                              /* XND */    
dx[12] = 1.0/vol*(u[14]*(u[12]-x[12]))+reac13;                                              /* SALK */    

dx[13] = 0.0;                                                                               /* TSS */

dx[14] = 0.0;                                                                               /* Flow */

if (tempmodel < 0.5)                                                                        /* Temp */    
   dx[15] = 0.0;                                  
else 
   dx[15] = 1.0/vol*(u[14]*(u[15]-x[15]));  
  
dx[16] = 1.0/vol*(u[14]*(u[16]-x[16]))+reac16;                                              /*SNO2 */

dx[17] = 1.0/vol*(u[14]*(u[17]-x[17]))+reac17 + Kla_NO * ((P_NO_air * 14 / H_NO) -  x[17]); /*SNO  */
dx[18] = 1.0/vol*(u[14]*(u[18]-x[18]))+reac18 + Kla_N2O *((P_N2O_air * 28 / H_N2O) - x[18]);/*SN2O */
dx[19] = 1.0/vol*(u[14]*(u[19]-x[19]))+reac19 + Kla_N2 * ((P_N2_air * 28 / H_N2) - x[19]);  /*SN2   */

dx[20] = 1.0/vol*(u[14]*(u[20]-x[20]))+reac20;                                              /*XBA2 */

dx[21] = 1.0/vol*(u[14]*(u[21]-x[21]))+reac21;                                              /*D1 */
dx[22] = 1.0/vol*(u[14]*(u[22]-x[22]))+reac22;                                              /*D2 */
dx[23] = 1.0/vol*(u[14]*(u[23]-x[23]))+reac23;                                              /*D3 */
dx[24] = 1.0/vol*(u[14]*(u[24]-x[24]))+reac24;                                              /*D4 */
dx[25] = 1.0/vol*(u[14]*(u[25]-x[25]))+reac25;                                              /*D5 */

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


