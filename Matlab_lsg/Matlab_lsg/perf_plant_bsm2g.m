% Effluent pollutant concentration discharge limits
totalCODemax = 100;
totalNemax = 18;
SNHemax = 4;
TSSemax = 30;
BOD5emax = 10;

% Pollutant weighting factors, effluent pollutants
BSS=2; 
BCOD=1;
BNKj=30; 
BNO=10;
BBOD5=2;

% Pumping energy factors
PF_Qintr = 0.004; % kWh/m3, pumping energy factor, internal AS recirculation
PF_Qr  = 0.008;   % kWh/m3, pumping energy factor, AS sludge recycle 
PF_Qw  = 0.05;    % kWh/m3, pumping energy factor, AS wastage flow
PF_Qpu = 0.075;   % kWh/m3, pumping energy factor, pumped underflow from primary clarifier
PF_Qtu = 0.060;   % kWh/m3, pumping energy factor, pumped underflow from thickener
PF_Qdo = 0.004;   % kWh/m3, pumping energy factor, pumped underflow from dewatering unit

%Influent_bsm2g
%PrimarySettler_bsm2g
%SludgeDisposal_bsm2g
Effluent_bsm2g

% Influent and Effluent quality index
% Note: DUMMY variables should be added here if they have COD, BOD or N content
% In this case, dummy state are included becuase they have C and N
%inf_TSS = TSS_COD * (inf_X_S + inf_X_I + inf_X_P + inf_X_BA + inf_X_BH);
Qinvec = inf_Q.*timevector;
inf_COD = inf_S_I + inf_S_S + inf_X_S + inf_X_I + inf_X_P + inf_X_BA1 + inf_X_BA2 + inf_X_BH;
inf_BOD5 = BOD_COD_inf * (inf_S_S + inf_X_S + (1 - f_P) * (inf_X_BA1 + inf_X_BA2 + inf_X_BH));
inf_TKN = inf_S_NH + inf_S_ND + inf_X_ND + i_XB * (inf_X_BA1 + inf_X_BA2 + inf_X_BH) + i_XP * (inf_X_I + inf_X_P);
inf_S_NOx = inf_S_NO3;
inf_TN = inf_TKN + inf_S_NOx;
inf_TSS = (inf_X_BH + inf_X_BA1 + inf_X_BA2 + inf_X_I + inf_X_S + inf_X_P) * TSS_COD;
inf_TN_load_aver = sum(inf_TN .*Qinvec)/(totalt*1000); %kgN/d

%eff_TSS = TSS_COD * (eff_X_S + eff_X_I + eff_X_P + eff_X_BA + eff_X_BH);
Qevec = eff_Q.*timevector;
eff_COD = eff_S_I + eff_S_S + eff_X_S + eff_X_I + eff_X_P + eff_X_BA1 + eff_X_BA2 + eff_X_BH;
eff_BOD5 = BOD_COD_eff * (eff_S_S + eff_X_S + (1 - f_P) * (eff_X_BA1 + eff_X_BA2 + eff_X_BH)); % WDK20101019
eff_TKN = eff_S_NH + eff_S_ND + eff_X_ND + i_XB * (eff_X_BA1 + eff_X_BA2 + eff_X_BH) + i_XP * (eff_X_I + eff_X_P);
eff_S_NOx = eff_S_NO3;
eff_TN = eff_TKN + eff_S_NOx;
eff_TKN_av = sum(eff_TKN.*Qevec)/Qetot;
eff_TN_av = sum(eff_TN.*Qevec)/Qetot;
eff_COD_av = sum(eff_COD.*Qevec)/Qetot;
eff_BOD5_av = sum(eff_BOD5.*Qevec)/Qetot;

EQIvec = (b_TSS * eff_TSS + b_COD * eff_COD + b_BOD5 * eff_BOD5 + b_TKN * eff_TKN + b_NO * eff_S_NOx).*Qevec;
IQIvec = (b_TSS * inf_TSS + b_COD * inf_COD + b_BOD5 * inf_BOD5 + b_TKN * inf_TKN + b_NO * inf_S_NOx).*Qinvec;

IQI=sum(IQIvec)/(totalt*1000); %kg poll.units/day
EQI=sum(EQIvec)/(totalt*1000); %kg poll.units/day

% exceedances length and frequency
% lag was not defined > changed to totalt by WDK20100906
vec = eff_BOD5;
th = max_BOD5;
BSM2_exceed
ex_d_BOD5 = sexc;
ex_p_BOD5 = sexc / totalt * 100; % exceedance length in percentage of total period
ex_f_BOD5 = freq; % exceedance frequency

vec = eff_COD;
th = max_COD;
BSM2_exceed
ex_d_COD = sexc;
ex_p_COD = sexc / totalt * 100; % exceedance length in percentage of total period
ex_f_COD = freq; % exceedance frequency

vec = eff_TSS;
th = max_TSS;
BSM2_exceed
ex_d_TSS = sexc;
ex_p_TSS = sexc / totalt * 100; % exceedance length in percentage of total period
ex_f_TSS = freq; % exceedance frequency

vec = eff_S_NH;
th = max_NH4;
BSM2_exceed
ex_d_NH4 = sexc;
ex_p_NH4 = sexc / totalt * 100; % exceedance length in percentage of total period
ex_f_NH4 = freq; % exceedance frequency

vec = eff_TN;
th = max_TN;
BSM2_exceed
ex_d_TN = sexc;
ex_p_TN = sexc / totalt * 100; % exceedance length in percentage of total period
ex_f_TN = freq; % exceedance frequency

% percentiles (95%)
prctl_eff_NH4 = prctile(eff_S_NH, 95);
prctl_eff_TN = prctile(eff_TN, 95);
prctl_eff_TSS = prctile(eff_TSS, 95);

% effluent violations (to save in XLS)
eff_viol_vars = [prctl_eff_NH4; prctl_eff_TN; prctl_eff_TSS; 0 ; ex_d_TN; ex_p_TN; ex_f_TN; 0; ex_d_COD; ex_p_COD; ex_f_COD; 0;...
    ex_d_NH4; ex_p_NH4; ex_f_NH4; 0; ex_d_TSS; ex_p_TSS; ex_f_TSS; 0; ex_d_BOD5; ex_p_BOD5; ex_f_BOD5];
