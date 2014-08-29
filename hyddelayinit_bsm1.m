% This initialization file was created//modified by Xavier Flores-Alsina
% BSM2 including GHG emissions 
% modelEAU, Department of Civil Engineering, Laval University, Quebec,(Canada)
% IEA, Div of Industrial Electrical Engineering and Automation, Lund University, Lund (Sweden)
% December 2010



% hyddelay prior to AS reactors
S_I_ASin =  2.9056e+06;
S_S_ASin =  1.2920e+06;
X_I_ASin =  1.5864e+08;
X_S_ASin =  7.1895e+06;
X_BH_ASin = 2.3118e+08;
X_BA1_ASin = 1.7266e+07;
X_P_ASin =  9.9818e+07;
S_O_ASin =  1.1391e+05;
S_NO3_ASin = 7.6188e+05;
S_NH_ASin = 7.4430e+05;
S_ND_ASin = 1.6232e+05;
X_ND_ASin = 4.1923e+05;
S_ALK_ASin = 5.3815e+05;
TSS_ASin = 3.8557e+08;

Q_ASin = 1.03531e+05;
T_ASin = 14.8581;

S_NO2_ASin = 1.1391e+05;
S_NO_ASin = 1.1391e+05;
S_N2O_ASin = 1.1391e+05;
S_N2_ASin = 1.1391e+05;
X_BA2_ASin = 1.7266e+07;

S_D1_ASin = 0;
S_D2_ASin = 0;
S_D3_ASin = 0;
X_D4_ASin = 0;
X_D5_ASin = 0;


XINITDELAY = [ S_I_ASin  S_S_ASin  X_I_ASin  X_S_ASin  X_BH_ASin  X_BA1_ASin  X_P_ASin  S_O_ASin  S_NO3_ASin  S_NH_ASin  S_ND_ASin  X_ND_ASin  S_ALK_ASin TSS_ASin Q_ASin T_ASin    S_NO2_ASin S_NO_ASin S_N2O_ASin S_N2_ASin X_BA2_ASin    S_D1_ASin S_D2_ASin S_D3_ASin X_D4_ASin X_D5_ASin ];

% time constant for hyddelay (days)
T = 0.0001;