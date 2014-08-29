%
% Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
%
% initiates parameters for all sensors and actuators in use
% this file works together with reginit.m

% This initialization file was created//modified by Xavier Flores-Alsina
% BSM2 including GHG emissions 
% modelEAU, Department of Civil Engineering, Laval University, Quebec,(Canada)
% IEA, Div of Industrial Electrical Engineering and Automation, Lund University, Lund (Sweden)
% December 2010


% actuator limitations
% aeration equipment capacity (also used by controller for anti-winup)
KLa1_max = 360; %maximum possible KLa value to reactor1
KLa2_max = 360; %maximum possible KLa value to reactor2
KLa3_max = 360; %maximum possible KLa value to reactor3
KLa4_max = 360; %maximum possible KLa value to reactor4
KLa5_max = 360; %maximum possible KLa value to reactor5

% external carbon flow addition capacity (also used by controller for anti-winup)
carb1_max = 5; %maximum possible external carbon flow rate to reactor1
carb2_max = 5; %maximum possible external carbon flow rate to reactor2
carb3_max = 5; %maximum possible external carbon flow rate to reactor3
carb4_max = 5; %maximum possible external carbon flow rate to reactor4
carb5_max = 5; %maximum possible external carbon flow rate to reactor5

% pumping equipment capacity (also used by controllers for anti-windup)
Qintr_max = 5*Qin0; %maximum pump capacity for Qintr (= 103240 m3/d)
Qw_max = 0.1*Qin0; %maximum pump capacity for Qw (= 2064.8 m3/d)
Qr_max = 2*Qin0; %maximum pump capacity for Qr (= 41296 m3/d)
Qstorage_max = 1500; %maximum pump capacity for Qstorage
QwT = T*10; % time delay for artifiial Qw actuator (first-order filter)

% actuator definitions BSM2 default strategy
% KLa3 actuator, according to BSM definition, T90=4 min, n=2
T90_KLa3 = 4; %minutes
T_KLa3 = T90_KLa3/(24*60)/3.89;
useideal_KLa3 = 0; %select ideal actuator or not (0=non-ideal, 1=ideal (for testing)) for KLa3
% KLa4 actuator, according to BSM definition, T90=4 min, n=2
T90_KLa4 = 4; %minutes
T_KLa4 = T90_KLa4/(24*60)/3.89;
useideal_KLa4 = 0; %select ideal actuator or not (0=non-ideal, 1=ideal (for testing)) for KLa4
% KLa5 actuator, according to BSM definition, T90=4 min, n=2
T90_KLa5 = 4; %minutes
T_KLa5 = T90_KLa5/(24*60)/3.89;
useideal_KLa5 = 0; %select ideal actuator or not (0=non-ideal, 1=ideal (for testing)) for KLa5

%sensor definitions BSM2 deafult strategy, SO4 sensor
% DO sensor, according to BSM definition (class A), T90=1 min, n=2
min_SO4 = 0; %lower measurement limit, mg/l
max_SO4 = 10; %upper measurement limit, mg/l
T90_SO4 = 1; %response time in minutes
T_SO4 = T90_SO4/(24*60)/3.89;
std_SO4 = 0.025; %standard deviation of noise
NOISEDATA_SO4 = SENSORNOISE(:, [1 2]); %define which column in SENSORNOISE to use, column 1 = time
noiseseed_SO4 = 1; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
noisesource_SO4 = 0; %select noise source: 0=random generator, 1=predefined noisefile
usenoise_SO4 = 1; %select noise or not (0=no noise, 1=use noise) for DO sensor
useideal_SO4 = 0; %select ideal sensor or not (0=non-ideal, 1=ideal (for testing)) for DO sensor, overrides usenoise_SO5

% SNO sensor, according to BSM1 definition (class B0), T90=10 min, n=8
% min_SNO2 = 0; %lower measurement limit, 0 mg N/l
% max_SNO2 = 20; %upper measurement limit, 20 mg N/l
% T90_SNO2 = 10; %response time in minutes
% T_SNO2 = T90_SNO2/(24*60)/11.7724;
% std_SNO2 = 0.025; %standard deviation of noise
% NOISEDATA_SNO2 = SENSORNOISE(:, [1 3]); %define which column in SENSORNOISE to use, column 1 = time
% noiseseed_SNO2 = 2; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
% noisesource_SNO2 = 1; %select noise source: 0=random generator, 1=predefined noisefile
% usenoise_SNO2 = 1; %select noise or not (0=no noise, 1=use noise) for SNO sensor
% useideal_SNO2 = 0; %select ideal sensor or not (0=non-ideal, 1=ideal (for testing)) for SNO sensor, overrides usenoise_SNO2

% for a KLa actuator use simulink model KLa_actuator, and change XXX in file and simulink
% T90_XXX = 4; %minutes
% T_XXX = T90_XXX/(24*60)/3.89;
% useideal_XXX = 0; %select ideal actuator or not (0=non-ideal, 1=ideal)

% for A class sensors use simulink model sensor A, and change XXX in file and simulink
% min_XXX = ?; %lower measurement limit
% max_XXX = ?; %upper measurement limit
% T90_XXX = 1; %respose time minutes
% T_XXX = T90_XXX/(24*60)/3.89;
% std_XXX = 0.025; %standard deviation of noise
% NOISEDATA_XXX = SENSORNOISE(:, [1 2]); %define which column in SENSORNOISE to use
% noiseseed_XXX = 2; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
% noisesource_XXX = 1; %select noise source: 0=random generator, 1=predefined noisefile
% usenoise_XXX = 1; %select noise or not (0=no noise, 1=use noise) for sensor
% useideal_XXX = 0; %select ideal sensor or not (0=non-ideal, 1=ideal) for sensor

% for B0 class sensors use simulink model sensor B0, and change XXX in file and simulink
% min_XXX = ?; %lower measurement limit
% max_XXX = ?; %upper measurement limit
% T90_XXX = 10; %response time minutes
% T_XXX = T90_XXX/(24*60)/11.7724;
% std_XXX = 0.025; %standard deviation of noise
% NOISEDATA_XXX = SENSORNOISE(:, [1 3]); %define which column in SENSORNOISE to use
% noiseseed_XXX = 2; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
% noisesource_XXX = 1; %select noise source: 0=random generator, 1=predefined noisefile
% usenoise_XXX = 1; %select noise or not (0=no noise, 1=use noise) for sensor
% useideal_XXX = 0; %select ideal sensor or not (0=non-ideal, 1=ideal) for sensor

% for C0 class sensors use simulink model sensor C0, and change XXX in file and simulink
% min_XXX = ?; %lower measurement limit
% max_XXX = ?; %upper measurement limit
% T90_XXX = 20; %response time minutes
% T_XXX = T90_XXX/(24*60)/11.7724;
% std_XXX = 0.025; %standard deviation of noise
% NOISEDATA_XXX = SENSORNOISE(:, [1 4]); %define which column in SENSORNOISE to use
% noiseseed_XXX = 3; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
% noisesource_XXX = 1; %select noise source: 0=random generator, 1=predefined noisefile
% usenoise_XXX = 1; %select noise or not (0=no noise, 1=use noise) for sensor
% useideal_XXX = 0; %select ideal sensor or not (0=non-ideal, 1=ideal) for sensor

% for B1 class sensors use simulink model sensor B1, and change XXX in file and simulink
% min_XXX = ?; %lower measurement limit
% max_XXX = ?; %upper measurement limit
% T90_XXX = 10; %response time minutes
% T_XXX = T90_XXX/(24*60)/11.7724;
% T0_XXX = 5/(24*60); %sample time (minutes = 5)
% std_XXX = 0.025; %standard deviation of noise
% NOISEDATA_XXX = SENSORNOISE(:, [1 3]); %define which column in SENSORNOISE to use
% noiseseed_XXX = 2; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
% noisesource_XXX = 1; %select noise source: 0=random generator, 1=predefined noisefile
% usenoise_XXX = 1; %select noise or not (0=no noise, 1=use noise) for sensor
% useideal_XXX = 0; %select ideal sensor or not (0=non-ideal, 1=ideal) for sensor

% for C1 class sensors use simulink model sensor C1, and change XXX in file and simulink
% min_XXX = ?; %lower measurement limit
% max_XXX = ?; %upper measurement limit
% T90_XXX = 20; %response time minutes
% T_XXX = T90_XXX/(24*60)/11.7724;
% T0_XXX = 5/(24*60); %sample time (minutes = 5)
% std_XXX = 0.025; %standard deviation of noise
% NOISEDATA_XXX = SENSORNOISE(:, [1 3]); %define which column in SENSORNOISE to use
% noiseseed_XXX = 2; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
% noisesource_XXX = 1; %select noise source: 0=random generator, 1=predefined noisefile
% usenoise_XXX = 1; %select noise or not (0=no noise, 1=use noise) for sensor
% useideal_XXX = 0; %select ideal sensor or not (0=non-ideal, 1=ideal) for sensor

% for D class sensors use simulink model sensor D, and change XXX in file and simulink
% min_XXX = ?; %lower measurement limit
% max_XXX = ?; %upper measurement limit
% T0_XXX = 30/(24*60); %sample time (minutes = 30)
% std_XXX = 0.025; %standard deviation of noise
% NOISEDATA_XXX = SENSORNOISE(:, [1 3]); %define which column in SENSORNOISE to use
% noiseseed_XXX = 2; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
% noisesource_XXX = 1; %select noise source: 0=random generator, 1=predefined noisefile
% usenoise_XXX = 1; %select noise or not (0=no noise, 1=use noise) for sensor
% useideal_XXX = 0; %select ideal sensor or not (0=non-ideal, 1=ideal) for
% sensor

