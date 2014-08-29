% clear start
% 4 places with "%" to chage for WEST2016
clear all
close all hidden
clc
set(0,'DefaulttextInterpreter','none');
format long

b_TSS = 2;              % wheight for PI
b_COD = 1;              % wheight for PI
b_BOD5 = 2;             % wheight for PI
b_TKN = 30;             % wheight for PI
b_NO = 10;              % wheight for PI
E_kWh = 0.1;            % €/kWh
E_m3c = 70;             % €/m3 for C-source
E_ts = 100;             % €/ton of DS for sludge disposal (capital + variable)
E_m3ad = 300;           % €/m3 for construction of AD
BOD_COD_inf = 0.65;     % BOD/COD in the influent
BOD_COD_eff = 0.25;     % BOD/COD in the efffluent
TSS_COD = 0.75;         % TSS/particulateCOD in the influent
f_P = 0.08;              % Fraction of biomass to particulate products
i_XB = 0.08;            % Fraction nitrogen in biomass
i_XP = 0.06;            % Fraction nitrogen in particulate products
max_NH4 = 4;            % maximum NH4 in the effluent [g/m3]
max_TN = 18;            % maximum TN in the effluent [g/m3]
max_BOD5 = 10;          % maximum BOD5 in the effluent [g/m3]
max_COD = 100;          % maximum COD in the effluent [g/m3]
max_TSS = 30;           % maximum TSS in the effluent [g/m3]
% AD_T_oper = 308.15;     % operating temperature of AD [°K]    WDK20101014: input instead of parameter (just to be sure, probably it's constant)
AD_V_liquid = 3400;     % liquid volume of AD [m3]
slc = 30;               % service life of construction
sle = 15;               % service life of equipment
ir = 0.04;              % interest rate
VAT = 0.16;             % VAT
PE_waste_Q = 0.05;      % kWh/m3, pumping energy factor, AS wastage flow
PE_PST_Q_under = 0.075; % kWh/m3, pumping energy factor, pumped underflow from primary clarifier
PE_sludge_rec_Q = 0.008;% kWh/m3, pumping energy factor, AS sludge recycle
PE_MLSS_rec_Q = 0.004;  % kWh/m3, pumping energy factor, internal AS recirculation
PE_dewat_Q_over = 0.004;% kWh/m3, pumping energy factor, pumped underflow from dewatering unit
PE_thick_Q_under = 0.06;% kWh/m3, pumping energy factor, pumped underflow from thickener
SOTE = 1.8;             % standard oxygen transfer efficiency [kgO2/kWh]
% C_source_Q = 2;       % C source flow rate [m3/d]
C_source_COD = 400000;  % C source COD content [g/m3]
ME_unit_ASU = 0.005;    % unic cost for mixing energy calculation [kW/m3]
ME_unit_AD = 0.005;     %0.01 kW/m3 (Keller and Hartley, 2003)
ME_C = 24;              % mixing energy coefficient
MP_C1 = 16;             % methane production coefficient % molecular weight of CH4 [g/mol]
MP_C2 = 1.013;          % methane production coefficient
MP_energy = 13.8928;    % from CH4 to kWh, =50014/3600;
R = 0.083145;           % gas law constant
truck_Q = 0;           % truck load [m3], not considering truck in Flores Alsina
truck_T = 300;         % truck load temparature [°K], not considering truck in Flores Alsina
H2O_rho = 1000;         % Water density in kg/m3
H2O_cw = 4.186;         % Specific heat capacity for water in (W * s)/(gram * degC)
HE_C1 = 24;             % heating energy coefficient 
HE_C2 = 86400;          % heating energy coefficient 

% define parameters
tini = 7;             % initial time [d]
tend = 14;             % final time [d] % 100 for version 376 and 99.995 for version 2012

%[file_name, pth] = uigetfile('*.txt','Select file: ','F:\west\Project\BSM2G_test'); %WEST376
[file_name, pth] = uigetfile('*.txt','Select file: ','F:\DHI\data\projects\BSM1G'); %WEST2016
name = char(strcat(pth,file_name));
[head,data] = hdrload(name);
[m,n] = size(data);

% time
t = data(:,1); %initial values are not displayed

starttime = tini; 
stoptime = tend;

startindex=min(find(t >= starttime));
stopindex=max(find(t <= stoptime));

time_eval=t(startindex:stopindex);
sampletime = time_eval(2)-time_eval(1);
totalt=time_eval(end)-time_eval(1);
timevector = time_eval(2:end)-time_eval(1:(end-1));
nvars = (n-1);

data = data(startindex:(stopindex-1),2:n); % data matrix without time column
lt = size(data,1);

%head = head(3,1:(size(head,2))); % WEST376
head = head(1,2:(size(head,2))); % WEST2016
        [trash,head] = strtok(head);
        var_names = cell(1);
        for i = 1:nvars % creates vector of variables names
            [waste,head] = strtok(head);
            %waste = waste(2:length(waste)); %WEST376
            waste = waste(1:length(waste)); %WEST2016
            var_names{i} = waste;
        end
        clear head
BSM1G_variables_west 

perf_plant_bsm2g
perf_GHG_bsm2g

eff_concentration = [Qeav,SIeav,SSeav,XIeav,XSeav,XBHeav,XBA1eav,XPeav,SOeav,SNO3eav,SNHeav,SNDeav,XNDeav,SALKeav,TSSeav,SNO2eav,SNOeav,SN2Oeav,SN2eav,XBA2eav,eff_TKN_av,eff_TN_av,eff_COD_av,eff_BOD5_av]';
eff_load_vec = [eff_S_I.*eff_Q,eff_S_S.*eff_Q,eff_X_I.*eff_Q,eff_X_S.*eff_Q,eff_X_BH.*eff_Q,eff_X_BA1.*eff_Q,eff_X_P.*eff_Q,eff_S_O.*eff_Q,eff_S_NO3.*eff_Q,eff_S_NH.*eff_Q,eff_S_ND.*eff_Q,eff_X_ND.*eff_Q,eff_S_ALK.*eff_Q,eff_TSS.*eff_Q,eff_S_NO2.*eff_Q,eff_S_NO.*eff_Q,eff_S_N2O.*eff_Q,eff_S_N2.*eff_Q,eff_X_BA2.*eff_Q,eff_TKN.*eff_Q,eff_TN.*eff_Q,eff_COD.*eff_Q,eff_BOD5.*eff_Q];
eff_load_av = mean(eff_load_vec,1)'/1000;

energy_vars = [airenergy_newperd, pump_energy_perday, mix_energy_perday]';  

other_vars = [IQI; EQI; waste_sludge_kg_perday; eff_sludge_kg_perday; tot_sludge_kg_perday];

OCI_vars = [TSScost;airenergy_newcost;pumpenergy_newcost;carbonmasscost;mixenergycost;OCI_new];
% Performance_GHG=[CO2_emission_biotreatment, CO2_decay,CO2_BODox,CO2_credit,CO2_equivalent_N2O]';
% GHG_from_sludge_processing=CO2_emission_digester;
% power=[CO2_emission_parasitic_power,OCI_energy];
% Chemical=[CO2_emission_chemical];
% sludge_disposal=[CO2_emission_sludge_reuse];
% Other=[CO2_emission_digester, carbondioxide_prod_perday_mean, CO2_equivalent_ch4_consumed, CO2_emission_power, CO2_emission_aeration_power, CO2_emission_non_parasitic_power, credit_power, CO2_emission_chemical,CO2_emission_sludge_reuse]';

% N2Onet_Heter=[N2Onet_Heter_ASU1;N2Onet_Heter_ASU2;N2Onet_Heter_ASU3;N2Onet_Heter_ASU4;N2Onet_Heter_ASU5];
% N2Oprod_Heter=[N2Oprod_Heter_ASU1;N2Oprod_Heter_ASU2;N2Oprod_Heter_ASU3;N2Oprod_Heter_ASU4;N2Oprod_Heter_ASU5];
% N2Oremove_Heter=[N2Oremove_Heter_ASU1;N2Oremove_Heter_ASU2;N2Oremove_Heter_ASU3;N2Oremove_Heter_ASU4;N2Oremove_Heter_ASU5];
% N2Oprod_AOB=[N2Oprod_AOB_ASU1;N2Oprod_AOB_ASU2;N2Oprod_AOB_ASU3;N2Oprod_AOB_ASU4;N2Oprod_AOB_ASU5];
% N2Ostrip=sum([N2Ostrip_ASU1.*timevector,N2Ostrip_ASU2.*timevector,N2Ostrip_ASU3.*timevector,N2Ostrip_ASU4.*timevector,N2Ostrip_ASU5.*timevector],1)/totalt;
% N2Ostrip=N2Ostrip';
N2Oemissions_av = [N2Oemission1ave; N2Oemission2ave; N2Oemission3ave; N2Oemission4ave; N2Oemission5ave; N2Oemittedperd];

% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',Performance_GHG,'E3:E7');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',Other,1,'E9:E17');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',Performance_plant,1,'E21:E45');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',N2Onet_Heter,1,'E49:E53');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',N2Oprod_AOB,1,'E55:E59');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',N2Ostrip,1,'E61:E65');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',eff_ave_concentration,1,'E67:E72');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',N2Ocon_ave_concentration,1,'E74:E78');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',N2Oprod_Heter,1,'E80:E84');
% xlswrite('F:\DHI\data\projects\BSM1G\data\out.xls',N2Oremove_Heter,1,'E86:E90');