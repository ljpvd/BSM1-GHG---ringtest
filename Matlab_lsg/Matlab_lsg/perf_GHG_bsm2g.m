%%% wastewater treatment and sludge processing data 

% N2O emissions 
N2Oemission1ave = sum(44/28*(N2Oemission_ASU1)./1000.*timevector)/totalt;
N2Oemission2ave = sum(44/28*(N2Oemission_ASU2)./1000.*timevector)/totalt;
N2Oemission3ave = sum(44/28*(N2Oemission_ASU3)./1000.*timevector)/totalt;
N2Oemission4ave = sum(44/28*(N2Oemission_ASU4)./1000.*timevector)/totalt;
N2Oemission5ave = sum(44/28*(N2Oemission_ASU5)./1000.*timevector)/totalt;
N2Oemissionvec = 44/28*(N2Oemission_ASU1 + N2Oemission_ASU2 + N2Oemission_ASU3 + N2Oemission_ASU4 + N2Oemission_ASU5)./1000;
N2Oemitted = sum(N2Oemissionvec.*timevector);
N2Oemittedperd = N2Oemitted/totalt;

% sludge
Qsvec = waste_Q.*timevector; % m3
TSSsvec = waste_TSS.*Qsvec; % g
TSSsload = sum(TSSsvec); % g
waste_sludge_kg = TSSsload / 1000; % g 
waste_sludge_kg_perday = waste_sludge_kg / totalt; % g/d

Qevec = eff_Q.*timevector;
TSSevec = eff_TSS.*Qevec;
TSSeload = sum(TSSevec);
eff_sludge_kg = TSSeload / 1000;
eff_sludge_kg_perday = eff_sludge_kg / totalt;

tot_sludge_kg = waste_sludge_kg + eff_sludge_kg;
tot_sludge_kg_perday = waste_sludge_kg_perday + eff_sludge_kg_perday;

% pumping energy [kWh/d]
sludge_rec_Q_m = sum(sludge_rec_Q.*timevector);
MLSS_rec_Q_m = sum(MLSS_rec_Q.*timevector);
waste_Q_m = sum(Qsvec);

pump_energy = sum((PE_waste_Q * waste_Q_m + PE_sludge_rec_Q * sludge_rec_Q + PE_MLSS_rec_Q * MLSS_rec_Q) .* timevector);% added sum(().* timevector) WDK20100906
pump_energy_perday = pump_energy / totalt;  %kWh/d

% aeration energy [kWh/d]
kla1newvec_new = SOTE*VOL1*ASU1_Kla.*O2sat1;
kla2newvec_new = SOTE*VOL2*ASU2_Kla.*O2sat2;
kla3newvec_new = SOTE*VOL3*ASU3_Kla.*O2sat3;
kla4newvec_new = SOTE*VOL4*ASU4_Kla.*O2sat4;
kla5newvec_new = SOTE*VOL5*ASU5_Kla.*O2sat5;

airenergyvec_new = (kla1newvec_new+kla2newvec_new+kla3newvec_new+kla4newvec_new+kla5newvec_new)/(1.8*1000);
airenergy_new = sum(airenergyvec_new.*timevector);
airenergy_newperd = airenergy_new/totalt; % for OCI

% aer_energy = sum((ASU1_OTR + ASU2_OTR + ASU3_OTR + ASU4_OTR + ASU5_OTR) .* timevector) / (SOTE*1000);  % changed * to .* WDK20100906
% aer_energy_perday = aer_energy / totalt; % standard BSM
       
% mixing energy cost
mixnumreac1 = length(find(ASU1_Kla<20));
mixnumreac2 = length(find(ASU2_Kla<20));
mixnumreac3 = length(find(ASU3_Kla<20));
mixnumreac4 = length(find(ASU4_Kla<20));
mixnumreac5 = length(find(ASU5_Kla<20));

mix_energy_ASU1 = mixnumreac1*ME_unit_ASU*VOL1;
mix_energy_ASU2 = mixnumreac2*ME_unit_ASU*VOL2;
mix_energy_ASU3 = mixnumreac3*ME_unit_ASU*VOL3;
mix_energy_ASU4 = mixnumreac4*ME_unit_ASU*VOL4;
mix_energy_ASU5 = mixnumreac5*ME_unit_ASU*VOL5;

mix_energy_perday = ME_C * ((mix_energy_ASU1 + mix_energy_ASU2 + mix_energy_ASU3 + mix_energy_ASU4...
    + mix_energy_ASU5)* sampletime)/totalt;

% operational cost index (to save in XLS) %kWh/d
TSScost = 5*tot_sludge_kg_perday;
airenergy_newcost = 1*airenergy_newperd; %updated BSM1
mixenergycost = 1*mix_energy_perday; %based on BSM2
pumpenergy_newcost = 1*pump_energy_perday; % based on BSM2
carbonmasscost = 0;
metalmasscost = 0;

OCI_new = TSScost + airenergy_newcost + mixenergycost + pumpenergy_newcost + carbonmasscost + metalmasscost;

%OCI = aer_energy_perday + pump_energy_perday + 3 * waste_sludge_kg_perday + 3 * C_source_perday + mix_energy_perday + heat_energy_net_perday - 6 * methane_prod_perday;  %WDK version
%OCI_energy = aer_energy_perday + pump_energy_perday + mix_energy_perday + heat_energy_net_perday - methane_energy_prod_perday; %Flores version  