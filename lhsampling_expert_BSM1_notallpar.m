% This m file creates, plots and saves LHS sampling matrices
% created July 17, 2007 
% updated on 05 August 2011
% Gurkan Sin @DTU Chemical Engineering

clear
close all
clc

%% Define a priori probability distribution of parameters.
% we use the posterior estimated by MLE
% load MLEour
% clearvars -except pmin p95 pcor
% pmin = [0.0050 0.0001 0.0020];
pcor = eye(50);

% specify no of LHS samples
nsample = 500;
nvar = length(pcor);
init_bsm1
idx = [1:50];
par=PAR1(idx); % get the reference parameter values /we assume we dont have any data for PE. otherwise we dont use Monte Carlo
d = ones(1,50);
inputunc=d*0.50; % expert input uncertainty indicates degree of uncertainty [0: Low , 1: High]
xlu(1,:)= par .* (ones(1,nvar)-inputunc);
xlu(2,:)= par .* (ones(1,nvar)+inputunc);
cor = eye(nvar);  % expert assumes no correlation
% cor = pcor;  % expert assumes  correlation from MLE

%% Below is the call for the lhs with Iman Conover rank correlation control
%% note this file should be edited by the user to specify the prior dist
%% function (see inside the code)
[Xcorr X CX]= lhs_imancon_u(nsample,nvar,xlu,cor);

hd1 = strcat(['LHS500_','ExpertNoCorr_',datestr(date,'ddmmmyy')]);
save(hd1) ;
par = {'b_{A1}', 'b_{A2}', 'b_H', 'D_{N2}',	'D_{N2O}',	'D_{NO}', 'D_{O2}','F_BOD_COD', 'f_P', 'F_TSS_COD', 'H_{N2}', 'H_{N2O}', 'H_{NO}', 'i_X_B', 'i_X_P', 'k_a', 'K_{FA}', 'K_{FNA}'...
    'k_h', 'K_{I10FA}',	'K_{I10FNA}', 'K_{I3NO}', 'K_{I4NO}', 'K_{I5NO}', 'K_{I9FA}', 'K_{I9FNA}', 'K_{N2O}', 'K_{NO}',	'K_{NO2}',	'K_{NO3}', 'K_{OA1}', 'K_{OA2}', 'K_{OH}', 'K_{OH1}'...
    'K_{OH2}', 'K_{OH3}', 'K_{OH4}', 'K_{OH5}', 'K_{S1}', 'K_{S2}', 'K_{S3}', 'K_{S4}',	'K_{S5}', 'K_X', 'n_{g2}', 'n_{g3}', 'n_{g4}',	'n_{g5}', 'n_h', 'n_Y', 'pH'...
    'P N2O air', 'P N2 air', 'P NO air', 'Y_{A1}', 'Y_{A2}', 'Y_H', 'X I2TSS', 'X S2TSS', 'X BH2TSS', 'X BA2TSS', 'X P2TSS', 'b_{Ratkowsky mu A1}',	'b_{Ratkowsky mu A2}'...
    'b_{Ratkowsky mu H}', 'c_{Ratkowsky mu A1}', 'c_{Ratkowsky mu A2}',	'c_{Ratkowsky mu H}', 'Temp Ref', 'theta b_{A1}', 'theta b_{A2}', 'theta b_H', 'theta_{kla}', 'theta k_a'...
    'theta k_h', 'Tmax Ratkowsky mu A1', 'Tmax Ratkowsky mu A2', 'Tmax Ratkowsky mu H', 'Tmin Ratkowsky mu A1', 'Tmin Ratkowsky mu A2',	'Tmin Ratkowsky mu H', 'K_{SNH aob1}'...
    'K_{SNH aob2}',	'K_{SNO2 aob}',	'K_{SNO aob}', 'K_{SO aob1}', 'K_{SO aob2}', 'K_{SO AOBden1}', 'K_{IO AOBden1}', 'K_{SO AOBden2}', 'K_{IO AOBden2}', 'n_{AOB}', 'n_Y_{AOB}'...
    'K_{FNA aob}', 'K_{FA aob}'};

%% view results
figure(1)
[h,ax,bax,P] = plotmatrix(Xcorr) ;
set(ax,'FontSize',14,'FontWeight','bold')
for i=1:62
ylabel(ax(i,1),par(i))
xlabel(ax(62,i),par(i))
end
title('Sampling without Correlation Control')
% 
% drn='figures';
% fg1 = strcat(['CompleteViewLHS_50',num2str(1),'_',datestr(date,'ddmmmyy')]);
% f1 = fullfile(pwd,drn,fg1) ;
% saveas(1,f1,'tiff') ;


