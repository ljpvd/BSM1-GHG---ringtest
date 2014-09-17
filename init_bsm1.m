%
% Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
%
%%
% This initialization file was created//modified by Xavier Flores-Alsina
% BSM2 including GHG emissions 
% modelEAU, Department of Civil Engineering, Laval University, Quebec,(Canada)
% IEA, Div of Industrial Electrical Engineering and Automation, Lund University, Lund (Sweden)
% December 2012
%
%%
% clear
% clc

load ./Influent_data/sensornoise_bsm1;
load ./Influent_data/constinfluent_bsm1;
load ./Influent_data/dyninfluent_bsm1;

load ./Influent_data/dryinfluent_bsm1;
load ./Influent_data/storminfluent_bsm1;
load ./Influent_data/raininfluent_bsm1;

asm1init_bsm1;
settler1dinit_bsm1_UJ;
hyddelayinit_bsm1;
reginit_bsm1;
sensorinit_bsm1;


% General parameter for all subsystems
% TEMPMODEL: 0 - influent wastewater temperature is just passed through process reactors 
%            1 - mass balance for the wastewater temperature is used in
%            process reactors
% Note: thickener and dewatering are ideal models, i.e. no impact since T_out =
% T_in, in flow splitters T_out = T_in and in flow combiners mass balance based heat balance is always used.

TEMPMODEL = [ 1 ]; 

% to activate calculation of  dummy states in settler, AS reactors and storage tank set ACTIVATE = 1
ACTIVATE = [ 0 ];

% settler
MODELTYPE = [ 2 ];
