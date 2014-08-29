%
% Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden
%            Xavier Flores-Alsina, IEA, Lund University, Sweden

% Default carbon addition to AS reactors
carb1 = 0; % external carbon flow rate to reactor 1
carb2 = 0; % external carbon flow rate to reactor 2
carb3 = 0; % external carbon flow rate to reactor 3
carb4 = 0; % external carbon flow rate to reactor 4
carb5 = 0; % external carbon flow rate to reactor 5
CARBONSOURCECONC = 400000; % external carbon source concentration = 400000 mg COD/l

% Default KLa values for AS reactors
KLa1 = 2;
KLa2 = 2;
KLa3 = 240;
KLa4 = 240;
KLa5 = 240;

% Default KLa values for N2O for AS reactors
a_KLaN2O1 = 0;
a_KLaN2O2 = 0;
a_KLaN2O3 = 0;
a_KLaN2O4 = 0;
a_KLaN2O5 = 0;

% Default closed loop control of Qw
Qw_high = 450;
Qw_low = 300;

% initiates parameters for all controllers in use
% this file works together with sensorinit_bsm2.m
%continuous PI O2-controller
KSO4 = 25;       %Amplification, 500 in BSM1 book
TiSO4 = 0.002;   %I-part time constant (d = 2.88 min)), integral time constant, 0.001 in BSM1 book
TtSO4 = 0.001;   %Antiwindup time constant (d), tracking time constant, 0.0002 in BSM1 book
SO4intstate = 0; %initial value of I-part
SO4awstate = 0;  %initial value of antiwindup I-part
SO4ref = 2;      %setpoint for controller, mg (-COD)/l
KLa4offset = 120; %reasonable offset value for control around SO4ref
                  % = controller output if the rest is turned off, (1/d)
useantiwindupSO4 = 1;  %0=no antiwindup, 1=use antiwindup for oxygen control
KLa3gain = 1.0; %gain for control signal to reactor 3
KLa5gain = 0.5; %gain for control signal to reactor 5

% to be used for a Qintr controller of BSM1 type, should then be
% sensorinit_bsm2 since it is an rudimentary 'actuator model' similar to QwT
%QintrT = T*10;

Qintrtype = 1;    % type 1, specified flow internally recycled, default Qw value in asm1init_bsm2.m