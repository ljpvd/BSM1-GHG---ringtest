% Plant performance module for the BSM2
% Updated 2008-07-29 by UJ according to latest BSM1 specifications

% startime = time where plant evaluation period start (integer)
% stoptime = time where plant evaluation period stops (integer)

% Cut away first and last samples, i.e. t=smaller than starttime and 
% t = larger than stoptime. The last 52 weeks of simulated data are used to
% evaluate the plant performance. Set plotflag = 1 to activate plotting.
%
% Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden

close all

start=clock; 
disp(' ')
disp('***** Plant evaluation of BSM1 system initiated *****')
disp(['Start time for BSM1 evaluation (hour:min:sec) = ', num2str(round(start(4:6)))]); %Display start time of evaluation
disp(' ')

plotflag = 1;

starttime = 7;
stoptime = 14;
startindex=max(find(t <= starttime));
stopindex=min(find(t >= stoptime));

time=t(startindex:stopindex);

sampletime = time(2)-time(1);
totalt=time(end)-time(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The plant performance is calculated using the updated BSM2
% plant performance evaluation criteria. The main criteria are:
% - EQI as for BSM1 (updated since Copp, 2002)
% - Operational Cost Index (OCI) = 1 * pumpenergyperd
%                                + 1 * airenergyperd
%                                + 1 * mixenergyperd
%                                  5 * TSSproducedperd
%                                + 3 * carbonmassperd
%                              
% - Violations including time in violation, number of violations and % of
% time in violation as for BSM1 (Copp, 2002)
% - 95 percentile effluent concentrations for SNH, TN and TSS
% - Risk index is calculated in a separate script (called at the end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Effluent pollutant concentration discharge limits
totalCODemax = 100;
totalNemax = 18;
SNHemax = 4;
TSSemax = 30;
BOD5emax = 10;

% Pollutant weighting factors, effluent pollutants
BTSS=2; 
BCOD=1;
BNKj=30; 
BNO=10;
BBOD5=2;

% Pumping energy factors
PF_Qintr = 0.004; % kWh/m3, pumping energy factor, internal AS recirculation
PF_Qr  = 0.008;   % kWh/m3, pumping energy factor, AS sludge recycle 
PF_Qw  = 0.05;    % kWh/m3, pumping energy factor, AS wastage flow

%cut out the parts of the files to be used
inpart=in(startindex:(stopindex-1),:);
ASinputpart=ASinput(startindex:(stopindex-1),:);
reac1part=reac1(startindex:(stopindex-1),:);
reac2part=reac2(startindex:(stopindex-1),:);
reac3part=reac3(startindex:(stopindex-1),:);
reac4part=reac4(startindex:(stopindex-1),:);
reac5part=reac5(startindex:(stopindex-1),:);
settlerpart=settler(startindex:(stopindex-1),:);
effluentpart=effluent(startindex:(stopindex-1),:);
sludgepart=sludge(startindex:(stopindex-1),:);
recpart=rec(startindex:(stopindex-1),:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Influent concentrations

%timevector = time_eval(2:end)-time_eval(1:(end-1));
timevector=time(2:end)-time(1:(end-1));

Qinvec = inpart(:,15).*timevector;
SIinvec = inpart(:,1).*Qinvec;
SSinvec = inpart(:,2).*Qinvec;     
XIinvec = inpart(:,3).*Qinvec;
XSinvec = inpart(:,4).*Qinvec;  
XBHinvec = inpart(:,5).*Qinvec;  
XBA1invec = inpart(:,6).*Qinvec;
XPinvec = inpart(:,7).*Qinvec;
SOinvec = inpart(:,8).*Qinvec;
SNO3invec = inpart(:,9).*Qinvec;
SNHinvec = inpart(:,10).*Qinvec;
SNDinvec = inpart(:,11).*Qinvec;
XNDinvec = inpart(:,12).*Qinvec;
SALKinvec = inpart(:,13).*Qinvec;
TSSinvec = inpart(:,14).*Qinvec;
Tempinvec = inpart(:,16).*Qinvec;
SNO2invec = inpart(:,17).*Qinvec;
SNOinvec = inpart(:,18).*Qinvec;
SN2Oinvec = inpart(:,19).*Qinvec;
SN2invec = inpart(:,20).*Qinvec;
XBA2invec = inpart(:,21).*Qinvec;
Qintot = sum(Qinvec);
Qinav = Qintot/totalt;


SIinload = sum(SIinvec);
SSinload = sum(SSinvec);
XIinload = sum(XIinvec);
XSinload = sum(XSinvec);
XBHinload = sum(XBHinvec);
XBA1inload = sum(XBA1invec);
XPinload = sum(XPinvec);
SOinload = sum(SOinvec);
SNO3inload = sum(SNO3invec);
SNHinload = sum(SNHinvec);
SNDinload = sum(SNDinvec);
XNDinload = sum(XNDinvec);
SALKinload = sum(SALKinvec);
TSSinload = sum(TSSinvec);
Tempinload = sum(Tempinvec);
SNO2inload = sum(SNO2invec);
SNOinload = sum(SNOinvec);
SN2Oinload = sum(SN2Oinvec);
SN2inload = sum(SN2invec);
XBA2inload = sum(XBA2invec);

SIinav= SIinload/Qintot;
SSinav= SSinload/Qintot;
XIinav= XIinload/Qintot;
XSinav= XSinload/Qintot;
XBHinav= XBHinload/Qintot;
XBA1inav= XBA1inload/Qintot;
XPinav= XPinload/Qintot;
SOinav= SOinload/Qintot;
SNO3inav= SNO3inload/Qintot;
SNHinav= SNHinload/Qintot;
SNDinav= SNDinload/Qintot;
XNDinav= XNDinload/Qintot;
SALKinav= SALKinload/Qintot;
TSSinav= TSSinload/Qintot;
Tempinav= Tempinload/Qintot;
SNO2inav= SNO2inload/Qintot;
SNOinav= SNOinload/Qintot;
SN2Oinav= SN2Oinload/Qintot;
SN2inav= SN2inload/Qintot;
XBA2inav= XBA2inload/Qintot;

TKNinav =        SNHinav+SNDinav+XNDinav+i_X_B*(XBHinav+XBA1inav +XBA2inav )+i_X_P*(XIinav+XPinav);%%% XBA2 is inclued in TKN
TNinav =        (SNO3inav + SNO2inav + SNOinav + SN2Oinav) + TKNinav; %%%%%%% the dissolved forms of N are included in NOx
TCODinav =       SIinav+SSinav+XIinav+XSinav+XBHinav+XBA1inav+ XBA2inav +XPinav;%%% XBA2 is inclued in TCOD
BOD5inav =       0.65*(SSinav+XSinav+(1-f_P)*(XBHinav+ XBA1inav + XBA2inav));%%%XBA2is inclued in BOD5

totalNKjinvec2= (SNHinvec+SNDinvec+XNDinvec+i_X_B*(XBHinvec+XBA1invec + XBA1invec)+i_X_P*(XPinvec+XIinvec ))./Qinvec;%%% XBA2 is inclued in TKN
totalNinvec2=   ((SNO3invec + SNO2invec + SNOinvec + SN2Oinvec))./Qinvec + totalNKjinvec2;%%%%%%% the dissolved forms of N are included in NOx
totalCODinvec2= (SIinvec+SSinvec+XIinvec+XSinvec+XBHinvec+XBA1invec+XPinvec + XBA2invec)./Qinvec;%%XBA2  is inclued in TCOD
SNHinvec2=       SNHinvec./Qinvec;
TSSinvec2=       TSSinvec./Qinvec;
BOD5invec2=     (0.65*(SSinvec+XSinvec+(1-f_P)*(XBHinvec+XBA1invec + XBA2invec)))./Qinvec;%%%XBA2 is inclued in BOD5

totalNKjinload=  SNHinload+SNDinload+XNDinload+i_X_B*(XBHinload+XBA1inload + XBA2inload )+i_X_P*(XPinload+XIinload); %%%XBA2 is inclued in TKN
totalNinload=   (SNO3inload + SNO2inload + SNOinload + SN2Oinload) + totalNKjinload;%%%%%%% the dissolved forms of N are included in NOx
totalCODinload= (SIinload+SSinload+XIinload+XSinload+XBHinload+XBA1inload+ XBA2inload+ XPinload);%XBA2 is inclued in TCOD
BOD5inload=     (0.65*(SSinload+XSinload+(1-f_P)*(XBHinload+XBA1inload + XBA2inload  )));%%%XBA2 is inclued in BOD5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effluent concentrations

Qevec = effluentpart(:,15).*timevector;
SIevec = effluentpart(:,1).*Qevec;
SSevec = effluentpart(:,2).*Qevec;     
XIevec = effluentpart(:,3).*Qevec;
XSevec = effluentpart(:,4).*Qevec;  
XBHevec = effluentpart(:,5).*Qevec;  
XBA1evec = effluentpart(:,6).*Qevec;
XPevec = effluentpart(:,7).*Qevec;
SOevec = effluentpart(:,8).*Qevec;
SNO3evec = effluentpart(:,9).*Qevec;
SNHevec = effluentpart(:,10).*Qevec;
SNDevec = effluentpart(:,11).*Qevec;
XNDevec = effluentpart(:,12).*Qevec;
SALKevec = effluentpart(:,13).*Qevec;
TSSevec = effluentpart(:,14).*Qevec;
Tempevec = effluentpart(:,16).*Qevec;
SNO2evec = effluentpart(:,17).*Qevec;
SNOevec = effluentpart(:,18).*Qevec;
SN2Oevec = effluentpart(:,19).*Qevec;
SN2evec = effluentpart(:,20).*Qevec;
XBA2evec = effluentpart(:,21).*Qevec;

Qetot = sum(Qevec);
Qeav = Qetot/totalt;


SIeload = sum(SIevec);
SSeload = sum(SSevec);
XIeload = sum(XIevec);
XSeload = sum(XSevec);
XBHeload = sum(XBHevec);
XBA1eload = sum(XBA1evec);
XPeload = sum(XPevec);
SOeload = sum(SOevec);
SNO3eload = sum(SNO3evec);
SNHeload = sum(SNHevec);
SNDeload = sum(SNDevec);
XNDeload = sum(XNDevec);
SALKeload = sum(SALKevec);
TSSeload = sum(TSSevec);
Tempeload = sum(Tempevec);
SNO2eload = sum(SNO2evec);
SNOeload = sum(SNOevec);
SN2Oeload = sum(SN2Oevec);
SN2eload = sum(SN2evec);
XBA2eload = sum(XBA2evec);


SIeav = SIeload/Qetot;
SSeav = SSeload/Qetot;
XIeav = XIeload/Qetot;
XSeav = XSeload/Qetot;
XBHeav = XBHeload/Qetot;
XBA1eav = XBA1eload/Qetot;
XPeav = XPeload/Qetot;
SOeav = SOeload/Qetot;
SNO3eav = SNO3eload/Qetot;
SNHeav = SNHeload/Qetot;
SNDeav = SNDeload/Qetot;
XNDeav = XNDeload/Qetot;
SALKeav = SALKeload/Qetot;
TSSeav = TSSeload/Qetot;
Tempeav = Tempeload/Qetot;
SNO2eav = SNO2eload/Qetot;
SNOeav = SNOeload/Qetot;
SN2Oeav = SN2Oeload/Qetot;
SN2eav = SN2eload/Qetot;
XBA2eav = XBA2eload/Qetot;


TKNeav = SNHeav+SNDeav+XNDeav+i_X_B*(XBHeav+ XBA1eav + XBA2eav)+i_X_P*(XIeav+XPeav); %%XBA2 is inclued in TKN
TNeav = (SNO3eav + SNO2eav + SN2Oeav+ SNOeav) + TKNeav;%%%%%%% the dissolved forms of N are included in NOx
TCODeav = SIeav+SSeav+XIeav+XSeav+XBHeav+XBA1eav+XPeav + XBA2eav;%%% XBA2 is inclued in TCOD

% special to handle different BOD factors in normal effluent and bypassed effluent
BOD5_SSload = 0.25*(SSeload);
BOD5_XSload = 0.25*(XSeload);
BOD5_XBHload = 0.25*(1-f_P)*(XBHeload);
BOD5_XBA1load = 0.25*(1-f_P)*(XBA1eload);
BOD5_XBA2load = 0.25*(1-f_P)*(XBA2eload);%%%% New variable included in the balance XBA2
BOD5_XPload = 0.25*(1-f_P)*(XPeload);

BOD5eav = (BOD5_SSload + BOD5_XSload + BOD5_XBHload + BOD5_XBA1load + BOD5_XBA2load)/Qetot;%%%% XBA2 is included in the whole BOD5 balance

totalNKjevec2=(SNHevec+SNDevec+XNDevec+i_X_B*(XBHevec+ XBA1evec + XBA2evec )+i_X_P*(XPevec+XIevec))./Qevec;%%XBA2 as dummy state 5 is inclued in TKN
totalNevec2=((SNO3evec + SNO2evec + SNOevec + SN2Oevec)+ totalNKjevec2 )./Qevec;%%%%%%% the dissolved forms of N are included in NOx
totalCODevec2=(SIevec+SSevec+XIevec+XSevec+XBHevec+XBA1evec + XPevec + XBA2evec)./Qevec;%%% XBA2 as dummy state 5 is inclued in TCOD
SNHevec2=SNHevec./Qevec;
TSSevec2=TSSevec./Qevec;

% special to handle different BOD factors in normal effluent and bypassed effluent
BOD5_SSloadvec = 0.25*(SSevec);
BOD5_XSloadvec = 0.25*(XSevec);
BOD5_XBHloadvec = 0.25*(1-f_P)*(XBHevec);
BOD5_XBA1loadvec = 0.25*(1-f_P)*(XBA1evec);
BOD5_XBA2loadvec = 0.25*(1-f_P)*(XBA2evec);%%%% New variable included in the balance XBA2
BOD5_XPloadvec = 0.25*(1-f_P)*(XPevec);

BOD5evec2 = (BOD5_SSloadvec + BOD5_XSloadvec + BOD5_XBHloadvec + BOD5_XBA1loadvec + BOD5_XBA2loadvec + BOD5_XPloadvec)./Qevec;%%%% XBA2 is included in the whole BOD5 balance

totalNKjeload=SNHeload+SNDeload+XNDeload+i_X_B*(XBHeload + XBA1eload + XBA2eload)+i_X_P*(XPeload+XIeload);%%XBA2 as dummy state 5 is inclued in TKN
totalNeload=(SNO3eload + SNO2eload+ SNOeload+ SN2Oeload) + totalNKjeload;%%%%%%% the dissolved forms of N are included in NOx
totalCODeload=(SIeload+SSeload+XIeload+XSeload+XBHeload+XBA1eload+ XPeload + XBA2eload );%%% XBA2 as dummy state 5 is inclued in TCOD
BOD5eload = (BOD5_SSload + BOD5_XSload + BOD5_XBHload + BOD5_XBA1load + BOD5_XBA2load );%%%% XBA2 is included in the whole BOD5 balance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sludge disposal concentrations

Qsvec = sludgepart(:,15).*timevector;
SIsvec = sludgepart(:,1).*Qsvec;
SSsvec = sludgepart(:,2).*Qsvec;     
XIsvec = sludgepart(:,3).*Qsvec;
XSsvec = sludgepart(:,4).*Qsvec;  
XBHsvec = sludgepart(:,5).*Qsvec;  
XBA1svec = sludgepart(:,6).*Qsvec;
XPsvec = sludgepart(:,7).*Qsvec;
SOsvec = sludgepart(:,8).*Qsvec;
SNO3svec = sludgepart(:,9).*Qsvec;
SNHsvec = sludgepart(:,10).*Qsvec;
SNDsvec = sludgepart(:,11).*Qsvec;
XNDsvec = sludgepart(:,12).*Qsvec;
SALKsvec = sludgepart(:,13).*Qsvec;
TSSsvec = sludgepart(:,14).*Qsvec;
Tempsvec = sludgepart(:,16).*Qsvec;
SNO2svec = sludgepart(:,17).*Qsvec;
SNOsvec = sludgepart(:,18).*Qsvec;
SN2Osvec = sludgepart(:,19).*Qsvec;
SN2svec = sludgepart(:,20).*Qsvec;
XBA2svec = sludgepart(:,21).*Qsvec;

Qstot = sum(Qsvec);
Qsav = Qstot/totalt;

SIsload = sum(SIsvec);
SSsload = sum(SSsvec);
XIsload = sum(XIsvec);
XSsload = sum(XSsvec);
XBHsload = sum(XBHsvec);
XBA1sload = sum(XBA1svec);
XPsload = sum(XPsvec);
SOsload = sum(SOsvec);
SNO3sload = sum(SNO3svec);
SNHsload = sum(SNHsvec);
SNDsload = sum(SNDsvec);
XNDsload = sum(XNDsvec);
SALKsload = sum(SALKsvec);
TSSsload = sum(TSSsvec);
Tempsload = sum(Tempsvec);
SNO2sload = sum(SNO2svec);
SNOsload = sum(SNOsvec);
SN2Osload = sum(SN2Osvec);
SN2sload = sum(SN2svec);
XBA2sload = sum(XBA2svec);


SIsav = SIsload/Qstot;
SSsav = SSsload/Qstot;
XIsav = XIsload/Qstot;
XSsav = XSsload/Qstot;
XBHsav = XBHsload/Qstot;
XBA1sav = XBA1sload/Qstot;
XPsav = XPsload/Qstot;
SOsav = SOsload/Qstot;
SNO3sav = SNO3sload/Qstot;
SNHsav = SNHsload/Qstot;
SNDsav = SNDsload/Qstot;
XNDsav = XNDsload/Qstot;
SALKsav = SALKsload/Qstot;
TSSsav = TSSsload/Qstot;
Tempsav = Tempsload/Qstot;
SNO2sav = SNO2sload/Qstot;
SNOsav = SNOsload/Qstot;
SN2Osav = SN2Osload/Qstot;
SN2sav = SN2sload/Qstot;
XBA2sav = XBA2sload/Qstot;

TKNsav = SNHsav+SNDsav+XNDsav+i_X_B*(XBHsav+XBA1sav + XBA2sav)+i_X_P*(XIsav + XPsav + XBA2sav);%%% XBA2  is inclued in TKN
TNsav = (SNO3sav + SNO2sav + SNOsav + SN2Osav) + TKNsav;%%%%%%% the dissolved forms of N are included in NOx
TCODsav = SIsav+SSsav+XIsav+XSsav+XBHsav+XBA1sav+XPsav + XBA2sav;%% XBA2  is inclued in TCOD
BOD5sav = 0.25*(SSsav+XSsav+(1-f_P)*(XBHsav + XBA1sav + XBA2sav));%% XBA2  is inclued in BOD5

totalNKjsvec2=(SNHsvec+SNDsvec+XNDsvec+i_X_B*(XBHsvec + XBA1svec + XBA2svec )+i_X_P*(XPsvec+XIsvec))./Qsvec;%%% XBA2  is inclued in TKN
totalNsvec2=((SNO3svec + SNO2svec + SNOsvec + SN2Osvec))./Qsvec + totalNKjsvec2;%%%%%%% the dissolved forms of N are included in NOx
totalCODsvec2=(SIsvec+SSsvec+XIsvec+XSsvec+XBHsvec+XBA1svec+XPsvec + XBA2svec)./Qsvec;%% XBA2is inclued in TCOD
SNHsvec2=SNHsvec./Qsvec;
TSSsvec2=TSSsvec./Qsvec;
BOD5svec2=(0.25*(SSsvec+XSsvec+(1-f_P)*(XBHsvec+ XBA1svec + XBA2svec)))./Qsvec;%% XBA2 is inclued in BOD5

totalNKjsload=SNHsload+SNDsload+XNDsload+i_X_B*(XBHsload + XBA1sload + XBA2sload)+i_X_P*(XPsload+XIsload);%%% XBA2  is inclued in TKN
totalNsload=(SNO3sload + SNO2sload + SNOsload + SN2Osload) + totalNKjsload;%%%%%%% the dissolved forms of N are included in NOx
totalCODsload=(SIsload+SSsload+XIsload+XSsload+XBHsload+XBA1sload+XPsload+XBA2sload);%% XBA2  is inclued in TCOD
BOD5sload=(0.25*(SSsload+XSsload+(1-f_P)*(XBHsload+XBA1sload+XBA2sload)));%% XBA2 is inclued in BOD5

% Influent and Effluent quality index
% Note: DUMMY variables should be added here if they have COD, BOD or N content
% In this case, dummy state are included becuase they have C and N
TSSin=inpart(:,14);
CODin=inpart(:,1)+inpart(:,2)+inpart(:,3)+inpart(:,4)+inpart(:,5)+inpart(:,6)+inpart(:,7) + inpart(:,21);%% XBA2 is inclued in COD
SNKjin=inpart(:,10)+inpart(:,11)+inpart(:,12)+i_X_B*(inpart(:,5)+inpart(:,6) + inpart(:,21))+i_X_P*(inpart(:,3)+inpart(:,7));%%% XBA2  is inclued in TKN
SNOin=inpart(:,9) + inpart(:,17)+ inpart(:,18)+ inpart(:,19);%%%%%% the dissolved forms of N are included in NOx
BOD5in=0.65*(inpart(:,2)+inpart(:,4)+(1-f_P)*(inpart(:,5)+inpart(:,6) + inpart(:,21)));%% XBA2 is inclued in BOD

TSSe=effluentpart(:,14);
CODe=effluentpart(:,1)+effluentpart(:,2)+effluentpart(:,3)+effluentpart(:,4)+effluentpart(:,5)+effluentpart(:,6)+effluentpart(:,7)+ effluentpart(:,21);%% XBA2  is inclued in COD
SNKje=effluentpart(:,10)+effluentpart(:,11)+effluentpart(:,12)+i_X_B*(effluentpart(:,5)+effluentpart(:,6)+effluentpart(:,21))+i_X_P*(effluentpart(:,3)+effluentpart(:,7));%%% XBA2 is inclued in TKN
SNOe=effluentpart(:,9) + effluentpart(:,17) + effluentpart(:,18) + effluentpart(:,19);%% the dissolved forms of N are included in NOx
BOD5e=BOD5evec2; % XBA2 is inclued in BOD

EQIvecinst=(BTSS*TSSe+BCOD*CODe+BNKj*SNKje+BNO*SNOe+BBOD5*BOD5e).*effluentpart(:,15);

IQIvec=(BTSS*TSSin+BCOD*CODin+BNKj*SNKjin+BNO*SNOin+BBOD5*BOD5in).*Qinvec;
IQI=sum(IQIvec)/(totalt*1000);
EQIvec=(BTSS*TSSe +BCOD*CODe +BNKj*SNKje +BNO*SNOe +BBOD5*BOD5e).*Qevec;
EQI=sum(EQIvec)/(totalt*1000);

% Sludge production (calculated as kg TSS produced for the complete evaluation period)

TSSreactors_start = (reac1part(1,14)*VOL1+reac2part(1,14)*VOL2+reac3part(1,14)*VOL3+reac4part(1,14)*VOL4+reac5part(1,14)*VOL5)/1000;
TSSreactors_end = (reac1part(end,14)*VOL1+reac2part(end,14)*VOL2+reac3part(end,14)*VOL3+reac4part(end,14)*VOL4+reac5part(end,14)*VOL5)/1000;

TSSsettler_start=(settlerpart(1,54)*DIM(1)*DIM(2)/10+settlerpart(1,55)*DIM(1)*DIM(2)/10+settlerpart(1,56)*DIM(1)*DIM(2)/10+settlerpart(1,57)*DIM(1)*DIM(2)/10+settlerpart(1,58)*DIM(1)*DIM(2)/10+settlerpart(1,59)*DIM(1)*DIM(2)/10+settlerpart(1,60)*DIM(1)*DIM(2)/10+settlerpart(1,61)*DIM(1)*DIM(2)/10+settlerpart(1,62)*DIM(1)*DIM(2)/10+settlerpart(1,63)*DIM(1)*DIM(2)/10)/1000;
TSSsettler_end=(settlerpart(end,54)*DIM(1)*DIM(2)/10+settlerpart(end,55)*DIM(1)*DIM(2)/10+settlerpart(end,56)*DIM(1)*DIM(2)/10+settlerpart(end,57)*DIM(1)*DIM(2)/10+settlerpart(end,58)*DIM(1)*DIM(2)/10+settlerpart(end,59)*DIM(1)*DIM(2)/10+settlerpart(end,60)*DIM(1)*DIM(2)/10+settlerpart(end,61)*DIM(1)*DIM(2)/10+settlerpart(end,62)*DIM(1)*DIM(2)/10+settlerpart(end,63)*DIM(1)*DIM(2)/10)/1000;

TSSsludgeconc=sludgepart(:,14)/1000;  %kg/m3
Qsludgeflow=sludgepart(:,15);         %m3/d

TSSsludgevec=TSSsludgeconc.*Qsludgeflow.*timevector;

TSSproduced=sum(TSSsludgevec)+TSSreactors_end+TSSsettler_end-TSSreactors_start-TSSsettler_start;
TSSproducedperd = TSSproduced/totalt; %for OCI

Sludgetoeff=TSSeload/1000;
Sludgetoeffperd=TSSeload/(1000*totalt);

Totsludgeprod=TSSproduced+TSSeload/1000;
Totsludgeprodperd=TSSproduced/totalt+TSSeload/(1000*totalt);

% Aeration energy (calculated as kWh consumed for the complete evaluation period)
kla1vec = kla1in(startindex:(stopindex-1),:);
kla2vec = kla2in(startindex:(stopindex-1),:);
kla3vec = kla3in(startindex:(stopindex-1),:);
kla4vec = kla4in(startindex:(stopindex-1),:);
kla5vec = kla5in(startindex:(stopindex-1),:);

kla1newvec = SOSAT1*VOL1*kla1vec;
kla2newvec = SOSAT2*VOL2*kla2vec;
kla3newvec = SOSAT3*VOL3*kla3vec;
kla4newvec = SOSAT4*VOL4*kla4vec;
kla5newvec = SOSAT5*VOL5*kla5vec;
airenergyvec = (kla1newvec+kla2newvec+kla3newvec+kla4newvec+kla5newvec)/(1.8*1000);
airenergy = sum(airenergyvec.*timevector);
airenergyperd = airenergy/totalt; % for OCI

% Mixing energy (calculated as kWh consumed for the complete evaluation period)
mixnumreac1 = length(find(kla1vec<20));
mixnumreac2 = length(find(kla2vec<20));
mixnumreac3 = length(find(kla3vec<20));
mixnumreac4 = length(find(kla4vec<20));
mixnumreac5 = length(find(kla5vec<20));

mixenergyunitreac = 0.005; %kW/m3

mixenergyreac1 = mixnumreac1*mixenergyunitreac*VOL1;
mixenergyreac2 = mixnumreac2*mixenergyunitreac*VOL2;
mixenergyreac3 = mixnumreac3*mixenergyunitreac*VOL3;
mixenergyreac4 = mixnumreac4*mixenergyunitreac*VOL4;
mixenergyreac5 = mixnumreac5*mixenergyunitreac*VOL5;

mixenergyreac = 24*(mixenergyreac1+mixenergyreac2+mixenergyreac3+mixenergyreac4+mixenergyreac5)*sampletime;

mixenergyunitAD = 0.005; %0.01 kW/m3 (Keller and Hartley, 2003)

mixenergy = mixenergyreac;
mixenergyperd = mixenergy/totalt;

% Pumping energy (calculated as kWh consumed for the complete evaluation period)
Qintrflow = recpart(:,15);
Qrflow = settlerpart(:,15);
Qwflow = settlerpart(:,27);

% should we add pumping energy for storage tank pumping?

pumpenergyvec = PF_Qintr*Qintrflow+PF_Qr*Qrflow+PF_Qw*Qwflow;
pumpenergy = sum(pumpenergyvec.*timevector);
pumpenergyperd=pumpenergy/totalt;

% Carbon source addition
carbon1vec = carbon1in(startindex:(stopindex-1),:);
carbon2vec = carbon2in(startindex:(stopindex-1),:);
carbon3vec = carbon3in(startindex:(stopindex-1),:);
carbon4vec = carbon4in(startindex:(stopindex-1),:);
carbon5vec = carbon5in(startindex:(stopindex-1),:);
Qcarbonvec = (carbon1vec+carbon2vec+carbon3vec+carbon4vec+carbon5vec);
carbonmassvec = Qcarbonvec*CARBONSOURCECONC/1000;
Qcarbon = sum(Qcarbonvec.*timevector); %m3
carbonmass = sum(carbonmassvec.*timevector); %kg COD
carbonmassperd = carbonmass/totalt; %for OCI

% Operational Cost Index, OCI
TSScost=5*TSSproducedperd;
airenergycost=1*airenergyperd;
mixenergycost=1*mixenergyperd;
pumpenergycost=1*pumpenergyperd;
carbonmasscost=3*carbonmassperd;

OCI=TSScost+airenergycost+mixenergycost+pumpenergycost+carbonmasscost;

% Calculate 95% percentiles for effluent SNH, TN and TSS
SNHeffprctile=prctile(SNHevec2,95);
TNeffprctile=prctile(totalNevec2,95);
TSSeffprctile=prctile(TSSevec2,95);

disp(' ')
disp(['Overall plant performance during time ',num2str(time(1)),' to ',num2str(time(end)),' days'])
disp('*****************************************************')
disp(' ')
disp('Influent average concentrations based on load')
disp('---------------------------------------------')
disp(['Influent average flow rate = ',num2str(Qinav),' m3/d'])
disp(['Influent average SI conc = ',num2str(SIinav),' g COD/m3'])
disp(['Influent average SS conc = ',num2str(SSinav),' g COD/m3'])
disp(['Influent average XI conc = ',num2str(XIinav),' g COD/m3'])
disp(['Influent average XS conc = ',num2str(XSinav),' g COD/m3'])
disp(['Influent average XBH conc = ',num2str(XBHinav),' g COD/m3'])
disp(['Influent average XBA1 conc = ',num2str(XBA1inav),' g COD/m3'])
disp(['Influent average XP conc = ',num2str(XPinav),' g COD/m3'])
disp(['Influent average SO conc = ',num2str(SOinav),' g (-COD)/m3'])
disp(['Influent average SNO3 conc = ',num2str(SNO3inav),' g N/m3'])
disp(['Influent average SNH conc = ',num2str(SNHinav),' g N/m3  (limit = 4 g N/m3)'])
disp(['Influent average SND conc = ',num2str(SNDinav),' g N/m3'])
disp(['Influent average XND conc = ',num2str(XNDinav),' g N/m3'])
disp(['Influent average SALK conc = ',num2str(SALKinav),' mol HCO3/m3'])
disp(['Influent average TSS conc = ',num2str(TSSinav),' g SS/m3  (limit = 30 g SS/m3)'])
disp(['Influent average Temperature = ',num2str(Tempinav),' degC'])
disp(['Influent average SNO2 conc = ',num2str(SNO2inav),' g N/m3'])
disp(['Influent average SNO conc = ',num2str(SNOinav),' g N/m3'])
disp(['Influent average SN2O conc = ',num2str(SN2Oinav),' g N/m3'])
disp(['Influent average SN2 conc = ',num2str(SN2inav),' g N/m3'])
disp(['Influent average XBA2 conc = ',num2str(XBA2inav),' g COD/m3'])

disp(' ')
disp(['Influent average Kjeldahl N conc = ',num2str(TKNinav),' g N/m3'])
disp(['Influent average total N conc = ',num2str(TNinav),' g N/m3  (limit = 18 g N/m3)'])
disp(['Influent average total COD conc = ',num2str(TCODinav),' g COD/m3  (limit = 100 g COD/m3)'])
disp(['Influent average BOD5 conc = ',num2str(BOD5inav),' g/m3  (limit = 10 g/m3)'])
disp(' ')
disp('Influent average load')
disp('---------------------')
disp(['Influent average SI load = ',num2str(SIinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average SS load = ',num2str(SSinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XI load = ',num2str(XIinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XS load = ',num2str(XSinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XBH load = ',num2str(XBHinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XBA1 load = ',num2str(XBA1inload/(1000*totalt)),' kg COD/day'])
disp(['Influent average XP load = ',num2str(XPinload/(1000*totalt)),' kg COD/day'])
disp(['Influent average SO load = ',num2str(SOinload/(1000*totalt)),' kg (-COD)/day'])
disp(['Influent average SNO3 load = ',num2str(SNO3inload/(1000*totalt)),' kg N/day'])
disp(['Influent average SNH load = ',num2str(SNHinload/(1000*totalt)),' kg N/day'])
disp(['Influent average SND load = ',num2str(SNDinload/(1000*totalt)),' kg N/day'])
disp(['Influent average XND load = ',num2str(XNDinload/(1000*totalt)),' kg N/day'])
disp(['Influent average SALK load = ',num2str(SALKinload/(1000*totalt)),' kmol HCO3/day'])
disp(['Influent average TSS load = ',num2str(TSSinload/(1000*totalt)),' kg SS/day'])
disp(['Influent average SNO2 load = ',num2str(SNO2inload/(1000*totalt)),' kg N/day'])
disp(['Influent average SNO load = ',num2str(SNOinload/(1000*totalt)),' kg N/day'])
disp(['Influent average SNO load = ',num2str(SN2Oinload/(1000*totalt)),' kg N/day'])
disp(['Influent average SN2O load = ',num2str(SN2inload/(1000*totalt)),' kg N/day'])
disp(['Influent average SN2 load = ',num2str(XBA2inload/(1000*totalt)),' kg COD/day'])

disp(' ')
disp(['Influent average Kjeldahl N load = ',num2str(totalNKjinload/(1000*totalt)),' kg N/d'])
disp(['Influent average total N load = ',num2str(totalNinload/(1000*totalt)),' kg N/d'])
disp(['Influent average total COD load = ',num2str(totalCODinload/(1000*totalt)),' kg COD/d'])
disp(['Influent average BOD5 load = ',num2str(BOD5inload/(1000*totalt)),' kg BOD5/d'])


disp(' ')
disp('Effluent average concentrations based on load')
disp('---------------------------------------------')
disp(['Effluent average flow rate = ',num2str(Qeav),' m3/d'])
disp(['Effluent average SI conc = ',num2str(SIeav),' g COD/m3'])
disp(['Effluent average SS conc = ',num2str(SSeav),' g COD/m3'])
disp(['Effluent average XI conc = ',num2str(XIeav),' g COD/m3'])
disp(['Effluent average XS conc = ',num2str(XSeav),' g COD/m3'])
disp(['Effluent average XBH conc = ',num2str(XBHeav),' g COD/m3'])
disp(['Effluent average XBA1 conc = ',num2str(XBA1eav),' g COD/m3'])
disp(['Effluent average XP conc = ',num2str(XPeav),' g COD/m3'])
disp(['Effluent average SO conc = ',num2str(SOeav),' g (-COD)/m3'])
disp(['Effluent average SNO3 conc = ',num2str(SNO3eav),' g N/m3'])
disp(['Effluent average SNH conc = ',num2str(SNHeav),' g N/m3  (limit = 4 g N/m3)'])
disp(['Effluent average SND conc = ',num2str(SNDeav),' g N/m3'])
disp(['Effluent average XND conc = ',num2str(XNDeav),' g N/m3'])
disp(['Effluent average SALK conc = ',num2str(SALKeav),' mol HCO3/m3'])
disp(['Effluent average TSS conc = ',num2str(TSSeav),' g SS/m3  (limit = 30 g SS/m3)'])
disp(['Effluent average Temperature = ',num2str(Tempeav),' degC'])
disp(['Effluent average SNO2 conc = ',num2str(SNO2eav),' g N/m3'])
disp(['Effluent average SNO conc = ',num2str(SNOeav),' g N/m3'])
disp(['Effluent average SN2O conc = ',num2str(SN2Oeav),' g N/m3'])
disp(['Effluent average SN2 conc = ',num2str(SN2eav),' g N/m3'])
disp(['Effluent average XBA2 conc = ',num2str(XBA2eav),' g COD/m3'])

disp(' ')
disp(['Effluent average Kjeldahl N conc = ',num2str(TKNeav),' g N/m3'])
disp(['Effluent average total N conc = ',num2str(TNeav),' g N/m3  (limit = 18 g N/m3)'])
disp(['Effluent average total COD conc = ',num2str(TCODeav),' g COD/m3  (limit = 100 g COD/m3)'])
disp(['Effluent average BOD5 conc = ',num2str(BOD5eav),' g/m3  (limit = 10 g/m3)'])
disp(' ')
disp('Effluent average load')
disp('---------------------')
disp(['Effluent average SI load = ',num2str(SIeload/(1000*totalt)),' kg COD/day'])
disp(['Effluent average SS load = ',num2str(SSeload/(1000*totalt)),' kg COD/day'])
disp(['Effluent average XI load = ',num2str(XIeload/(1000*totalt)),' kg COD/day'])
disp(['Effluent average XS load = ',num2str(XSeload/(1000*totalt)),' kg COD/day'])
disp(['Effluent average XBH load = ',num2str(XBHeload/(1000*totalt)),' kg COD/day'])
disp(['Effluent average XBA1 load = ',num2str(XBA1eload/(1000*totalt)),' kg COD/day'])
disp(['Effluent average XP load = ',num2str(XPeload/(1000*totalt)),' kg COD/day'])
disp(['Effluent average SO load = ',num2str(SOeload/(1000*totalt)),' kg (-COD)/day'])
disp(['Effluent average SNO3 load = ',num2str(SNO3eload/(1000*totalt)),' kg N/day'])
disp(['Effluent average SNH load = ',num2str(SNHeload/(1000*totalt)),' kg N/day'])
disp(['Effluent average SND load = ',num2str(SNDeload/(1000*totalt)),' kg N/day'])
disp(['Effluent average XND load = ',num2str(XNDeload/(1000*totalt)),' kg N/day'])
disp(['Effluent average SALK load = ',num2str(SALKeload/(1000*totalt)),' kmol HCO3/day'])
disp(['Effluent average TSS load = ',num2str(TSSeload/(1000*totalt)),' kg SS/day'])
disp(['Effluent average SNO2 load = ',num2str(SNO2eload/(1000*totalt)),' kg N/day'])
disp(['Effluent average SNO load = ',num2str(SNOeload/(1000*totalt)),' kg N/day'])
disp(['Effluent average SNO load = ',num2str(SN2Oeload/(1000*totalt)),' kg N/day'])
disp(['Effluent average SN2O load = ',num2str(SN2eload/(1000*totalt)),' kg N/day'])
disp(['Effluent average SN2 load = ',num2str(XBA2eload/(1000*totalt)),' kg COD/day'])

disp(' ')
disp(['Effluent average Kjeldahl N load = ',num2str(totalNKjeload/(1000*totalt)),' kg N/d'])
disp(['Effluent average total N load = ',num2str(totalNeload/(1000*totalt)),' kg N/d'])
disp(['Effluent average total COD load = ',num2str(totalCODeload/(1000*totalt)),' kg COD/d'])
disp(['Effluent average BOD5 load = ',num2str(BOD5eload/(1000*totalt)),' kg BOD5/d'])

disp(' ')
disp('Sludge for disposal average concentrations based on load')
disp('--------------------------------------------------------')
disp(['Sludge for disposal average flow rate = ',num2str(Qsav),' m3/d'])
disp(['Sludge for disposal average SI conc = ',num2str(SIsav),' g COD/m3'])
disp(['Sludge for disposal average SS conc = ',num2str(SSsav),' g COD/m3'])
disp(['Sludge for disposal average XI conc = ',num2str(XIsav),' g COD/m3'])
disp(['Sludge for disposal average XS conc = ',num2str(XSsav),' g COD/m3'])
disp(['Sludge for disposal average XBH conc = ',num2str(XBHsav),' g COD/m3'])
disp(['Sludge for disposal average XBA1 conc = ',num2str(XBA1sav),' g COD/m3'])
disp(['Sludge for disposal average XP conc = ',num2str(XPsav),' g COD/m3'])
disp(['Sludge for disposal average SO conc = ',num2str(SOsav),' g (-COD)/m3'])
disp(['Sludge for disposal average SNO3 conc = ',num2str(SNO3sav),' g N/m3'])
disp(['Sludge for disposal average SNH conc = ',num2str(SNHsav),' g N/m3'])
disp(['Sludge for disposal average SND conc = ',num2str(SNDsav),' g N/m3'])
disp(['Sludge for disposal average XND conc = ',num2str(XNDsav),' g N/m3'])
disp(['Sludge for disposal average SALK conc = ',num2str(SALKsav),' mol HCO3/m3'])
disp(['Sludge for disposal average TSS conc = ',num2str(TSSsav),' g SS/m3'])
disp(['Sludge for disposal average Temperature = ',num2str(Tempsav),' degC'])
disp(['Sludge for disposal average SNO2 conc = ',num2str(SNO2sav),' g N/m3'])
disp(['Sludge for disposal average SNO conc = ',num2str(SNOsav),' g N/m3'])
disp(['Sludge for disposal average SN2O conc = ',num2str(SN2Osav),' g N/m3'])
disp(['Sludge for disposal average SN2 conc = ',num2str(SN2sav),' g N/m3'])
disp(['Sludge for disposal average XBA2 conc = ',num2str(XBA2sav),' g COD/m3'])

disp(' ')
disp(['Sludge for disposal average Kjeldahl N conc = ',num2str(TKNsav),' g N/m3'])
disp(['Sludge for disposal average total N conc = ',num2str(TNsav),' g N/m3'])
disp(['Sludge for disposal average total COD conc = ',num2str(TCODsav),' g COD/m3'])
disp(['Sludge for disposal average BOD5 conc = ',num2str(BOD5sav),' g BOD5/m3'])
disp(' ')
disp('Sludge for disposal average load')
disp('--------------------------------')
disp(['Sludge for disposal average SI load = ',num2str(SIsload/(1000*totalt)),' kg COD/day'])
disp(['Sludge for disposal average SS load = ',num2str(SSsload/(1000*totalt)),' kg COD/day'])
disp(['Sludge for disposal average XI load = ',num2str(XIsload/(1000*totalt)),' kg COD/day'])
disp(['Sludge for disposal average XS load = ',num2str(XSsload/(1000*totalt)),' kg COD/day'])
disp(['Sludge for disposal average XBH load = ',num2str(XBHsload/(1000*totalt)),' kg COD/day'])
disp(['Sludge for disposal average XBA load = ',num2str(XBA1sload/(1000*totalt)),' kg COD/day'])
disp(['Sludge for disposal average XP load = ',num2str(XPsload/(1000*totalt)),' kg COD/day'])
disp(['Sludge for disposal average SO load = ',num2str(SOsload/(1000*totalt)),' kg (-COD)/day'])
disp(['Sludge for disposal average SNO3 load = ',num2str(SNO3sload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average SNH load = ',num2str(SNHsload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average SND load = ',num2str(SNDsload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average XND load = ',num2str(XNDsload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average SALK load = ',num2str(SALKsload/(1000*totalt)),' kmol HCO3/day'])
disp(['Sludge for disposal average TSS load = ',num2str(TSSsload/(1000*totalt)),' kg SS/day'])
disp(['Sludge for disposal average SNO2 load = ',num2str(SNO2sload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average SNO load = ',num2str(SNOsload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average SN2O load = ',num2str(SN2Osload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average SN2 load = ',num2str(SN2sload/(1000*totalt)),' kg N/day'])
disp(['Sludge for disposal average XBA2 load = ',num2str(XBA2sload/(1000*totalt)),' kg N/day'])

disp(' ')
disp(['Sludge for disposal average Kjeldahl N load = ',num2str(totalNKjsload/(1000*totalt)),' kg N/d'])
disp(['Sludge for disposal average total N load = ',num2str(totalNsload/(1000*totalt)),' kg N/d'])
disp(['Sludge for disposal average total COD load = ',num2str(totalCODsload/(1000*totalt)),' kg COD/d'])
disp(['Sludge for disposal average BOD5 load = ',num2str(BOD5sload/(1000*totalt)),' kg BOD5/d'])
disp(' ')
disp('Other effluent quality variables')
disp('--------------------------------')
disp(['Influent Quality Index (IQI) = ',num2str(IQI),' kg poll.units/d'])
disp(['Effluent Quality Index (EQI) = ',num2str(EQI),' kg poll.units/d'])
disp(' ')
disp(['Sludge production for disposal = ',num2str(TSSproduced),' kg SS'])
disp(['Average sludge production for disposal per day = ',num2str(TSSproducedperd),' kg SS/d'])
disp(['Sludge production released into effluent = ',num2str(Sludgetoeff),' kg SS'])
disp(['Average sludge production released into effluent per day = ',num2str(Sludgetoeffperd),' kg SS/d'])
disp(['Total sludge production = ',num2str(Totsludgeprod),' kg SS'])
disp(['Total average sludge production per day = ',num2str(Totsludgeprodperd),' kg SS/d'])

disp(' ')
disp(['Average aeration energy = ',num2str(airenergyperd),' kWh/d'])
disp(['Average pumping energy = ',num2str(pumpenergyperd),' kWh/d'])
disp(['Average carbon source addition = ',num2str(carbonmassperd),' kg COD/d'])
disp(['Average mixing energy = ',num2str(mixenergyperd),' kWh/d'])
disp(' ')

disp('Operational Cost Index')
disp('----------------------')
disp(['Sludge production cost index = ',num2str(TSScost)])
disp(['Aeration energy cost index = ',num2str(airenergycost)])
disp(['Pumping energy cost index = ',num2str(pumpenergycost)])
disp(['Carbon source dosage cost index = ',num2str(carbonmasscost)])
disp(['Mixing energy cost index = ',num2str(mixenergycost)])
disp(['Total Operational Cost Index (OCI) = ',num2str(OCI)])
disp(' ')
disp('N2O Emissions')
disp('-------------------')
GHG_reac1part = GHG_reac1(startindex:(stopindex-1),:);
GHG_reac2part = GHG_reac2(startindex:(stopindex-1),:);
GHG_reac3part = GHG_reac3(startindex:(stopindex-1),:);
GHG_reac4part = GHG_reac4(startindex:(stopindex-1),:);
GHG_reac5part = GHG_reac5(startindex:(stopindex-1),:);

N2Oemissionvec1 = (GHG_reac1part(:,2))./1000;
N2Oemissionvec2 = (GHG_reac2part(:,2))./1000;
N2Oemissionvec3 = (GHG_reac3part(:,2))./1000;
N2Oemissionvec4 = (GHG_reac4part(:,2))./1000;
N2Oemissionvec5 = (GHG_reac5part(:,2))./1000;

N2Oemissionvec1_t = sum(N2Oemissionvec1.*timevector);
N2Oemissionvec2_t = sum(N2Oemissionvec2.*timevector);
N2Oemissionvec3_t = sum(N2Oemissionvec3.*timevector);
N2Oemissionvec4_t = sum(N2Oemissionvec4.*timevector);
N2Oemissionvec5_t = sum(N2Oemissionvec5.*timevector);

N2Oemissionvec1_perd = N2Oemissionvec1_t/totalt;
N2Oemissionvec2_perd = N2Oemissionvec2_t/totalt;
N2Oemissionvec3_perd = N2Oemissionvec3_t/totalt;
N2Oemissionvec4_perd = N2Oemissionvec4_t/totalt;
N2Oemissionvec5_perd = N2Oemissionvec5_t/totalt;

N2Oemissionvec = (GHG_reac1part(:,2)+ GHG_reac2part(:,2)+ GHG_reac3part(:,2) + GHG_reac4part(:,2)+ GHG_reac5part(:,2))./1000;
N2Oemitted = sum(N2Oemissionvec.*timevector);
N2Oemittedperd = N2Oemitted/totalt;

disp(['N2O emission during nitrification/denitrification (ANOX1) = ',num2str(N2Oemissionvec1_perd),' kg N2O/d'])
disp(['N2O emission during nitrification/denitrification (ANOX2) = ',num2str(N2Oemissionvec2_perd),' kg N2O/d'])
disp(['N2O emission during nitrification/denitrification (AER1) = ',num2str(N2Oemissionvec3_perd),' kg N2O/d'])
disp(['N2O emission during nitrification/denitrification (AER2) = ',num2str(N2Oemissionvec4_perd),' kg N2O/d'])
disp(['N2O emission during nitrification/denitrification (AER3) = ',num2str(N2Oemissionvec5_perd),' kg N2O/d'])
disp(['N2O emission during nitrification/denitrification (total) = ',num2str(N2Oemittedperd),' kg N2O/d'])
disp(' ')
disp('Effluent violations')
disp('-------------------')
disp(['95% percentile for effluent SNH (Ammonia95) = ',num2str(SNHeffprctile),' g N/m3'])
disp(['95% percentile for effluent TN (TN95) = ',num2str(TNeffprctile),' g N/m3'])
disp(['95% percentile for effluent TSS (TSS95) = ',num2str(TSSeffprctile),' g SS/m3'])
disp(' ')


%output=[Effluentav; Effload; IQI; EQI; TSSproduced; TSSproducedperd; Sludgetoeff; Sludgetoeffperd; Totsludgeprod; Totsludgeprodperd; airenergyperd; pumpenergyperd; carbonmassperd; mixenergyperd; Heatenergyperd; Methaneprodperd; TSScost; airenergycost; pumpenergycost; carbonmasscost; mixenergycost; Heatenergycost; EnergyfromMethaneperdcost; OCI; SNHeffprctile; TNeffprctile; TSSeffprctile];

%save results output;

Nviolation=find(totalNevec2>totalNemax);
CODviolation=find(totalCODevec2>totalCODemax);
SNHviolation=find(SNHevec2>SNHemax);
TSSviolation=find(TSSevec2>TSSemax);
BOD5violation=find(BOD5evec2>BOD5emax);

noofNviolation = 1;
noofCODviolation = 1;
noofSNHviolation = 1;
noofTSSviolation = 1;
noofBOD5violation = 1;

if not(isempty(Nviolation))
  disp('The maximum effluent total nitrogen level (18 g N/m3) was violated')
  disp(['during ',num2str(min(totalt,length(Nviolation)*sampletime)),' days, i.e. ',num2str(min(100,length(Nviolation)*sampletime/totalt*100)),'% of the operating time.'])
  Nviolationtime=min(totalt,length(Nviolation)*sampletime);
  Nviolationtimepercent=min(100,length(Nviolation)*sampletime/totalt*100);
  for i=2:length(Nviolation)
    if Nviolation(i-1)~=(Nviolation(i)-1)
      noofNviolation = noofNviolation+1;
    end
  end
  disp(['The limit was violated at ',num2str(noofNviolation),' different occasions.'])
  disp(' ')
  %output=[output; Nviolationtime; Nviolationtimepercent; noofNviolation];
end

if not(isempty(CODviolation))
  disp('The maximum effluent total COD level (100 g COD/m3) was violated')
  disp(['during ',num2str(min(totalt,length(CODviolation)*sampletime)),' days, i.e. ',num2str(min(100,length(CODviolation)*sampletime/totalt*100)),'% of the operating time.'])
  CODviolationtime=min(totalt,length(CODviolation)*sampletime);
  CODviolationtimepercent=min(100,length(CODviolation)*sampletime/totalt*100);
  for i=2:length(CODviolation)
    if CODviolation(i-1)~=(CODviolation(i)-1)
      noofCODviolation = noofCODviolation+1;
    end
  end
  disp(['The limit was violated at ',num2str(noofCODviolation),' different occasions.'])
  disp(' ')
  %output=[output; CODviolationtime; CODviolationtimepercent; noofCODviolation];
end

if not(isempty(SNHviolation))
  disp('The maximum effluent ammonia nitrogen level (4 g N/m3) was violated')
  disp(['during ',num2str(min(totalt,length(SNHviolation)*sampletime)),' days, i.e. ',num2str(min(100,length(SNHviolation)*sampletime/totalt*100)),'% of the operating time.'])
  SNHviolationtime=min(totalt,length(SNHviolation)*sampletime);
  SNHviolationtimepercent=min(100,length(SNHviolation)*sampletime/totalt*100);
  for i=2:length(SNHviolation)
    if SNHviolation(i-1)~=(SNHviolation(i)-1)
      noofSNHviolation = noofSNHviolation+1;
    end
  end
  disp(['The limit was violated at ',num2str(noofSNHviolation),' different occasions.'])
  disp(' ')
  %output=[output; SNHviolationtime; SNHviolationtimepercent; noofSNHviolation];
end

if not(isempty(TSSviolation))
  disp('The maximum effluent total suspended solids level (30 g SS/m3) was violated')
  disp(['during ',num2str(min(totalt,length(TSSviolation)*sampletime)),' days, i.e. ',num2str(min(100,length(TSSviolation)*sampletime/totalt*100)),'% of the operating time.'])
  TSSviolationtime=min(totalt,length(TSSviolation)*sampletime);
  TSSviolationtimepercent=min(100,length(TSSviolation)*sampletime/totalt*100);
  for i=2:length(TSSviolation)
    if TSSviolation(i-1)~=(TSSviolation(i)-1)
      noofTSSviolation = noofTSSviolation+1;
    end
  end
  disp(['The limit was violated at ',num2str(noofTSSviolation),' different occasions.'])
  disp(' ')
  %output=[output; TSSviolationtime; TSSviolationtimepercent; noofTSSviolation];
end

if not(isempty(BOD5violation))
  disp('The maximum effluent BOD5 level (10 mg/l) was violated')
  disp(['during ',num2str(min(totalt,length(BOD5violation)*sampletime)),' days, i.e. ',num2str(min(100,length(BOD5violation)*sampletime/totalt*100)),'% of the operating time.'])
  BOD5violationtime=min(totalt,length(BOD5violation)*sampletime);
  BOD5violationtimepercent=min(100,length(BOD5violation)*sampletime/totalt*100);
  for i=2:length(BOD5violation)
    if BOD5violation(i-1)~=(BOD5violation(i)-1)
      noofBOD5violation = noofBOD5violation+1;
    end
  end
  disp(['The limit was violated at ',num2str(noofBOD5violation),' different occasions.'])
  disp(' ')
  %output=[output; BOD5violationtime; BOD5violationtimepercent; noofBOD5violation];
end

%Inputs

% if plotflag==1
%     disp(' ')
%     disp('Plotting of BSM2 evaluation results has been initiated')
%     disp('******************************************************')
%     disp(' ')
%     movingaveragewindow = 96; % even number
%     timeshift = movingaveragewindow/2;
%     b = ones(1,movingaveragewindow)./movingaveragewindow;
%     
%     figure(1)
%     plot(time_eval(1:(end-1)),totalNevec2,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,totalNevec2);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     plot([time_eval(1) time_eval(end-1)],[totalNemax totalNemax],'r','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('TN in effluent (g N/m^3)','FontSize',10,'FontWeight','bold')
%     title('Effluent total nitrogen (raw and filtered) and limit value','FontSize',10,'FontWeight','bold')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(2)
%     plot(time_eval(1:(end-1)),totalCODevec2,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,totalCODevec2);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     plot([time_eval(1) time_eval(end-1)],[totalCODemax totalCODemax],'r','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Total COD in effluent (g COD/m^3)','FontSize',10,'FontWeight','bold')
%     title('Effluent total COD (raw and filtered) and limit value','FontSize',10,'FontWeight','bold')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(3)
%     plot(time_eval(1:(end-1)),SNHevec2,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,SNHevec2);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     plot([time_eval(1) time_eval(end-1)],[SNHemax SNHemax],'r','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Ammonia in effluent (g N/m^3)','FontSize',10,'FontWeight','bold')
%     title('Effluent total ammonia (raw and filtered) and limit value','FontSize',10,'FontWeight','bold')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(4)
%     plot(time_eval(1:(end-1)),TSSevec2,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,TSSevec2);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     plot([time_eval(1) time_eval(end-1)],[TSSemax TSSemax],'r','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Suspended solids in effluent (g SS/m^3)','FontSize',10,'FontWeight','bold')
%     title('Effluent suspended solids (raw and filtered) and limit value','FontSize',10,'FontWeight','bold')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(5)
%     plot(time_eval(1:(end-1)),BOD5evec2,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,BOD5evec2);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     plot([time_eval(1) time_eval(end-1)],[BOD5emax BOD5emax],'r','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('BOD5 in effluent (g/m^3)','FontSize',10,'FontWeight','bold')
%     title('Effluent BOD5 (raw and filtered) and limit value','FontSize',10,'FontWeight','bold')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(6)
%     plot(time_eval(1:(end-1)),EQIvecinst./1000,'b','LineWidth',1)
%     hold on;
%     filteredout=filter(b,1,EQIvecinst./1000);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Instantaneous EQ index (kg poll.units/d)','FontSize',10,'FontWeight','bold')
%     title('Effluent Quality Index (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(7)
%     plot(time_eval(1:(end-1)),TSSsludgeconc.*Qsludgeflow,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,TSSsludgeconc.*Qsludgeflow);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Instantaneous sludge wastage rate (kg SS/d)','FontSize',10,'FontWeight','bold')
%     title('Sludge wastage rate (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
%     
%     figure(8)
%     plot(time_eval(1:(end-1)),pumpenergyvec,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,pumpenergyvec);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Instantaneous pumping energy (kWh/d)','FontSize',10,'FontWeight','bold')
%     title('Pumping energy (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(9)
%     plot(time_eval(1:(end-1)),airenergyvec,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,airenergyvec);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Instantaneous aeration energy (kWh/d)','FontSize',10,'FontWeight','bold')
%     title('Aeration energy (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(10)
%     plot(time_eval(1:(end-1)),carbonmassvec,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,carbonmassvec);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Instantaneous carbon source addition (kg COD/d)','FontSize',10,'FontWeight','bold')
%     title('Carbon source addition (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
%    
%     figure(11)
%     plot(time_eval(1:(end-1)),Methaneflowvec,'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,Methaneflowvec);
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Instantaneous methane production (kg CH_4/d)','FontSize',10,'FontWeight','bold')
%     title('Methane production (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
%        
%     figure(12)
%     plot(time_eval(1:(end-1)),digesteroutpart(:,51),'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,digesteroutpart(:,51));
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Instantaneous total gas flow (Nm^3/d)','FontSize',10,'FontWeight','bold')
%     title('Total gas flow from AD normalized to P-atm (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
%     
%     figure(13)
%     plot(time_eval(1:(end-1)),storagepart(:,22),'b','LineWidth',1)
%     hold on
%     filteredout=filter(b,1,storagepart(:,22));
%     filteredout=filteredout(movingaveragewindow:end);
%     plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1.5)
%     xlabel('time (days)','FontSize',10,'FontWeight','bold')
%     ylabel('Liquid volume (m^3)','FontSize',10,'FontWeight','bold')
%     title('Liquid volume in storage tank (raw and filtered)')
%     hold off
%     xlim([starttime stoptime])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
%     
%     % Plot of the SNH, TN and TSS curves
%     SNHeffsort=sort(SNHevec2);
%     TNeffsort=sort(totalNevec2);
%     TSSeffsort=sort(TSSevec2);
%     n=size(SNHevec2,1);
%     xvalues=[1:n].*(100/n);
% 
%     figure(14)
%     plot(xvalues,SNHeffsort,'b','LineWidth',1)
%     hold on
%     plot([0 95],[SNHeffprctile SNHeffprctile],'k--','LineWidth',1.5)
%     hold on
%     plot([95 95],[0 SNHeffprctile],'k--','LineWidth',1.5)
%     xlabel('Ordered S_N_H effluent concentrations (%)','FontSize',10,'FontWeight','bold')
%     ylabel('S_N_H effluent concentrations (g N/m^3)','FontSize',10,'FontWeight','bold')
%     title('Ordered effluent S_N_H concentrations with 95% percentile','FontSize',10,'FontWeight','bold')
%     xlim([0 105])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(15)
%     plot(xvalues,TNeffsort,'b','LineWidth',1)
%     hold on
%     plot([0 95],[TNeffprctile TNeffprctile],'k--','LineWidth',1.5)
%     hold on
%     plot([95 95],[0 TNeffprctile],'k--','LineWidth',1.5)
%     xlabel('Ordered TN effluent concentrations (%)','FontSize',10,'FontWeight','bold')
%     ylabel('TN effluent concentrations (g N/m^3)','FontSize',10,'FontWeight','bold')
%     title('Ordered effluent TN concentrations with 95% percentile','FontSize',10,'FontWeight','bold')
%     xlim([0 105])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
% 
%     figure(16)
%     plot(xvalues,TSSeffsort,'b','LineWidth',1)
%     hold on
%     plot([0 95],[TSSeffprctile TSSeffprctile],'k--','LineWidth',1.5)
%     hold on
%     plot([95 95],[0 TSSeffprctile],'k--','LineWidth',1.5)
%     xlabel('Ordered TSS effluent concentrations (%)','FontSize',10,'FontWeight','bold')
%     ylabel('TSS effluent concentrations (g SS/m^3)','FontSize',10,'FontWeight','bold')
%     title('Ordered effluent TSS concentrations with 95% percentile','FontSize',10,'FontWeight','bold')
%     xlim([0 105])
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
%     
%     disp('Plotting of BSM2 evaluation results has been completed')
%     disp('******************************************************')
%     disp(' ')
% 
% end


stop=clock;
disp('***** Plant evaluation of BSM2 system successfully finished *****')
disp(['End time (hour:min:sec) = ', num2str(round(stop(4:6)))]); %Display simulation stop time
disp(' ')