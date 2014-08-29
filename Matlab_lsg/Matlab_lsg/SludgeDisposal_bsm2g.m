%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sludge disposal concentrations

Qsvec = s_Q.*timevector;
SIsvec = s_S_I.*Qsvec;
SSsvec = s_S_S.*Qsvec;     
XIsvec = s_X_I.*Qsvec;
XSsvec = s_X_S.*Qsvec;  
XBHsvec = s_X_BH.*Qsvec;  
XBAsvec = s_X_BA.*Qsvec;
XPsvec = s_X_P.*Qsvec;
SOsvec = s_S_O.*Qsvec;
SNOsvec = s_S_NO.*Qsvec;
SNHsvec = s_S_NH.*Qsvec;
SNDsvec = s_S_ND.*Qsvec;
XNDsvec = s_X_ND.*Qsvec;
SALKsvec = s_S_ALK.*Qsvec;
TSSsvec = s_TSS.*Qsvec;
%Tempsvec = s_Temp.*Qsvec;
% if (ACTIVATE > 0.5)
    DUMMY1svec =  s_S_NO3.*Qsvec;
    DUMMY2svec =  s_S_NO2.*Qsvec;
    DUMMY3svec =  s_S_N2O.*Qsvec;
    DUMMY4svec = s_S_N2.*Qsvec;
    DUMMY5svec = s_X_BA2.*Qsvec;
% end
Qstot = sum(Qsvec);
Qsav = Qstot/totalt;

SIsload = sum(SIsvec);
SSsload = sum(SSsvec);
XIsload = sum(XIsvec);
XSsload = sum(XSsvec);
XBHsload = sum(XBHsvec);
XBAsload = sum(XBAsvec);
XPsload = sum(XPsvec);
SOsload = sum(SOsvec);
SNOsload = sum(SNOsvec);
SNHsload = sum(SNHsvec);
SNDsload = sum(SNDsvec);
XNDsload = sum(XNDsvec);
SALKsload = sum(SALKsvec);
TSSsload = sum(TSSsvec);
% Tempsload = sum(Tempsvec);
% if (ACTIVATE > 0.5)
    DUMMY1sload = sum(DUMMY1svec);
    DUMMY2sload = sum(DUMMY2svec);
    DUMMY3sload = sum(DUMMY3svec);
    DUMMY4sload = sum(DUMMY4svec);
    DUMMY5sload = sum(DUMMY5svec);
% end

SIsav = SIsload/Qstot;
SSsav = SSsload/Qstot;
XIsav = XIsload/Qstot;
XSsav = XSsload/Qstot;
XBHsav = XBHsload/Qstot;
XBAsav = XBAsload/Qstot;
XPsav = XPsload/Qstot;
SOsav = SOsload/Qstot;
SNOsav = SNOsload/Qstot;
SNHsav = SNHsload/Qstot;
SNDsav = SNDsload/Qstot;
XNDsav = XNDsload/Qstot;
SALKsav = SALKsload/Qstot;
TSSsav = TSSsload/Qstot;
% Tempsav = Tempsload/Qstot;
% if (ACTIVATE > 0.5)
    DUMMY1sav = DUMMY1sload/Qstot;
    DUMMY2sav = DUMMY2sload/Qstot;
    DUMMY3sav = DUMMY3sload/Qstot;
    DUMMY4sav = DUMMY4sload/Qstot;
    DUMMY5sav = DUMMY5sload/Qstot;
% end

TKNsav = SNHsav+SNDsav+XNDsav+i_XB*(XBHsav+XBAsav)+i_XP*(XIsav + XPsav + DUMMY5sav);%%% XBA2 as dummy state 5 is inclued in TKN
TNsav = (SNOsav + DUMMY1sav + DUMMY2sav + DUMMY3sav) + TKNsav;%%%%%%% the dissolved forms of N are included in NOx
TCODsav = SIsav+SSsav+XIsav+XSsav+XBHsav+XBAsav+XPsav + DUMMY5sav;%% XBA2 as dummy state 5 is inclued in TCOD
BOD5sav = 0.25*(SSsav+XSsav+(1-f_P)*(XBHsav + XBAsav + DUMMY5sav));%% XBA2 as dummy state 5 is inclued in BOD5

%Sludgeav=[Qsav SIsav SSsav XIsav XSsav XBHsav XBAsav XPsav SOsav SNOsav SNHsav SNDsav XNDsav SALKsav TSSsav Tempsav TKNsav TNsav TCODsav BOD5sav]';
% if (ACTIVATE > 0.5)
    Sludgeav=[Qsav SIsav SSsav XIsav XSsav XBHsav XBAsav XPsav SOsav SNOsav SNHsav SNDsav XNDsav SALKsav TSSsav TKNsav TNsav TCODsav BOD5sav DUMMY1sav DUMMY2sav DUMMY3sav DUMMY4sav DUMMY5sav ]';
% end

totalNKjsvec2=(SNHsvec+SNDsvec+XNDsvec+i_XB*(XBHsvec + XBAsvec + DUMMY5svec )+i_XP*(XPsvec+XIsvec))./Qsvec;%%% XBA2 as dummy state 5 is inclued in TKN
totalNsvec2=((SNOsvec + DUMMY1svec + DUMMY2svec + DUMMY3svec)+ SNHsvec+SNDsvec+XNDsvec+i_XB*(XBHsvec+XBAsvec+DUMMY5svec )+i_XP*(XPsvec+XIsvec))./Qsvec;%%%%%%% the dissolved forms of N are included in NOx
totalCODsvec2=(SIsvec+SSsvec+XIsvec+XSsvec+XBHsvec+XBAsvec+XPsvec + DUMMY5svec)./Qsvec;%% XBA2 as dummy state 5 is inclued in TCOD
SNHsvec2=SNHsvec./Qsvec;
TSSsvec2=TSSsvec./Qsvec;
BOD5svec2=(0.25*(SSsvec+XSsvec+(1-f_P)*(XBHsvec+ XBAsvec + DUMMY5svec)))./Qsvec;%% XBA2 as dummy state 5 is inclued in BOD5

totalNKjsload=SNHsload+SNDsload+XNDsload+i_XB*(XBHsload + XBAsload + DUMMY5sload)+i_XP*(XPsload+XIsload);%%% XBA2 as dummy state 5 is inclued in TKN
totalNsload=(SNOsload + DUMMY1sload + DUMMY2sload + DUMMY3sload) + totalNKjsload;%%%%%%% the dissolved forms of N are included in NOx
totalCODsload=(SIsload+SSsload+XIsload+XSsload+XBHsload+XBAsload+XPsload+DUMMY5sload);%% XBA2 as dummy state 5 is inclued in TCOD
BOD5sload=(0.25*(SSsload+XSsload+(1-f_P)*(XBHsload+XBAsload)));%% XBA2 as dummy state 5 is inclued in BOD5

%Sludgeload=[SIsload SSsload XIsload XSsload XBHsload XBAsload XPsload SOsload SNOsload SNHsload SNDsload XNDsload SALKsload TSSsload Tempsload totalNKjsload totalNsload totalCODsload BOD5sload]'./(1000*totalt);
% if (ACTIVATE > 0.5)
    Sludgeload=[SIsload SSsload XIsload XSsload XBHsload XBAsload XPsload SOsload SNOsload SNHsload SNDsload XNDsload SALKsload TSSsload totalNKjsload totalNsload totalCODsload BOD5sload DUMMY1sload DUMMY2sload DUMMY3sload DUMMY4sload DUMMY5sload ]'./(1000*totalt);
% end