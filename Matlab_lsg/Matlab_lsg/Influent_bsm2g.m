%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Influent concentrations

Qinvec = inf_Q.*timevector; %m3
SIinvec = inf_S_I.*Qinvec; %g
SSinvec = inf_S_S.*Qinvec;     
XIinvec = inf_X_I.*Qinvec;
XSinvec = inf_X_S.*Qinvec;  
XBHinvec = inf_X_BH.*Qinvec;  
XBAinvec = inf_X_BA.*Qinvec;
XPinvec = inf_X_P.*Qinvec;
SOinvec = inf_S_O.*Qinvec;
SNOinvec = inf_S_NO.*Qinvec;
SNHinvec = inf_S_NH.*Qinvec;
SNDinvec = inf_S_ND.*Qinvec;
XNDinvec = inf_X_ND.*Qinvec;
SALKinvec = inf_S_ALK.*Qinvec;
TSSinvec = inf_TSS.*Qinvec;
%Tempinvec = inf_Temp.*Qinvec;
% if (ACTIVATE > 0.5)
    DUMMY1invec = inf_S_NO3.*Qinvec;
    DUMMY2invec = inf_S_NO2.*Qinvec;
    DUMMY3invec = inf_S_N2O.*Qinvec;
    DUMMY4invec = inf_S_N2.*Qinvec;
    DUMMY5invec = inf_X_BA2.*Qinvec;
% end
Qintot = sum(Qinvec); %m3
Qinav = Qintot/totalt; %m3/d

SIinload = sum(SIinvec); %g
SSinload = sum(SSinvec);
XIinload = sum(XIinvec);
XSinload = sum(XSinvec);
XBHinload = sum(XBHinvec);
XBAinload = sum(XBAinvec);
XPinload = sum(XPinvec);
SOinload = sum(SOinvec);
SNOinload = sum(SNOinvec);
SNHinload = sum(SNHinvec);
SNDinload = sum(SNDinvec);
XNDinload = sum(XNDinvec);
SALKinload = sum(SALKinvec);
TSSinload = sum(TSSinvec);
%Tempinload = sum(Tempinvec);
% if (ACTIVATE > 0.5)
    DUMMY1inload = sum(DUMMY1invec);
    DUMMY2inload = sum(DUMMY2invec);
    DUMMY3inload = sum(DUMMY3invec);
    DUMMY4inload = sum(DUMMY4invec);
    DUMMY5inload = sum(DUMMY5invec);
% end

SIinav= SIinload/Qintot; %g/m3
SSinav= SSinload/Qintot;
XIinav= XIinload/Qintot;
XSinav= XSinload/Qintot;
XBHinav= XBHinload/Qintot;
XBAinav= XBAinload/Qintot;
XPinav= XPinload/Qintot;
SOinav= SOinload/Qintot;
SNOinav= SNOinload/Qintot;
SNHinav= SNHinload/Qintot;
SNDinav= SNDinload/Qintot;
XNDinav= XNDinload/Qintot;
SALKinav= SALKinload/Qintot;
TSSinav= TSSinload/Qintot;
%Tempinav= Tempinload/Qintot;
% if (ACTIVATE > 0.5)
    DUMMY1inav= DUMMY1inload/Qintot;
    DUMMY2inav= DUMMY2inload/Qintot;
    DUMMY3inav= DUMMY3inload/Qintot;
    DUMMY4inav= DUMMY4inload/Qintot;
    DUMMY5inav= DUMMY5inload/Qintot;
% end

TKNinav =        SNHinav+SNDinav+XNDinav+i_XB*(XBHinav+XBAinav +DUMMY5inav )+i_XP*(XIinav+XPinav);%%% XBA2 as dummy state 5 is inclued in TKN
TNinav =        (SNOinav + DUMMY1inav + DUMMY2inav + DUMMY3inav) + TKNinav; %%%%%%% the dissolved forms of N are included in NOx
TCODinav =       SIinav+SSinav+XIinav+XSinav+XBHinav+XBAinav+ DUMMY5inav +XPinav;%%% XBA2 as dummy state 5 is inclued in TCOD
BOD5inav =       0.65*(SSinav+XSinav+(1-f_P)*(XBHinav+ XBAinav + DUMMY5inav));%%%XBA2 as dummy state 5 is inclued in BOD5

totalNKjinvec2= (SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+XBAinvec + DUMMY5invec)+i_XP*(XPinvec+XIinvec ))./Qinvec;%%% XBA2 as dummy state 5 is inclued in TKN
totalNinvec2=   ((SNOinvec + DUMMY1invec + DUMMY2invec + DUMMY3invec) + SNHinvec+SNDinvec+XNDinvec+i_XB*(XBHinvec+ XBAinvec + DUMMY5invec)+i_XP*(XPinvec+XIinvec))./Qinvec;%%%%%%% the dissolved forms of N are included in NOx
totalCODinvec2= (SIinvec+SSinvec+XIinvec+XSinvec+XBHinvec+XBAinvec+XPinvec + DUMMY5invec)./Qinvec;%%XBA2 as dummy state 5 is inclued in TCOD
SNHinvec2=       SNHinvec./Qinvec;
TSSinvec2=       TSSinvec./Qinvec;
BOD5invec2=     (0.65*(SSinvec+XSinvec+(1-f_P)*(XBHinvec+XBAinvec + DUMMY5invec)))./Qinvec;%%%XBA2 as dummy state 5 is inclued in BOD5

totalNKjinload=  SNHinload+SNDinload+XNDinload+i_XB*(XBHinload+XBAinload + DUMMY5inload )+i_XP*(XPinload+XIinload); %%%XBA2 as dummy state 5 is inclued in TKN
totalNinload=   (SNOinload + DUMMY1inload + DUMMY2inload + DUMMY3inload) + totalNKjinload;%%%%%%% the dissolved forms of N are included in NOx
totalCODinload= (SIinload+SSinload+XIinload+XSinload+XBHinload+XBAinload+XPinload);%XBA2 as dummy state 5 is inclued in TCOD
BOD5inload=     (0.65*(SSinload+XSinload+(1-f_P)*(XBHinload+XBAinload + DUMMY5inload  )));%%%XBA2 as dummy state 5 is inclued in BOD5