%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effluent concentrations

Qevec = eff_Q.*timevector;
SIevec = eff_S_I.*Qevec;
SSevec = eff_S_S.*Qevec;     
XIevec = eff_X_I.*Qevec;
XSevec = eff_X_S.*Qevec;  
XBHevec = eff_X_BH.*Qevec;  
XBA1evec = eff_X_BA1.*Qevec;
XPevec = eff_X_P.*Qevec;
SOevec = eff_S_O.*Qevec;
SNOevec = eff_S_NO.*Qevec;
SNHevec = eff_S_NH.*Qevec;
SNDevec = eff_S_ND.*Qevec;
XNDevec = eff_X_ND.*Qevec;
SALKevec = eff_S_ALK.*Qevec;
TSSevec = eff_TSS.*Qevec;
%Tempevec = eff_Temp.*Qevec;
% if (ACTIVATE > 0.5)
    SNO3evec = eff_S_NO3.*Qevec;
    SNO2evec = eff_S_NO2.*Qevec;
    SN2Oevec = eff_S_N2O.*Qevec;
    SN2evec = eff_S_N2.*Qevec;
    XBA2evec = eff_X_BA2.*Qevec;
% end
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
SNOeload = sum(SNOevec);
SNHeload = sum(SNHevec);
SNDeload = sum(SNDevec);
XNDeload = sum(XNDevec);
SALKeload = sum(SALKevec);
TSSeload = sum(TSSevec);
%Tempeload = sum(Tempevec);
% if (ACTIVATE > 0.5)
    SNO3eload = sum(SNO3evec);
    SNO2eload = sum(SNO2evec);
    SN2Oeload = sum(SN2Oevec);
    SN2eload = sum(SN2evec);
    XBA2eload = sum(XBA2evec);
% end

SIeav = SIeload/Qetot;
SSeav = SSeload/Qetot;
XIeav = XIeload/Qetot;
XSeav = XSeload/Qetot;
XBHeav = XBHeload/Qetot;
XBA1eav = XBA1eload/Qetot;
XPeav = XPeload/Qetot;
SOeav = SOeload/Qetot;
SNOeav = SNOeload/Qetot;
SNHeav = SNHeload/Qetot;
SNDeav = SNDeload/Qetot;
XNDeav = XNDeload/Qetot;
SALKeav = SALKeload/Qetot;
TSSeav = TSSeload/Qetot;
%Tempeav = Tempeload/Qetot;
% if (ACTIVATE > 0.5)
    SNO3eav = SNO3eload/Qetot;
    SNO2eav = SNO2eload/Qetot;
    SN2Oeav = SN2Oeload/Qetot;
    SN2eav = SN2eload/Qetot;
    XBA2eav = XBA2eload/Qetot;
% end