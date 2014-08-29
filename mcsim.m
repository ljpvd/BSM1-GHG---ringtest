%% perform monte-carlo simulations
clear all
close all
clc

% load sampling matrix
flnm = 'LHS500_ExpertNoCorr_17Dec13';

load(flnm)
clearvars -except Xcorr fln
X=Xcorr;
[n m] = size(X) ;

% specify time (in hr) % initialize:
%t = [0:(2.5/100):2.5];
t = [25:(1/96):50];
benchmarkinit
idx = [74:132,178:180]';
% run Monte Carlo simulations
for i=1:n
    i
    % update the uncertain parameters
    PAR1(idx)=X(i,:);
    PAR2=PAR1;
    PAR3=PAR1;
    PAR4=PAR1;
    PAR5=PAR1;
    PAR6=PAR1;
    PAR7=PAR1;
    PAR1
    % Solution of the model
    y = feval(@BSM1sim2,t,PAR1,XINIT1);
    y1(:,i)=y(:,1); %record simulations
    y2(:,i)=y(:,2);
    y3(:,i)=y(:,3);
    y4(:,i)=y(:,4);
    y5(:,i)=y(:,5);
    y6(:,i)=y(:,6);
    y7(:,i)=y(:,7);
    y8(:,i)=y(:,8);
end

hd1 = strcat(['MCsims20_',flnm]);
save(hd1) ;
% mc_repuncertainty
