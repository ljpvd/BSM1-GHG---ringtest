%% In this mfile, a linear regression is performed on MonteCarlo
%% simulations for calculating standardized regression coefficients for
%% parameters which provides a regionalized sensitivity measure (somewhat
%% global)
% Source: Saltelli et al 2009
% Sim = Monte Carlo output matrix
% Simsc = scaled Sim
% X = input parameter matrix (i.e. LHS sampling)
% thetasc = scaled X

% Gurkan Sin, PhD @ DTU, 21 August 2007
% Modified at 19th Feb 2008

clc
clear all
close all

%% load the LHS simulations
%%Do you want to save the figures?
load MCsims20_LHS500_ExpertNoCorr_17Dec13
% clear unwanted information
clearvars -except t y1 y2 y3 y4 y5 y6 y7 y8 X
% X=X(1:20,:);
% a=X(:,1);
% b=X(:,4);
% c=X(:,5);
% d=X(:,6);
% e=X(:,7);
% f=X(:,8);
% g=X(:,9);
% h=X(:,10);
% i=X(:,11);
% 
% X=[a b c d e f g h i];

flagfig = 1 ; % 0 : No 1 : Yes
par = {'nH_Y', 'nA_Y'...
    'nPAO_Y', 'nPHA_Y', 'nH_g2', 'nH_g3', 'nH_g4', 'nH_g5', 'n1_AOB', 'n2_PAO', 'n2_AOB', 'n3_PAO', 'n4_PAO', 'n5_PAO', 'KA_2', 'KA_3', 'KA_4'...
    'KA_5', 'KF_2', 'KF_3', 'KF_4', 'KF_5', 'KH2_O2', 'KH3_O2', 'KH4_O2', 'KH5_O2', 'KH_NO2', 'KP_NO2', 'KH_NO', 'KP_NO', 'KH_N2O', 'KP_N2O'...
    'KH_I3NO', 'KP_I3NO', 'KH_I4NO', 'KP_I4NO', 'KH_I5NO', 'KP_I5NO', 'KH_O_Fe', 'KH_O_Ly', 'Mu_A1', 'Mu_A2', 'b_A1', 'b_A2', 'Y_A1', 'Y_A2', 'K_FA'...
    'K_FNA', 'K_FNA_AOBden', 'K_SNH_AOBden', 'K_SNO_AOBden', 'K_SO_AOBden', 'K_IO_AOBden', 'KOA1', 'KOA2', 'KI9FA', 'KI9FNA', 'KI10FA', 'KI10FNA' ...
        'K_NOx_Fe', 'K_NOx_Hy', 'K_NOx_Ly'};

%% this method requires scalar outputs, y:
% hence one needs to specify a meaningful property of time-series data: Let us focus on time=0.3 hr
% time=t;
% ii = find(time == 8) ;
% y1s = y1(ii,:);
% y2s = y2(ii,:);
% y3s = y3(ii,:);


% % alternatively, let us focus on the mean values:
y1s = mean(y1);
y2s = mean(y2);
y3s = mean(y3);
y4s = mean(y4);
y5s = mean(y5);
y6s = mean(y6);
y7s = mean(y7);
y8s = mean(y8);


% compile an model output matrix
Sim = [y1s' y2s' y3s' y4s' y5s' y6s' y7s' y8s'];
[n m]= size(X);

%% start of the script for building linear models on the MC outputs

%1. Make the data dimensionless (auto scaling)
msim = mean(Sim);
ssim = std(Sim);
mx = mean(X);
sx = std(X);

% scale model outputs
mus = ones(n,1)*msim ;
sds = ones(n,1)*ssim ;
Simsc = (Sim - mus) ./ sds ; 

% scale parameters
mu = ones(n,1)*mx ;
sd = ones(n,1)*sx ;
thetasc = (X - mu) ./ sd  ;
[jj k] = size(Sim);

var = {'S_{NO3}','S_{NH}','S_{NO2}','S_{NO}','S_{N2O}','S_{N2}','X_{BA1}','X_{BA2}'} ;

%2. perform linear regression FOR EACH MODEL OUTPUT
src = []; % standardized regression coefficients
srtx = []; % sorted parameters
ii = [];
rank = []; % ranking of parameters (according higher order of significance)
cov_ = []; % covariance of the coefficients of the linear model
cor_ = []; % correlation btw the coefficients of the linear model
lb=ones(m,1)* -1 ;% SRC can take a value between -1 and 1
ub=ones(m,1) * 1 ;%
for j = 1:k

    disp(['SRC number:',num2str(j)])
    [b,resnorm,residual] = lsqlin(thetasc,Simsc(:,j),[],[],[],[],lb,ub);
    y = thetasc * b ;
    R = corr2(Simsc(:,j),y);
    Rsq(j) = R^2 ;
    src(:,j) = b ;
    [srtx(:,j) ii(:,j)] = sort(abs(b),'descend');
    rank(ii(:,j),j) = 1:m;
    % below provides an estimate of covariance matrix also the
    % correlation matrix
    s2 = (sum(resnorm)/(n-m));
    cova = (sum(resnorm)/(n-m))*(thetasc'*thetasc)^-1; % covariance matrix
    for kk=1:m;
        for ll = 1:m
            core(kk,ll) = cova(kk,ll) / (sqrt(cova(kk,kk))*sqrt(cova(ll,ll)));
        end
    end
    cov_(:,:,j) =cova;
    cor_(:,:,j) =core;
    %% the sum of all the src should be 1 (incase parameters are independent), lets check
    unity(j) = b' * b;

     figure(j)
    plot(Simsc(:,j),y,'.')
    %plot(1:n,Simsc(:,j),'ro',1:n,y,'b-')
    %legend('original model outputs','linear model outputs')
    imax = find(max(y) == y);% Find the index of the min and max
    text(mean(y),mean(y),['R^2 = ',num2str(Rsq(j))])
    ylabel('Monte Carlo')
    xlabel('linear model')
    disp('')

end
yl={'Sensitivity measure (Absolute SRC)'};
xl={'Parameter significance ranking'};
ip = m;
figure
subplot(3,1,1)
plot(srtx(:,1),'.')
text(1:ip,srtx(1:ip,1),par(ii(1:ip,1)),'FontSize',8)
title(var(1))
legend('boxoff')
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')

subplot(3,1,2)
plot(srtx(:,2),'.')
%text(1:m,srtx(:,2),num2str(ii(:,2)),'FontSize',8)
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
text(1:ip,srtx(1:ip,2),par(ii(1:ip,2)),'FontSize',8)
title(var(2))
legend('boxoff')

subplot(3,1,3)
plot(srtx(:,3),'.')
text(1:ip,srtx(1:ip,3),par(ii(1:ip,3)),'FontSize',8)
title(var(3))
legend('boxoff')
ylabel(yl)
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')

% subplot(4,1,4)
% plot(srtx(:,4),'.')
% text(1:ip,srtx(1:ip,4),par(ii(1:ip,4)),'FontSize',8)
% title(var(4))
% legend('boxoff')
% xlabel(xl)
% set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')

if flagfig == 1
    flnm1 = strcat('SRCoefficient','.mat') ;
    save(flnm1,'src','srtx','Rsq','unity','rank','cov_','cor_')

    cd('Figures')
    for z=1:gcf
        f = strcat('Linearfit_mean_',num2str(z));
        saveas(z,f,'tiff');
    end

end

cd ..

yl={'SNO3','SNH4','SNO2','SNO','SN2O','SN2','XBA1','XBA2'};
disp('')
disp(yl)
disp(src)