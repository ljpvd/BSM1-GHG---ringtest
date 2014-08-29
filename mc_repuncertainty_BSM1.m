%% Representation of uncertainty of Monte Carlo outputs
%% plotting empirical CDFs
% close all
% clear all
% clc
% 
% dr0 = pwd;
% ylb = 'Cumulative probability';
% load MCsims_corr
%% one needs to use scalar values for CDF hence, we need to focus on a
%% meaningful property of time-series data: Let us focus on time=0.3 hr
ts = [0:(1/96):50];
time=ts;
starttime = 25;
stoptime = 50;
startindex=max(find(ts <= starttime));
stopindex=min(find(ts >= stoptime));
t=ts(startindex:(stopindex-1));

ii = find(time == 7) ;

y1s = y1(ii,:);

y2s = y2(ii,:);

y3s = y3(ii,:);

y4s = y4(ii,:);

y5s = y5(ii,:);

y6s = y6(ii,:);

y7s = y7(ii,:);

y8s = y8(ii,:);

mk = 'k';
lw = 2.0 ;
%LineStyleOrder
lso = {'k-','r:','b--','ko','ro','bo','k:','r-'};
kk=1;

figure(1)
subplot(2,2,1)
y = y1s;
mu(1)=mean(y);
var(1)=cov(y);
st(1)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx); % find the mean value!
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SNO3')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')
ylabel(gca,'Prob x<X')

subplot(2,2,2)
y = y2s;
mu(2,kk)=mean(y);
var(2,kk)=cov(y);
st(2,kk)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx);
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SNH4')
ylabel(gca,'Prob x<X')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')

subplot(2,2,3)
y = y3s;
mu(3,kk)=mean(y);
var(3,kk)=cov(y);
st(3,kk)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx);
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SNO2')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')
ylabel(gca,'Prob x<X')

subplot(2,2,4)
y = y4s;
mu(1)=mean(y);
var(1)=cov(y);
st(1)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx); % find the mean value!
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SNO')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')
ylabel(gca,'Prob x<X')

figure(2)
subplot(2,2,1)
y = y5s;
mu(2,kk)=mean(y);
var(2,kk)=cov(y);
st(2,kk)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx);
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SN2O')
ylabel(gca,'Prob x<X')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')

subplot(2,2,2)
y = y6s;
mu(3,kk)=mean(y);
var(3,kk)=cov(y);
st(3,kk)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx);
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'SN2')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')
ylabel(gca,'Prob x<X')

subplot(2,2,3)
y = y7s;
mu(1)=mean(y);
var(1)=cov(y);
st(1)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx); % find the mean value!
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'XBA1')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')
ylabel(gca,'Prob x<X')

subplot(2,2,4)
y = y8s;
mu(2,kk)=mean(y);
var(2,kk)=cov(y);
st(2,kk)=std(y);
[f,x] = ecdf(y);
xx =mean(y);
ff = interp1q(x,f,xx);
plot(x,f,mk,'LineWidth',lw)
text(xx,ff,['\rightarrow','\mu']) %,'=',num2str(xx,'%10.2f')
grid('on')
xlabel(gca,'XBA2')
ylabel(gca,'Prob x<X')
set(gca,'YTick',[0:0.25:1.0])
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold')


figure(3)
%plot raw data
subplot(2,2,1)
plot(t,y1)
ylabel('SNO3')
xlabel('Time (days)')
subplot(2,2,2)
plot(t,y2)
ylabel('SNH4')
xlabel('Time (days)')
xlim([7 time(end)])
% ylim([min(y2(:,1)) max(y2(:,1))])
subplot(2,2,3)
plot(t,y3)
ylabel('SNO2')
xlabel('Time (days)')
subplot(2,2,4)
plot(t,y3)
ylabel('SNO')
xlabel('Time (days)')

figure(4)
%plot raw data
subplot(2,2,1)
plot(t,y1)
ylabel('SN2O')
xlabel('Time (days)')
subplot(2,2,2)
plot(t,y2)
ylabel('SN2')
xlabel('Time (days)')
xlim([7 time(end)])
% ylim([min(y2(:,1)) max(y2(:,1))])
subplot(2,2,3)
plot(t,y3)
ylabel('XBA1')
xlabel('Time (days)')
subplot(2,2,4)
plot(t,y3)
ylabel('XBA2')
xlabel('Time (days)')

pcr=5; % probability at which to evaluate the distributions in % (i.e. 100% is 1)
figure(5)
subplot(2,2,1)
yy = [prctile(y1',100-pcr); mean(y1'); prctile(y1',pcr);std(y1')]';
%h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3),t,yy(:,2)+yy(:,4),t,yy(:,2)-yy(:,4));
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('SNO3')
xlabel('Time (days)')
xlim([5 t(end)])
legend('95^{th} percentile','mean','5^{th} percentile')
legend('boxoff')
%set(gca,'XTick',[])
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')


subplot(2,2,2)
yy = [prctile(y2',100-pcr); mean(y2'); prctile(y2',pcr);std(y2')]';
%plot(tt(:,1:3),yy)
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('SNH4')
xlabel('Time (days)')
%set(gca,'XTick',[])
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
xlim([5 t(end)])

subplot(2,2,3)
yy = [prctile(y3',100-pcr); mean(y3'); prctile(y3',pcr);std(y3')]';
% plot(tt(:,1:3),yy)
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('SNO2')
xlabel('Time (days)')
%set(gca,'XTick',0:0.5)
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
xlim([5 t(end)])

subplot(2,2,4)
yy = [prctile(y4',100-pcr); mean(y4'); prctile(y4',pcr);std(y4')]';
% plot(tt(:,1:3),yy)
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('SNO')
xlabel('Time (days)')
%set(gca,'XTick',0:0.5)
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
xlim([5 t(end)])

figure(6)
subplot(2,2,1)
yy = [prctile(y5',100-pcr); mean(y5'); prctile(y5',pcr);std(y5')]';
%h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3),t,yy(:,2)+yy(:,4),t,yy(:,2)-yy(:,4));
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('SN2O')
xlabel('Time (days)')
xlim([5 t(end)])
legend('95^{th} percentile','mean','5^{th} percentile')
legend('boxoff')
%set(gca,'XTick',[])
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')


subplot(2,2,2)
yy = [prctile(y6',100-pcr); mean(y6'); prctile(y6',pcr);std(y6')]';
%plot(tt(:,1:3),yy)
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('SN2')
xlabel('Time (days)')
%set(gca,'XTick',[])
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
xlim([5 t(end)])

subplot(2,2,3)
yy = [prctile(y7',100-pcr); mean(y7'); prctile(y7',pcr);std(y7')]';
% plot(tt(:,1:3),yy)
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('XBA1')
xlabel('Time (days)')
%set(gca,'XTick',0:0.5)
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
xlim([5 t(end)])

subplot(2,2,4)
yy = [prctile(y8',100-pcr); mean(y8'); prctile(y8',pcr);std(y8')]';
% plot(tt(:,1:3),yy)
h = plot(t,yy(:,1),t,yy(:,2),t,yy(:,3));
set(h(1),'Color','r','LineWidth',lw,'LineStyle','--')
set(h(2),'Color','k','LineWidth',lw,'LineStyle','-')
set(h(3),'Color','r','LineWidth',lw,'LineStyle','-.')
ylabel('XBA2')
xlabel('Time (days)')
%set(gca,'XTick',0:0.5)
set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold')
xlim([5 t(end)])
% 
% for j=[1,2,3,4,5,6]
%     f = strcat(hd1,'_50Fig',num2str(j));
%     saveas(j,f,'tiff');
% end
% cd ..