% Test Prediction
f1=[0.0500    0.1500    0.3000    0.4000    0.5500    0.7000    0.9000    1.1000    1.3500    1.7500    2.2000];
load('SET1.mat');
load('FRF_All.mat');
FRF_=[FRF_DEC;FRF_IC;FRF_EM];
%FRF_=FRF_*2
%BuildFRFsamp
figure(1)

subplot(2,1,2)

PH=angle(FRF_')';
PH=recompactUp(PH);
plot(f1,PH,'linewidth',1.5)
xlim([0.05 2.2]);
ylabel('phase (rad)')

%colororder(K);

subplot(2,1,1)
plot(f1,abs(SET1),'color',[0.5 0.5 0.5])
title('EC FRFs (Body sway)');
%xlabel('freqency (Hz)')
ylabel('gain (°/°)')
hold on;
subplot(2,1,2)
xlim([0.05 2.2]);
PH=angle(SET1')';
PH=recompactUp(PH);
plot(f1,PH,'color',[0.5 0.5 0.5])
xlabel('freqency (Hz)')
hold on;


subplot(2,1,1)
plot(f1,abs(FRF_),'linewidth',1.5)
xlim([0.05 2.2]);
%xlabel('freqency (Hz)')
K = [1 0 1; 0 1 1; 1 0.5 0.5]; % Just random colors.
colororder(K);


subplot(2,1,2)

PH=angle(FRF_')';
PH=recompactUp(PH);
plot(f1,PH,'linewidth',1.5)
xlim([0.05 2.2]);
ylabel('phase (rad)')

colororder(K);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print(gcf,'DataPrediction','-dpdf','-r0')


SET=SET1; % tautologic, it will be generalized in cross validation

%% Bootstrap
[avg,sigma,band,Cp,chist,values] = FRF_PredictionBand(SET,f1,1/22,0.95,4);

%% plots
ntrials=size(SET,1);
ns=441;

yt=zeros(ntrials,ns);
sf=22;

for i=1:ntrials %for the plot
    [x,t]=pseudopulse(SET(i,:),f1,sf);
    yt(i,:)=x;
end


% Histogram
figure(11)
%STAT=sort(STAT);
%H=histogram(STAT,1000,'Normalization','cdf','Edgecolor','none')
plot(values(1:end-1),chist)


xlabel('$C_c$','interpreter','latex')
ylabel('$\max \left( \left| \hat{x}_b - \hat{x}) \right| / \hat{\sigma_b} \right)$','interpreter','latex','Fontsize',11)

hold on;
yline(alpha,'Linewidth',1.5,'Linestyle','-.','Color',[1 0 0],'Label','\alpha = 95%','Fontsize',10)
hold on;
xline(Cp,'Linewidth',1.5,'Linestyle','-.','Color',[1 0 1],'Label','C_p','Fontsize',10,'LabelVerticalAlignment','middle')


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print(gcf,'HistogramPrediction','-dpdf','-r0')



UPPER=band(1,:);
LOWER=band(2,:);
res=zeros(1,length(UPPER));
res(UPPER<0)=UPPER(UPPER<0);
res(LOWER>0)=LOWER(LOWER>0);

st=t(2)-t(1);
ns=length(res);
fd=fft(res); %back to frequency domain
fp=(0:(ns-1))/(st*(ns-1));
RF=zeros(1,11);
for f=1:11
    RF(f)=2*fd(fp==f1(f))/ns;
end


figure(12)

plot(t,avg,'color',[0 0 0],'linewidth',1.5); hold on
plot(t,UPPER,'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');
plot(t,yt,'color',[0.8 0.8 0.8]); hold on
plot(t,avg,'color',[0 0 0],'linewidth',1.5);
plot(t,UPPER,'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');
plot(t,LOWER,'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');
yline(0,'k-.');

legend({'$\hat{x}(t)$','$95\%$','$x_i(t)$'},'interpreter','latex');


xlabel('time (s)');
ylabel('PIR (°)');

title('95% Prediction Band')

% tested samples in FRF_t,yx

ntrials=size(FRF_,1);
ns=441;

yx=zeros(ntrials,ns);

for i=1:ntrials
    [x,t]=pseudopulse(FRF_(i,:),f1,sf);
    yx(i,:)=x;
end
hold on;
plot(t,yx,'linewidth',1.5);
C = [0 1 1; 1 0.5 0.5;1 0 1]; % Just random colors.
colororder(C);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print(gcf,'TimeDomainPrediction','-dpdf','-r0')