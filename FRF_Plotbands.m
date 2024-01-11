function h=FRF_Plotbands(time,avg,band)
h=figure;
plot(time,xm,'color',[0 0 0],'linewidth',1.5); hold on;
plot(time,band(1,:),'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');
plot(time,band(2,:),'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');