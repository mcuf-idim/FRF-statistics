function [avg,sigma,band,Cc,chist,values] = FRF_ConfidenceBandDifference(FRF1,FRF2,phi,sample_time,alpha,B,Bs)
%[AVG,SIGMA,band,CC,CHIST,VALUES] =
%FRF_CONFIDENCEBANDDIFFERENCE(FRFS1,FRFS2,PHI,SAMPLE_TIME,B)
% Confidence bands on the difference between the means of two groups FRF1
% andd FRF2.
%
% B is the number of bootstrap repetitions
% Bs is the number of bootstrap repetitions used to estimate STD
%
% avg is the difference between average PIRs of the groups, sigma is the
% measure of variation ?ˆx(t), band a two row matrix with the boundaries of
% the band. Cp is the threshold constant obtained by the bootstrap.  FRFS
% is a matrix where each row represents a FRF of the set, phi is the vector
% of frequencies  and SAMPLE_TIME is the sample time of the PIRs. Chist is
% a vector representing the cumulative histogram for the values  returned
% in VALUES.

sf=1/sample_time;

N1=size(FRF1,1); %number of FRFs

x1=FRF_pseudoimpulse(FRF1(1,:),phi,sf);
ns=length(x1);
y1=zeros(N1,ns);
y1(1,:)=x1;

for i=2:N1
    [x,t]=FRF_pseudoimpulse(FRF1(i,:),phi,sf);
    y1(i,:)=x;
end

N2=size(FRF2,1); %number of FRFs

x1=FRF_pseudoimpulse(FRF2(1,:),phi,sf);
ns=length(x1);
y2=zeros(N1,ns);
y2(1,:)=x1;

for i=2:N2
    [x,t]=FRF_pseudoimpulse(FRF2(i,:),phi,sf);
    y2(i,:)=x;
end



xm=mean(y1)-mean(y2);

sSTAT=zeros(Bs,ns);
for b=1:Bs
    resamp1=randi(N1,1,N1);
    resamp2=randi(N2,1,N2);
    yb1=y1(resamp1,:);
    yb2=y2(resamp2,:);
    
    sSTAT(b,:)=mean(yb1)-mean(yb2);
end
sx=std(sSTAT);


STAT=zeros(1,B);


for b=1:B %GENERATE THE HISTOGRAM
    resamp1=randi(N1,1,N1);
    resamp2=randi(N2,1,N2);
    yb1=y1(resamp1,:);
    yb2=y2(resamp2,:);
    
    xb=mean(yb1)-mean(yb2);
    
    for b2=1:Bs
        resampb1=randi(N1,1,N1);
        resampb2=randi(N2,1,N2);
        ybb1=yb1(resampb1,:);
        ybb2=yb2(resampb2,:);
        
        sSTAT(b2,:)=mean(ybb1)-mean(ybb2);
    end
    
    sb=std(sSTAT);
    db=max(abs(xm-xb)./sb);
    STAT(b)=db;
end


%% Histogram

STAT=sort(STAT);
%H=histogram(STAT,1000,'Normalization','cdf','Edgecolor','none')
[chist, values] = histcounts(STAT,1000,'Normalization','cdf');
Cc=values(find(chist>alpha,1,'first'));

avg=xm;
sigma=sx;

band=[avg+Cc*sigma;avg-Cc*sigma];

end

