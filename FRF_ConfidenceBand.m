function [avg,sigma,band,Cc,chist,values] = FRF_ConfidenceBand(FRFs,phi,sample_time,alpha,B)
%[AVG,SIGMA,band,CC,CHIST,VALUES] =
%FRF_CONFIDENCEBAND(FRFS,PHI,SAMPLE_TIME,B)
% where avg is average PIR, sigma is the measure of variation ?ˆx(t), band
% a two row matrix with the boundaries of the band. Cp is the threshold
% constant obtained by the bootstrap.  FRFS is a matrix where each row
% represents a FRF of the set, phi is the vector of frequencies  and
% SAMPLE_TIME is the sample time of the PIRs. Chist is a vector
% representing the cumulative histogram for the values  returned in VALUES.

N=size(FRFs,1); %number of FRFs

sf=1/sample_time;

x1=FRF_pseudoimpulse(FRFs(1,:),phi,sf);
ns=length(x1);
yt=zeros(N,ns);
yt(1,:)=x1;

for i=2:N
    [x,t]=FRF_pseudoimpulse(FRFs(i,:),phi,sf);
    yt(i,:)=x;
end

    xm=mean(yt);
    sx=std(yt);


STAT=zeros(1,B);


for b=1:B %GENERATE THE HISTOGRAM
    resamp=randi(N,1,N);
    yb=yt(resamp,:);
    
    xb=mean(yb);
    sb=std(yb);
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

