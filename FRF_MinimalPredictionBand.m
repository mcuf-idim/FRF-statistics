function [avg,sigma,band,Cp,chist,values,alpha] = FRF_MinimalPredictionBand(X,FRFs,phi,sample_time,B)
%[AVG,SIGMA,VOL,CP,CHIST,VALUES,ALPHA] =
%FRF_MINIMALPREDICTIONBAND(X,FRFS,PHI,SAMPLE_TIME,B)
% Computes the minimum Cp that includes the tested FRF X and the empirical
% ALPHA associated with it.
% where avg is average PIR, sigma is the measure of variation ?ˆx(t), band
% a two row matrix with the boundaries of the band. Cp is the threshold
% constant obtained by the bootstrap.  FRFS is a matrix where each row
% represents a FRF of the set, phi is the vector of frequencies  and
% SAMPLE_TIME is the sample time of the PIRs. Chist is a vector
% representing the cumulative histogram for the values  returned in VALUES.

N=size(FRFs,1); %number of FRFs
ns=441;

yt=zeros(N,ns);
sf=1/sample_time;

xt=FRF_pseudoimpulse(X,phi,sf);

for i=1:N
    [x,t]=FRF_pseudoimpulse(FRFs(i,:),phi,sf);
    yt(i,:)=x;
end

    xm=mean(yt);
    sx=std(yt);


STAT=zeros(1,B*N);

s=1;
for b=1:B %GENERATE THE HISTOGRAM
    resamp=randi(N,1,N);
    yb=yt(resamp,:);
    
    xb=mean(yb);
    sb=std(yb);
    for n=1:N
     db=max(abs(yt(n,:)-xb)./sb);
     STAT(s)=db;
     s=s+1;
    end
end


%% Histogram

STAT=sort(STAT);
%H=histogram(STAT,1000,'Normalization','cdf','Edgecolor','none')
[chist, values] = histcounts(STAT,1000,'Normalization','cdf'); 
%Cp=values(find(chist>alpha,1,'first'));

[Cp,idx]=max(abs(xt-xm)./sx);

alpha = chist(find(values>Cp,1,'first'));

avg = xm;
sigma = sx;

band=[avg+Cp*sigma;avg-Cp*sigma];
end

