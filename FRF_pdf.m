function [cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(X,FRFs,phi,sample_time,B)


N=size(FRFs,1); %number of FRFs

Dx=max(1,fix(N/20)); %make it a parameter?


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



STAT=zeros(1,B);
dSTAT=zeros(1,B);

for b=1:B %GENERATE THE HISTOGRAM
    resamp=randi(N,1,N);
    yb=yt(resamp,:);
   
    xb=mean(yb);
    
    es=sort(sum((yb-xb).^2,2)); %"error" for each sample
    et=sum((xt-xb).^2) %error on the test sample
 
    idx=find(es>et,1,'first');
    if isempty(idx)
        idx=N;
    end
    STAT(b)= idx/N;
        
    i1=idx-Dx;
    i2=idx+Dx;
    
    if i1<1
        i1=1;
        i2=1+Dx;
    end
    
    if i2>N
        i2=N;
        i1=N-Dx;
    end
    
    dSTAT(b)=(Dx)/(es(i2)-es(i1));
end

cdf=mean(STAT);
sigma_cdf=std(STAT);

pdf=mean(dSTAT);
sigma_pdf=std(dSTAT);
