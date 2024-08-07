function [cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(varargin)
% [cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(X,FRFs,phi,sample_time,B,metric)
%
% estimates the cumulative densitiy function and the density function
% associated to the sample X and a STD on their estimation. 
% the input METRIC defines the measure used to define the distance
% between X and the mean of the sample FRFs. By default it is the sum of
% squared residuals
% distance, if METRIC is a function handle it is applied directly. The
% following strings can be specified:
% - 'squared' sum of squared residuals
% - 'max' maximum difference between two samples
% 

if nargin==5
   metric=@(x,y) sum((x-y).^2,2);
elseif nargin==6
    metric = varargin{6};
    if isa(metric,'function_handle')
       disp(''); % so far do nothing and use it straightforward 
    elseif isa(metric,'string') || isa(metric,'char')
        if metric== "squared"
            metric=@(x,y) sum((x-y).^2,2);
        elseif metric== "max"
            metric=@(x,y) max(abs(x-y),[],2);
        else
            error([metric,' is not a valid metric']);
        end
    else
        error('METRIC must be a string or a function handle');
    end
    
else
    error('Input arguments must be 5 or 6');
end
X=varargin{1};
FRFs=varargin{2};
phi=varargin{3};
sample_time=varargin{4};
B=varargin{5};


N=size(FRFs,1); %number of FRFs

Ds=max(1,fix(N/20)); %make it a parameter?



sf=1/sample_time;

xt=FRF_pseudoimpulse(X,phi,sf);

ns=length(xt); %generalize!

yt=zeros(N,ns);

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
    
    %es=sort(sum((yb-xb).^2,2)); %"error" for each sample
    %et=sum((xt-xb).^2); %error on the test sample
    
    es=sort(metric(yb,xb));
    et=metric(xt,xb);
    
    idx=find(es>et,1,'first');
    if isempty(idx)
        idx=N;
    end
    STAT(b)= idx/N;
        
    i1=idx-Ds;
    i2=idx+Ds;
    
    if i1<1
        i1=1;
        i2=1+Ds;
    end
    
    if i2>N
        i2=N;
        i1=N-Ds;
    end
    
    dSTAT(b)=(Ds)/(N*(es(i2)-es(i1)));
end

cdf=mean(STAT);
sigma_cdf=std(STAT);

pdf=mean(dSTAT);
sigma_pdf=std(dSTAT);
