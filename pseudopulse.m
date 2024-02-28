function [x,t]=pseudopulse(y,F,sf)
% manually computes the inverse of the FRF chosing the frequencies (11
% components, maurer's version)
bf=double(gcd(sym(F))); %this defines the period
hf=double(lcm(sym(F)));
%sf=(10*F(end)); %sampling frequency that is at least twice the highest one represented in FRF

step=(1/sf);
t=0:step:(1/bf);
n=length(F);
x=0;
for i=1:n
    %x=x+abs(y(i))*cos(2*pi*F(i)*t+phase(y(i)));
    x=x+real(y(i))*cos(2*pi*F(i)*t)-imag(y(i))*sin(2*pi*F(i)*t);
end
