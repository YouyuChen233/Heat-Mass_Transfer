function [psubl, Tf]=p_subl(pars)
T=pars.T;
n(1)=-5866.6426;
n(2)=22.32870244;
n(3)=1.39387003E-2;
n(4)=-3.4262402E-5;
n(5)=2.7040955E-8;
n(6)=0.67063522;
tmp=0;
for k=1:5
    tmp=tmp+n(k)*T^(k-3);
end
tmp=tmp+n(6)*log(T);
psubl=exp(tmp);
tmp=0;
tmp1=0;
n(1)=212.57969;
n(2)=-10.264612;
n(3)=0.14354796;

m(1)=1;
m(2)=-8.2871619E-2;
m(3)=2.3540411E-3;
m(4)=-2.4363951E-5;

for k=1:3
    tmp=tmp+n(k)*(log(psubl)^(k-1));    
end
for i=1:4
    tmp1=tmp1+m(k)*(log(psubl)^(k-1));
end
Tf=tmp/tmp1;
end