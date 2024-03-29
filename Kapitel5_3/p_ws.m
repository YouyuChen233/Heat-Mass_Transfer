function [pws, Td]=p_ws(pars)
T=pars.T;
n(1)=-2836.5744;
n(2)=-6028.076559;
n(3)=19.54263612;
n(4)=-2.737830188E-2;
n(5)=1.6261698E-5;
n(6)=7.0229056E-10;
n(7)=-1.8680009E-13;
n(8)=2.7150305;
tmp=0;
for k=1:7
    tmp=tmp+n(k)*T^(k-3);
end
tmp=tmp+n(8)*log(T);
pws=exp(tmp);
tmp=0;
tmp1=0;
n(1)=207.98233;
n(2)=-20.156028;
n(3)=0.46778925;
n(4)=-9.2288067E-6;
m(1)=1;
m(2)=-0.13319669;
m(3)=5.6577518E-3;
m(4)=-7.5172865E-5;
for k=1:4
    tmp=tmp+n(k)*(log(pws)^(k-1));
    tmp1=tmp1+m(k)*(log(pws)^(k-1));
end
Td=tmp/tmp1;
end