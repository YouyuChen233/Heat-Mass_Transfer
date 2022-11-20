clc;clear;
load("nk.mat");
pars.nk52=nk52;
pars.T=358;
pars.nu=0.83;
pars.R=8.3144;
pars.t=84.85;
y(1)=0.7812;
y(2)=0.0092;
y(3)=0.2096;
M(1)=28.01348;
M(2)=39.948;
M(3)=31.9988;
Mw=18.015268;
Ma=0;
for i=1:3
    Ma=Ma+y(i)*M(i);
end
X=0.35;
pars.yw=X/(Mw/Ma+X);
pars.tag_VC_NS=true;
pars.tag_VC_Hy=not(pars.tag_VC_NS);

