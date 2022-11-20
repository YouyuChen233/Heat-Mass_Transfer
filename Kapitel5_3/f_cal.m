function ret = f_cal(yws,pars)
ik=pars.nk52(:,2);
jk=pars.nk52(:,3);
nk=pars.nk52(:,4);
yas=1-yws;
T=pars.T;
R=pars.R;
p=p_cal(pars);
t=84.85;
[pws, Td]=p_ws(pars);

if pars.tag_VC_NS
    [Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww]=VC_NS(pars);
elseif pars.tag_VC_Hy
    [Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww]=VC_Hy(pars);
end
beta_H=0;
n(1)=23.5199;
n(2)=-0.60277;
n(3)=1.1518;
n(4)=-0.12556E-3;
n(5)=0.74217E-6;
n(6)=0.17765E-8;
for k=1:6
    beta_H=beta_H+1e-6*n(k)*t^(k-1);
end
gamma_tau=0;
gamma_tau2=0;
Pi=138600/16.53e6;  %%spezifisch bitte achten
tau=1386/358;       %%spezifisch bitte achten
for k=1:34
    gamma_tau=gamma_tau-nk(k)*ik(k)*(7.1-Pi)^(ik(k)-1)*(tau-1.222)^jk(k);
end
for k=1:34
    gamma_tau2=gamma_tau2+nk(k)*ik(k)*(ik(k)-1)*(7.1-Pi)^(ik(k)-2)*(tau-1.222)^jk(k);
end
nu_ws=R*T/pws*Pi*(gamma_tau);
k_T=-1/p*Pi*gamma_tau2*gamma_tau^(-1);
g=nu_ws*(1-k_T*(p-pws)/2);
ret=g+R*T*log(11-beta_H*yas*p)+Baa*yas^2*p-Bww*(p-pws-yas^2*p)-2*Baw*yas^2*p...
    +Caaa*yas^3*p^2/R/T+Caaw*3*yas^2*(1-2*yas)*p^2-Caww*3*yas^2*yws*p^(2)...
    -Cwww*((1+2*yas)*yws^2*p^2-pws^2)/2/R/T-Baa*Bww*yas^2*(1-3*yas)*yws*p^(2)...
    -Baa*Baw*2*yas^3*(2-3*yas)*p^2/R/T+Bww*Baw*6*yas^2*yws^2*p^2/R/T...
    -Baa^2*3*yas^4*p^2/2/R/T-Baw^2*2*yas^2*yws*(1-3*yas)*p^2/R/T...
    -Bww^2*(pws^2-(1+3*yas)*yws^3*p^2)/2/R/T;
end

