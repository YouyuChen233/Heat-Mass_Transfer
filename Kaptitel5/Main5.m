clc;clear;
load("nk.mat");
%% 5.1 %%
R_=8.31451;
T_a=132.6312;
rho_a=10447.7;
T=0.1;%???
rho=50;%???
delta=rho/rho_a;
tau=T_a/T;
dalphaid_dtau=0;
dalphaid_dtau2=0;
for k=1:5
    dalphaid_dtau=dalphaid_dtau+(k-4)*nk51(k,2)*tau^(k-5)+1.5*nk51(6,2)+nk51(7,2)/tau+...
        nk51(8,2)*nk51(11,2)/(exp(nk51(11,2)*tau)-1)+nk51(9,2)*nk51(12,2)/(exp(nk51(12,2)*tau)-1)...
        +nk51(10,2)*nk51(13,2)/(2/3*exp(-nk51(13,2)*tau)+1);
    dalphaid_dtau2=dalphaid_dtau2+(k-4)*(k-5)*nk51(k,2)*tau^(k-6)+0.75*nk51(6,2)*tau^(-0.5)...
        -nk51(7,2)*tau^(-2)-nk51(8,2)*nk51(11,2)^2*exp(nk51(11,2)*tau)/(exp(nk51(11,2)*tau)-1)^2....
        -nk51(9,2)*nk51(12,2)^2*exp(nk51(12,2)*tau)/(exp(nk51(12,2)*tau)-1)^2+...
        2/3*nk51(10,2)*nk51(13,2)^2*exp(-nk51(13,2)*tau)/(2/3*exp(-nk51(13,2)*tau)+1)^2;        
end
dalphares_ddelta=0;
for k=1:10
    dalphares_ddelta=dalphares_ddelta+nk51(k,3)*nk51(k,6)*delta^(nk51(k,3)-1)*tau^(nk51(k,4));
end
for k=11:19
    dalphares_ddelta=dalphares_ddelta+(nk51(k,3)-nk51(k,5)*delta^(nk51(k,5)))*nk51(k,6)...
        *delta^(nk51(k,3)-1)*tau^(nk51(k,4)-1)*exp(-delta^(nk51(k,5)));
end
dalphares_ddelta2=0;
for k=1:10
    dalphares_ddelta2=dalphares_ddelta2+nk51(k,3)*(nk51(k,3)-1)*nk51(k,6)*delta^(nk51(k,3)-2)*tau^(nk51(k,5));
end
for k=11:19
    dalphares_ddelta2=dalphares_ddelta2+((nk51(k,3)-nk51(k,5)*delta^(nk51(k,5)))*(nk51(k,3)-1-nk51(k,5)*delta^(nk51(k,5)))...
        -nk51(k,5)^2*delta^(nk51(k,5)))*nk51(k,6)*delta^(nk51(k,3)-2)*tau^(nk51(k,4))*exp(-delta^(nk51(k,5)));
end
dalphares_dtau=0;
for k=1:10
    dalphares_dtau=dalphares_dtau+nk51(k,4)*nk51(k,6)*delta^(nk51(k,3))*tau^(nk51(k,4)-1);
end
for k=11:19
    dalphares_dtau=dalphares_dtau+nk51(k,4)*nk51(k,6)*delta^(nk51(k,3))*tau^(nk51(k,4)-1)*exp(-delta^(nk51(k,5)));
end
dalphares_dtau2=0;
for k=1:10
    dalphares_dtau2=dalphares_dtau2+nk51(k,4)*(nk51(k,4)-1)*nk51(k,6)*delta^(nk51(k,3))*tau^(nk51(k,4)-2);
end
for k=11:19
    dalphares_dtau2=dalphares_dtau2+nk51(k,4)*(nk51(k,4)-1)*nk51(k,6)*delta^(nk51(k,3))*tau^(nk51(k,4)-2)*exp(-delta^(nk51(k,5)));
end
dalphares_dtau_ddelta=0;
for i=1:10
    dalphares_dtau_ddelta=dalphares_dtau_ddelta+nk51(k,3)*nk51(k,4)*nk51(k,6)*delta^(nk51(k,3)-1)*tau^(nk51(k,4)-1);
end
for i=11:19
    dalphares_dtau_ddelta=dalphares_dtau_ddelta+(nk51(k,3)-nk51(k,5)*delta^(nk51(k,5)))*nk51(k,4)*nk51(k,6)*delta^(nk51(k,3)-1)*tau^(nk51(k,4)-1)*exp(-delta^(nk51(k,5)));
end
Enthalpie51=tau*(dalphaid_dtau+dalphares_dtau)+delta*(dalphares_dtau)+1;
isochore51=-tau^2*(dalphaid_dtau2+dalphares_dtau2);
isobare51=isochore51+(1+delta*dalphares_ddelta-delta*tau*dalphares_dtau_ddelta)^(2)...
/(1+2*delta*dalphares_ddelta+delta^2*dalphares_dtau_ddelta)^2;
%% 5.2 %%
p_stern=16.53e6;
T_stern=1386;

p=1;%???
T=1;%???
Pi=p/p_stern;
tau=T_stern/T;
dgamma_dtau=0;
dgamma_dtau2=0;
for k=1:34
    dgamma_dtau=dgamma_dtau+nk52(k,4)*(7.1-Pi)^(nk52(k,2))*nk52(k,3)*(tau-1.222)^(nk52(k,3)-1);
    dgamma_dtau2=dgamma_dtau2+nk52(k,4)*(7.1-Pi)^(nk52(k,2))*nk52(k,3)*(nk52(k,3)-1)*(tau-1.222)^(nk52(k,3)-2);
end
Enthalpie52=tau*dgamma_dtau;
isobare52=-tau^2*dgamma_dtau2;

%% 5.3/5.4 %%
p=1;%???

p_stern=1e6;

Pi=p/p_stern;

gammaid=0;
gammaid2=0;
for i=1:9
    %gammaid=gammaid+nk53(k,3)*nk53(k,2)*tau^(nk53(k,2)-1);
    %gammaid2=gammaid2+nk53(k,3)*nk53(k,2)*(nk53(k,2)-1)*tau^(nk53(k,2)-2);
end
gammares=0;
gammares2=0;
for i=1:43
    gammares=gammares+nk54(k,4)*Pi^(nk54(k,2))*nk54(k,3)*(tau-0.5)^(nk54(k,3)-1);
    gammares2=gammares2+nk54(k,4)*Pi^(nk54(k,2))*nk54(k,3)*(nk54(k,3)-1)*(tau-0.5)^(nk54(k,3)-2);
end
Enthalpie53_4=tau*(gammaid+gammares);
isobare53_4=-tau^2*(gammaid2+gammares2);

%% 5.7 %%
rho=1;%???
T=1;%???
rho_stern=322;
T_stern=647.096;
R=461.526;
delta=rho/rho_stern;
tau=T_stern/T;
alpha_tau=0;
alpha_tau2=0;
alpha_delta=nk57(1,4)*delta^(-1);
alpha_delta2=-nk57(1,4)*delta^(-2);
alpha_delta_tau=0;
for k=2:40
    alpha_tau=alpha_tau+nk57(k,4)*delta^(nk57(k,2))*nk57(k,3)*tau^(nk57(k,3)-1);
    alpha_tau2=alpha_tau2+nk57(k,4)*delta^(nk57(k,2))*nk57(k,3)*(nk57(k,3)-1)*tau^(nk57(k,3)-2);
    alpha_delta=alpha_delta+nk57(k,4)*nk57(k,2)*delta^(nk57(k,2)-1)*tau^(nk57(k,3));
    alpha_delta2=alpha_delta2+nk57(k,4)*nk57(k,2)*(nk57(k,2)-1)*delta^(nk57(k,2)-2)*tau^(nk57(k,3));
    alpha_delta_tau=alpha_delta_tau+nk57(k,4)*nk57(k,2)*delta^(nk57(k,2)-1)*nk57(k,3)*tau^(nk57(k,3)-1);
end
Enthalpie57=tau*alpha_tau+delta*alpha_delta;
isobare57=-tau^2*alpha_tau2+(delta*alpha_delta-delta*alpha_delta_tau)^2/(2*delta*alpha_delta+delta^2*alpha_delta2);

%% 5.8 %%
p=1;%???
T=1;%???
p_stern=1e6;
T_stern=1;
beta=(p/p_stern)^0.25;
theta=T/T_stern+nk58(2,2)/(T/T_stern-nk58(10,2));
AA=theta^2+nk58(1,2)*theta+nk58(2,2);
BB=nk58(3,2)*theta^2+nk58(4,2)*theta+nk58(5,2);
CC=nk58(6,2)*theta^2+nk58(7,2)*theta+nk58(8,2);
ps_p_stern=(2*CC/(BB^2-4*AA*CC)-BB)^4;

%% 5.9 %%
T=1;%???
p_stern=611.657;
T_stern=273.16;
theta=T/T_stern;
ln_psubl_p_stern=1/theta;
for i=1:3
    ln_psubl_p_stern=ln_psubl_p_stern+ab59(i,2)*theta^(ab59(i,3));
end
clearvars -except ln_psubl_p_stern Enthalpie* iso* 