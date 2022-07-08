clc; clear;
Parameter;
i=1;
% Teil 1
beta_feL=alpha_feL/pho_feL/c_p_fel*Le^(2/3); %Gl 7
dm_kond=pho_feL*beta_feL*A*log((1-y_w_int)/(1-y));% y not a array, Gl 6
G_kond=dm_kond/delta_x/delta_y; %Gl 8.1
Co=c_p_fel*G_kond/alpha_feL; %Gl 8
T_m_feL=0.5*(T_feL(i)+T_feL(i+1)); %Gl 9 Array??

dQ_feL=dm_kond*delta_h+alpha_feL*A*Co/(1-exp(1)^(-Co))*(T_m_feL-T_int); %Gl 5
%Init Array
dm_feL=zeros(ndim,1);
dm_w=zeros(ndim,1);
h_feL=zeros(ndim,1);
h_w=zeros(ndim,1);
%rekursive berechnen
for i=1:ndim-1
    dm_feL(i+1)=dm_feL(i)-dm_kond;
    dm_w(i+1)=dm_w(i)+dm_kond;
    h_feL(i+1)=c_p_trL*(T_aus-T_TP)+X*(delta_h0+c_p_H2O*(T_aus-T_TP));
    h_w(i+1)=c_p_H2O*(T_aus-T_TP);
end

% Teil 2
pp=[-pho_w^2/3/mu_w*g*sin(theta) pho_w/2/mu_w*tau_delta 0 -dm_w(ndim)];
delta_film=roots(pp);


clearvars -except dm_feL dm_w h_w dm_kond delta_film
