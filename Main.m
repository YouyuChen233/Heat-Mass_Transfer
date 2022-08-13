clc; clear;
Parameter;
i=1;
%% Kapitel1
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
tau_delta=zeros(ndim,1);
delta_p_feL=zeros(ndim,1);
delta_p_acc=zeros(ndim,1);
delta_p_stat=zeros(ndim,1);
delta_p_r=zeros(ndim,1);
p_feL=zeros(ndim,1);
pho_feL=0.5*ones(ndim,1); %%%%% im Teil2 -> Array
for i=1:ndim-1
    f_feL_m=0.5*(f_feL(i)+f_feL(i+1));
    G_feL_m=0.5*(G_feL(i)+G_feL(i+1));
    pho_feL_m=0.5*(pho_feL(i)+pho_feL(i+1));

    delta_p_r(i)=2*f_feL_m*G_feL_m^2*delta_x/pho_feL_m/d_h;
    delta_p_stat(i)=pho_feL_m*g*sin(theta);
    delta_p_acc(i)=G_feL(i+1)^2/pho_feL(i+1)-G_feL(i)^2/pho_feL(i);

    delta_p_feL(i)=delta_p_r(i)+delta_p_stat(i)+delta_p_acc(i);
    tau_delta(i)=d_h*0.25/delta_x*delta_p_feL(i);
    p_feL(i+1)=p_feL(i)-delta_p_feL(i);
end


pp=zeros(ndim,4);
delta_film=zeros(ndim,3);
for i=1:ndim-1
    pp(i,:)=[-pho_w^2/3/mu_w*g*sin(theta) pho_w/2/mu_w*tau_delta(i) 0 -0.5*(dm_w(i)+dm_w(i+1))];
    delta_film(i,:)=roots(pp(i,:))';
end

%% Kapitel2
% feuchte Luft
X_feL_aus=dm_feL_aus/dm_feL_ein*(1+X_feL_ein)-1;
X_feL_m=0.5*(X_feL_aus+X_feL_ein);
T_feL_m=0.5*(T_feL_aus+T_feL_ein);

% Molanteile in feuchten Luft
M_trL=0;
for i=1:length(y_)
    M_trL=M_trL+y_(i)*M(i); % Molmasse von trockenen Luft
end
y_D=X_feL_m/M_D/(1/M_trL+X_feL_m/M_D); % Molanteil des Wassers in feuchten Luft
y_trL=1-y_D; % Molanteil der trockenen Luft in feuchten Luft

% Dichte der Feuchtluft
rho_feL=p*(0.5*X_feL_ein+X_feL_aus)+xi_O2+xi_N2;
rho_feL=rho_feL/T_feL_ein/(xi_O2*R_O2+xi_N2*R_N2+0.5*(X_feL_ein+X_feL_aus)*R_D);

% Dichte der trockenen Luft
rho_trL=rho_feL/(1+X_feL_m);

% Dichte von Wasserdampf
gamma_pi_o=p_stern/p;
gamma_pi_r=0;
for i=1:size(Wasserdampf,1)
    gamma_pi_r=gamma_pi_r+Wasserdampf(i,4)*Wasserdampf(i,2)*(p/p_stern)^(Wasserdampf(i,2)-1)*(T_stern/T_feL_m-0.5)^(Wasserdampf(i,3));
end
nu_D=R_D*T_feL_m*(gamma_pi_r+gamma_pi_o)/p_stern;
rho_D=1/nu_D;

%% 1.Dynamische Viskosität
% 1.1 Dynamische Viskosität der trockenen Luft
tmp=0;
for i=1:length(b)
    tmp=tmp+b(i)*log(T_stern)^(i-1);
end
omega=exp(tmp);
eta_0=0.0266958*sqrt(M_air*T_feL_m)/(sigma_air^2*omega);
eta_r=0;
tau=T_c/T_feL_m;
delta=rho_trL/rho_c;
for i=1:size(Residual,1)
    if(abs(Residual(i,4))<eps)
        gamma=0;
    else
        gamma=1;
    end
    eta_r=eta_r+Residual(i,2)*tau^(Residual(i,2))*delta^(Residual(i,3))*exp(-gamma*delta^(Residual(i,4)));
end
eta_trL=eta_0+eta_r;
% 1.2Dynamische Viskosität von Wasser und Wasserdampf
T_stern=647.096;
rho_stern=322;
eta_stern=1e-6;
T_=T_feL_m/T_stern;
rho_=rho_D/rho_stern;
eta2=1;
tmp=0;
for i=1:6
    for j=1:7
        tmp=tmp+ (1/T_-1)^(i-1)*H(i,j)*(rho_-1)^(j-1);
    end
end
eta1=exp(rho_*tmp);
tmp=0;
for i=1:4
    tmp=tmp+Hi(i)/T_^(i-1);
end
eta0=100*sqrt(T_)/tmp;
eta_=eta0+eta1+eta2;
eta_D=eta_*eta_stern;

% Dynamische Viskosität von feuchten Luft
phi_LD=(1+(eta_trL/eta_D)^0.5*(M_D/M_trL)^0.25)^2/(8*(1+M_trL/M_D))^0.5;
phi_DL=(1+(eta_D/eta_trL)^0.5*(M_trL/M_D)^0.25)^2/(8*(1+M_D/M_trL))^0.5;
eta_feL=y_trL*eta_trL/(y_trL+y_D*phi_LD)+y_D*eta_D/(y_D+y_trL*phi_DL);

%% 2.Wärmeleitfähigkeit
% 2.1 Wärmeleitfähigkeit der trockenen Luft
lambda_0=N(1)*eta_0+N(2)*tau^(t(2))+N(3)*tau^(t(3));
lambda_r=0;
lambda_c=0;
tau=T_c/T_feL_m;
delta=rho_trL/rho_c;
for i=4:size(Therm,1)
    if(abs(l(i))<eps)
        gamma=0;
    else
        gamma=1;
    end
    lambda_r=lambda_r+N(i)*tau^(t(i))*delta^(d(i))*exp(-gamma*delta^(l(i)));
end

lambda_trL=lambda_0+lambda_r+lambda_c;

% 2.2 Wärmeleitfähigkeit von Wasser und Wasserdampf

tmp=0;
for i=1:5
    tmp=tmp+Lk(i)/T_^(i-1);
end
lambda0_=sqrt(T_)/tmp;
tmp=0;
tmp1=0;
for j=1:6
    tmp=tmp+Lij(i,j)*(rho_-1)^(j-1);
end
for i=1:5
    tmp1=tmp1+tmp*(1/T_-1)^(i-1);
end
lambda1_=exp(rho_*tmp1);
lambda2_=0;
lambda_=lambda0_*lambda1_+lambda2_;
lambda_D=lambda_*lambda_stern;
% 2.3 Wärmeleitfähigkeit von feuchten Luft
lambda_feL=y_trL*lambda_trL/(y_trL+y_D*phi_LD)+y_D*lambda_D/(y_D+y_trL*phi_DL);
%% 3 pezifische Enthalpie der feuchten Luft
h_feL_ein=c_p_trL*T_feL_ein+X_feL_ein*(deltah0_v+c_p_D*T_feL_ein);
h_feL_aus=c_p_trL*T_feL_aus+X_feL_aus*(deltah0_v+c_p_D*T_feL_aus);
clearvars -except dm_feL dm_w h_w dm_kond delta_film eta_feL eta_D eta_trL lambda_trL h_feL_ein h_feL_aus lambda_feL
