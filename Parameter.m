clc;clear;
eps=1e-9;
%% Kapitel1
%Dimension
ndim=20;
%Parameter
alpha_feL=0.1;
pho_feL=0.1;
c_p_fel=0.1;
Le=0.87;
G_kond=0.1;
T_feL=0.15*ones(ndim,1);
A=0.1;
y_w_int=0.78;
y=0.1;
delta_h=0.1;
T_int=0.5;
c_p_trL=0.8;
T_aus=0.9;
T_TP=0.7;
X=0.8;
delta_h0=0.3;
c_p_H2O=0.99;
delta_x=0.5;
delta_y=0.4;
pho_w=1;
mu_w=4;
g=9.8;
theta=0.5;
f_feL=0.5*ones(ndim,1);
G_feL=0.5*ones(ndim,1);
d_h=0.1;

%% Kapitel2
T_feL_ein=358;
X_feL_ein=0.35;
dm_feL_ein=1.5;
T_feL_aus=348;
dm_feL_aus=1.4;
p=1.5;
% Molanteile
y_=[0.7812 0.0092 0.2096];
M=[28.01348 39.948 31.9988];
M_D=18.015268;
% Dichte der Feuchtluft
xi_O2=0.245;
xi_N2=0.755;
R_O2=259.8;
R_N2=296.8;
R_D=461;
p_stern=1;
T_stern=540;
Wasserdampf=importdata("Dichte_Wasserdampf.txt");
Hij=importdata("Hij.txt");
H=zeros(5,6);
for i=1:6
    for j=1:7
        for ii=1:size(Hij,1)
            if(abs(Hij(ii,1)-i+1)<eps&&abs(Hij(ii,2)-j+1)<eps)
                H(i,j)=Hij(ii,3);
            end
        end
    end
end
%% 1.1 Dynamische Viskosität der trockenen Luft
M_air=28.9586;
sigma_air=0.36;
b=[0.431 -0.4623 0.08406 0.005341 -0.00331];
Residual=[10.72 0.2 1 0;
    1.122 0.05 4 0;
    0.002019 0.05 4 0;
    -8.876 0.6 1 1;
    -0.02916 3.6 8 1;];
T_c=132.6312;
rho_c=10.4477;
Hi=[1.67752 2.220462 0.6366564 -0.241605];
Therm=importdata("ThermCond.txt");
N=Therm(:,2);
t=Therm(:,3);
d=Therm(:,4);
l=Therm(:,5);
Lij=importdata("Lij.txt");
Lk=[2.443221e-3 1.323095e-2 6.770357e-3 -3.454586e-3 4.096266e-4];
lambda_stern=1e-3;
c_p_D=1.86;
deltah0_v=2500;