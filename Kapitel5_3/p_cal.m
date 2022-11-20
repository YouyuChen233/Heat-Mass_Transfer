function ret=p_cal(pars)
R=pars.R;
T=pars.T;
nu=pars.nu;
yw=pars.yw;
ya=1-yw;
if pars.tag_VC_NS
    [Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww]=VC_NS(pars);
elseif pars.tag_VC_Hy
    [Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww]=VC_Hy(pars);
end
Bm=ya^2*Baa+2*ya*yw*Baw+yw^2*Bww;
Cm=ya^3*Caaa+3*ya^2*yw*Caaw+3*ya*yw^2*Caww+yw^3*Cwww;
ret=R*T/nu*(1+Bm/nu+Cm/nu^2);

end