clc;clear;
Parameter;
[pws, Td]=p_ws(pars);
p=p_cal(pars);
yws=1*pws/p_cal(pars);
for iteration=1:3
    f=f_cal(yws,pars);
    yws=f*pws/p;
end

