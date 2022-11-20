function [Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww]=VC_NS(pars)
T=pars.T;
Baa=0;
Caaa=0;
Bww=0;
Cwww=0;
Baw=0;
Caaw=0;
Caww=0;
n(1)=31.831763;
n(2)=-719.51195;
n(3)=-6538137.0;
n(4)=1.5929828E9;
n(5)=-2.5588842E11;
n(6)=2.2300382E13;
n(7)=-8.2793465E14;
for k=1:7
    Baa=Baa+n(k)*T^(1-k);
end
n(1)=1297.5378;
n(2)=46021.328;
n(3)=40813154.0;
n(4)=-3.202391E9;
n(5)=2.2964785E11;
n(6)=-5.3683467E12;
n(7)=-2.1183915E14;
for k=1:7
    Caaa=Caaa++n(k)*T^(1-k);
end
n(1)=-4965.8164;
n(2)=1.81918860E7;
n(3)=-2.92013995058E10;
n(4)=2.7032989E13;
n(5)=-1.6045262E16;
n(6)=6.3750397E18;
n(7)=-1.7206027E21;
n(8)=3.1222306E23;
n(9)=-3.6643847E25;
n(10)=2.5256562E27;
n(11)=-7.929348E28;
for k=1:11
    Bww=Bww+n(k)*T^(1-k);
end
n(1)=-6.566276606;
n(2)=0.3894679516;
n(3)=-0.0034281020537;
n(4)=1.333924918E-5;
n(5)=-2.726404078E-8;
n(6)=2.839369136E-11;
n(7)=-1.189114330E-14;
for k=1:7
    Cwww=Cwww+n(k)*T^(1-k);
end
Cwww=-1e6*exp(Cwww);
n(1)=32.366097;
n(2)=-14113.8;
n(3)=-1244535;
n(4)=0;
n(5)=-0.2348789E10;
for k=1:5
    Baw=Baw+n(k)*T^(1-k);
end
n(1)=482.737;
n(2)=105678;
n(3)=-6.56394E7;
n(4)=2.94442E10;
n(5)=-3.19317E12;
for k=1:5
    Caaw=Caaw+n(k)*T^(1-k);
end
n(1)=-10.72887;
n(2)=3478.04;
n(3)=-383383;
n(4)=33406000;
for k=1:4
    Caww=Caww+n(k)*T^(1-k);
end
Caww=-1e6*exp(Caww);
end