clear
clc
%% Semiamplitud
a = 0.0291;     %Semieje mayor [ua]
T = 2.6440;     %Periodo [dias]
K = 9.8702e-6;  %Semiamplitud [ua/d]
e = 0.13827;    %Excentricidad
G = 2.9597e-4;   %CGU [UA3/Msd^2]
i = deg2rad(86.44);      %Ángulo de inclinación [grados]
Ms = 0.4412;    %Masa estelar [Masas solares]

% A = (T*K^3*(1-e^2)^(2/3))/(2*pi*G);
% Mp = roots([sin(i)^3 -A -2*A*Ms -A*Ms^2]);

Mp = (K*T^(1/3)*(1-e^2)^(1/2)*Ms^(2/3))/(2*pi*G)^(1/3);
Mp*1047.57