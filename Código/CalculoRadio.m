clear
clc
syms R1(A,R2);
R1(A,R2) = sqrt(A*R2^2);
DR1A = diff(R1,A);
DR1R2 = diff(R1,R2); 

A = 0.0065; %Relación entre el cuadrado del radio del planeta y el cuadrado del radio de la estrella
dAS = 9e-4; %Error máximo de la relación
dAI = 9e-4; %Error mínimo de la relación

R2 = 0.464; %Radio de la estrella
dRsS = 0.009; %Error máximo del radio de la estrella
dRsI = 0.011; %Error mínimo del radio de la estrella

dR1Ss = double(DR1A(A,R2)*dAS+DR1R2(A,R2)*dRsS); %Error máximo radios solares
dR1Is = double(DR1A(A,R2)*dAI+DR1R2(A,R2)*dRsI); %Error mínimo radios solares

RjRs = 0.10049;
dR1Sj = dR1Ss*(1/RjRs) %Error máximo radios jovianos
dR1Ij = dR1Is*(1/RjRs) %Error mínimo radios jovianos

Rpj = sqrt(A*R2^2)/RjRs %Radio del exoplaneta en radios jovianos