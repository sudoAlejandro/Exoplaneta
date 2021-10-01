clear
clc
data = xlsread('Radial_Velocity_GJ436.xlsx');
P = 2.6440;
tiem = mod(data(:,1),P);
vel = data(:,2);
plot(tiem,vel,'.');
table = sortrows([tiem, vel]);
