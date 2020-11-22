%% Radial Velocity
clc
clear
time = [2457963.961008, 2457971.918824, 2457988.847288, 2457993.852429, 2458000.953924, 2458002.974541, 2458019.831732, 2458020.821615, 2458028.847305, 2.457940387743000e+06, 2457933.984132, 2.457940980813000e+06];
data = [81.27, -43.37, -54.42, 13.85, 13.49, -39.62, -44.74, -77.22, 56.14, -15.15, 114.31, 5.57];
for i = 1:length(time)
    while time(i) > 2457942.388
        time(i) = time(i)-16.3396;
    end
end
time'
min(time)
max(time)
plot(time, data,'o')

tiempo2 = linspace(2.457928494029000e+06,2.457940980813000e+06,100);

[p, S, mu] = polyfit(time,data,4);
f = polyval(p, tiempo2, [], mu);
hold on
plot(tiempo2, f)

% Propiedades de los ejes
xticks([2457928 2457935 2457942])
xticklabels({'-0.5','x = 0','0.5'})
xlabel('Fase')
ylabel('Velocidad [m/s]')
title('Velocidad radial de la estrella Kepler 643')