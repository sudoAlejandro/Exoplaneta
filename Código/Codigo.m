%% Carga de datos
clc
clear
%---Cargar archivo y mostrar contenido---%
file = 'tess2020050191121-s0022-s0022-0000000138819293-00309_dvt.fits'; 
data = fitsread(file, 'binarytable');
fitsdisp(file)
fase = data{4};
flujo = data{5};
modelo = data{9};
tiempo = data{1};
%% Tiempo de tránsito
%---Limpiar escenario---%
clf
cla
%---Diagrama---%
plot(fase,flujo,'.',fase,modelo,'.')
%---Propiedades de los ejes---%
axis([-0.05 0.05 -9e-3 4e-3])
title('Flujo de la estrella GJ 436');
xlabel('Fase [Días]');
ylabel('Flujo relativo [ppm]');
legend('Datos obtenidos','Modelo');
f = gcf;
f.Color = [1 1 1];
grid on

%% Periodo
%---Limpiar escenario---%
clf
cla
%---Diagrama---%
plot(tiempo,flujo,'.',tiempo,modelo,'-')
%---Propiedades de los ejes---%
axis([1901 1906 -8e-3 4e-3])
title('Dos tránsitos del exoplaneta GJ 436B');
xlabel('Tiempo BJD-2457000');
ylabel('Flujo relativo [ppm]');
f = gcf;
f.Color = [1 1 1];
grid on
%Calculo
Mean = mean(diff(tiempo(find(islocalmin(data{9})))))