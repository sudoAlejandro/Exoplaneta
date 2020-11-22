clear
clc
bindata = fitsread('3.fits','binarytable');
fitsdisp('3.fits')
plot(bindata{4},bindata{8},'.')
%%
% clc
% paso = []
% tiempos = []
% for i = 0:87
%     time = bindata{4}(473+800*i:1272+800*i);
%     values = bindata{8}(473+800*i:1272+800*i);
%     paso = [paso values];
%     tiempos = [tiempos time];
% end
% 
% for i = 1:87
%    hold on
%    plot(tiempos(:,i),paso(:,i),'.') 
%    pause(0.01)
% end

%% quitar NaN
clc
valores = bindata{8};
aeliminar = find(isnan(valores));
times = bindata{4};
times(aeliminar) = [];
valores = rmmissing(valores);
% plot(times, valores, '.');

%% Recorte
clc
recorteT = find(times > 1 | times <-1);
recorteT = []
times(recorteT) = [];
valores(recorteT) = [];
valores = valores+1
plot(times,valores,'.')
hold on
minline = -1.35e-3+1;
maxline = -1.98e-3+1;
meanline = (minline+maxline)/2
line([-1,1],[minline,minline],'color','red')
line([-1,1],[maxline,maxline],'color','red')
line([-1,1],[meanline,meanline],'color','green')
title('Curva de luz kepler 643B')
xlabel('Fase')
ylabel('Brillo relativo')
%% Promedio
clc
puntos = 100;
minimum = min(times)
maximum = max(times)
rangos = linspace(minimum, maximum, puntos);

intervalos = zeros(1,puntos-1);
for i = 1:puntos-1
%     tiemposS = find(times > rangos(i) & times < rangos(i+1))
%    rango = rangos(i):rangos(i+1);
   intervalos(i) = mean(valores(find(times > rangos(i) & times < rangos(i+1))));
end
% tiempos = times(floor(rangos(1:end-1)));
% plot(tiempos,intervalos,'.');
plot(intervalos)
min(intervalos)
%% Conversor radio
clc
rs = sqrt(0.00152)*2.52;
rs = 0.1027
km = 696340*rs;
rj = km/69911


