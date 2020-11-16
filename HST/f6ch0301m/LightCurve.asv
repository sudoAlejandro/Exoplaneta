clear
clc
bindata = fitsread('3.fits','binarytable');
fitsdisp('3.fits')
subplot(2,2,1)
plot(bindata{4},bindata{5},'.')
subplot(2,2,2)
plot(bindata{4},bindata{9})
subplot(2,2,3)
plot(bindata{4},bindata{7},'.')
subplot(2,2,4)
plot(bindata{4},bindata{8},'.')
%%
plot(bindata{4})
%%
clc
paso = []
tiempos = []
for i = 0:87
    time = bindata{4}(473+800*i:1272+800*i);
    values = bindata{8}(473+800*i:1272+800*i);
    paso = [paso values];
    tiempos = [tiempos time];
end

for i = 1:87
%    hold on
   plot(tiempos(:,i),paso(:,i),'.') 
   pause(0.01)
end
%% Promediando
clc
promedio = []
paso(:,3)
for i = 1:length(paso(:,3))
    disp("a")
end
%% Conversor radio
clc
rs = sqrt(0.0015)*2.52
km = 696340*rs
rj = km/69911


