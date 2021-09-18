clc
clear all
%importacion de los datos del exoploneta
data = fitsread('tess2020050191121-s0022-s0022-0000000138819293-00309_dvt.fits','binarytable')
fitsdisp('tess-s0022-1-1_175.543792_26.708547_10x15_astrocut.fits')

%%Periodo
flux = data{5}
time = data{1}
points = [0 0]

for i = 1:1:2
    fluxi = min(flux)
    ind = find(flux==fluxi)
    points(i)=data{1}(ind)
    for j = 0:1:2
        fluxi = min(flux)
        flux(flux==fluxi)=[]
    end
end
t= points(1)-points(2)
cut = find(data{5}==min(data{5}))
raw_1 = data{5}
raw_2 = data{1}
plot(raw_2,raw_1,'.')
hold on
xlabel('BJD time -2457000')
ylabel('flux (ppm)')
title('BJD time vs Flux')
line = data{1}(find(data{5}==min(data{5})))
for i = 0:1:2
    xline(line-t*i)
    hold on 
end
x1=[]
y1=[]
t
for time = 0:0.01:t
    x1 = [x1 time]
    y1 = [y1 -0.005]
end
for i=0:1:1
    p=points(1)-t*i
    hold on
    plot(-x1+p,y1)
end
legend('','','period','','',sprintf('t=%f',t))
grid('on')