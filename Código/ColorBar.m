clear
clc
heigth = 50;
waveLengths = 400:699;
redLayer = zeros(heigth,300);
greenLayer = zeros(heigth,300);
blueLayer = zeros(heigth,300);

for i = waveLengths
   RGB = wavelength2color(i);
   redLayer(:,i-399) = RGB(1);
   greenLayer(:,i-399) = RGB(2);
   blueLayer(:,i-399) = RGB(3);
end

colorBar = redLayer;
colorBar(:,:,2) = greenLayer;
colorBar(:,:,3) = blueLayer;

imshow(colorBar,'InitialMagnification',1000, 'Interpolation', 'bilinear')