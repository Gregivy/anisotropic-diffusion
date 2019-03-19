clear all
close all

I = imread('./imgs/smallpart.png');

I = rgb2gray(I);
I = double(I);

dI = pa_diffusion(I,30,2);

figure
imshow(uint8(I));
figure
imshow(uint8(dI));
%imwrite(uint16(reshape(x,m,n)'),'./article/final/filtered.bmp');