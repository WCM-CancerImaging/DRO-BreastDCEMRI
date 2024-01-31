function [S0,temp]=disp_mask(mask,S0,flag)


temp = zeros(320,320,3);
temp = fill_RGBmask(temp,mask.liver,[255,255,0]);
temp = fill_RGBmask(temp,mask.heart,[255,128,0]);
temp = fill_RGBmask(temp,mask.glandular,[153,255,153]);
% temp = fill_RGBmask(temp,mask.glandular_blood,[240,100,100]);
temp = fill_RGBmask(temp,mask.malignant,[255,0,0]);
temp = fill_RGBmask(temp,mask.benign,[255,0,0]);
temp = fill_RGBmask(temp,mask.vascular,[255,105,180]);
%temp = fill_RGBmask(temp,mask.fat,[229,204,255]);
temp = fill_RGBmask(temp,mask.skin,[0,0,255]);
temp = fill_RGBmask(temp,mask.muscle,[153,0,76]);
temp = im2double(temp);
if flag
% close all;
figure;set(gcf,'Position',[200,200,1400,500]);
subplot(1,2,1);
bck_img = abs(S0);
imshow(bck_img,[0,quantile(bck_img(:),0.98)],'InitialMagnification',1500);
colormap(gray);hold on; freezeColors();
% imshow(bck_img,[],'InitialMagnification',1500)

subplot(1,2,2);
imshow(abs(S0),[],'InitialMagnification',1500);
colormap(gray)
hold on;
freezeColors();
h = imshow(temp/255);
alpha = sum(temp,3)>0;
set(h,'AlphaData',alpha);
end