clc
clear
load('ProjReal1_8.mat')
projclean0=ProjReal1_8(:,:,7);

%% FBP

MagFactor = 1;
Proj1 = reshape(fliplr(projclean0),512,1,360);

Reconreal = Recon_offcenter_fullscan_pixsize_MTF_rot(Proj1,-1,1,-1,1,158,255-158,2*0.055,2*0.055,0.05,400,0,360,0);


figure,imshow(flipud(rot90(Reconreal,1)), [])

FBP1_8_1=flipud(rot90(Reconreal,1))/10;






CH7FBPREAL=flipud(rot90(Reconreal,1));
save('CH7FBPREAL.mat','CH7FBPREAL')

ImgFBP(:,:,1)=CH1FBPREAL;
ImgFBP(:,:,2)=CH2FBPREAL;
ImgFBP(:,:,3)=CH3FBPREAL;
ImgFBP(:,:,4)=CH4FBPREAL;
ImgFBP(:,:,5)=CH5FBPREAL;
ImgFBP(:,:,6)=CH6FBPREAL;
ImgFBP(:,:,7)=CH7FBPREAL;
ImgFBP(:,:,8)=CH8FBPREAL;
