%% 
clc
clear all
close all 

%%% Loading image
CHL=1;

load('FBP_Noisy_5e3_8_simpleI.mat')
load('FBP_Orig_5e3_8_simpleI.mat')

z1=FBP_Noisy_5e3_8_simpleI(:,:,CHL);

maxx=max(z1(:));

noisyimg1=z1./maxx*255;
noisyimg=noisyimg1;
noisyimg(noisyimg1<0)=0;

original=FBP_Orig_5e3_8_simpleI(:,:,CHL);
original=original./maxx*255;
original(original<0)=0;
%imshow(original,[])

% sigmaReal(CHL)=sqrt(sum((original(:)-noisyimg(:)).^2)/(512*512)) ;
% [sigmaReal(CHL) noiselevel(CHL)]
% % % % % load('lena.mat')
% % % % % Image1=double(lena);
% % % % % noisyimg=Image1 + 20*randn(size(Image1));

noisyimg=zeros(512,512);
FBP_Noisy_5e3_8_simpleI=FBP_Noisy_5e3_8_simpleI*max(max(FBP_Noisy_5e3_8_simpleI(:,:,1)));
for i=1:8
    noisyimg=noisyimg+FBP_Noisy_5e3_8_simpleI(:,:,i);
end

%% noise estimation
% percentage=10;
% sigmame=NoiseEstimationMoriiExclude( noisyimg, percentage );
% sigma=sigmame*4.90
sigma = noisestYANBO(noisyimg);
%sigmahat1=NoiseEstimationSimple( noisyimg )
%sigma=sigma*2.30;

%%
profile='np';print_to_screen=1;
%[PSNR, y_est,y_hat]        = BM3D([], noisyimg, sigma, 'np', 1);% 'np','lc','high'

[PSNR, y_est,y_hat] = BM3DTestMoey([], noisyimg, sigma, 'np', 1);
%%
pSNR0=psnr(noisyimg/255,original/255) 
peaksnrBASIC = psnr(original/255,double(y_hat)) ;
peaksnrFinal = psnr(original/255,double(y_est)) ;


SSIM0=ssim(original/255,noisyimg/255) ;
SSIM1=ssim(original/255,double(y_hat));
SSIM2=ssim(original/255,double(y_est));

PSNR=[pSNR0 peaksnrBASIC peaksnrFinal]
SSIMM=[SSIM0 SSIM1 SSIM2]
% 
% figure(3)
% imshow([noisyimg original y_est*255],[]);title('Authors 1st and second step')
% figure(4)
imshow([noisyimg (y_hat)*255 y_est*255],[0 150],'border','tight');title('Authors 1st and second step')
% figure(1)
% subplot(121);imshow(z,[])
% subplot(122);imshow(y_hat,[]);title('y after 1st step')

% figure(2)
% subplot(121)
% imshow(y_hat,[]);title('Authors 1st step')
% subplot(122);imshow(y_est);title('Authors second step')


