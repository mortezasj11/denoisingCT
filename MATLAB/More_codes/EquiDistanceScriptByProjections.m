clc
clear all
load('ProjEnergy_noisy4e4.mat')

downsample=1;
% ProjData

D = 132;
D_sd=180;
L_Detector=0.1;
imgSize   =512;%512
R_obj = imgSize*0.075/2; %15???????????????????????????
%function RecPic = EquiDistanceReconstraction(g,Params)
%g=ProjEnergy(:,:,1)';
g=ProjEnergy(:,:,1)';

%DecfanAngle = asin(1/5)*2.1;
N_detector  = size(g,2);
N_view      = size(g,1);

pi = 3.1415926536;
Deg2Rad = pi/180;

HalfUp_N_detector = round(N_detector/2);
uVector=((1:N_detector) - HalfUp_N_detector)*L_Detector; % du = [-255:256] * 0.45

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RAMP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RAMP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn=2^(nextpow2(N_detector)+1);%2048

Center = nn/2+1  ;  %1025
H_n = zeros(nn,1);  %2048*1

for i=1:nn/2    %1:1024         %for i= nn/2-299 :nn/2         
    H_n(i) = -(1-(-1)^(i-Center))./(2*pi^2*L_Detector^2*(i-Center)^2);
end

H_n(nn/2+1)=1/(4*L_Detector^2);  %1025

for i=nn/2+2:nn   %1026:2048     %for i=nn/2+2:nn/2+300
    
    H_n(i) = -(1-(-1)^(i-Center))./(2*pi^2*L_Detector^2*(i-Center)^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nn >=3*N_detector)
    nn2=nn;
else
    nn2=2*nn;
end

H_n_new = zeros(nn2,1);
H_n_new(1:nn/2) = H_n(nn/2+1:nn);    %  1:1024 <-- 1025:2048
H_n_new(nn2-nn/2+1:nn2) = H_n(1:nn/2);%  1024+1+2048:4096 <-- 1:1024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier of  RAMP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fourier__H_n=fft(H_n_new,nn2);% Ramp filter in fourier domain

%%%%%%%%%%%%%%%%%%%%%%%%%%% FINISHED  RAMP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% FINISHED  RAMP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Weighted g(u,beta)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Weighted g(u,beta)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gWeigthed=zeros(size(g));
for i=1:N_view
   gWeigthed(i,:)= g(i,:).*( D_sd ./ sqrt(D_sd^2 + uVector.^2 ));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fourier of  Weighted g  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = zeros(N_view,N_detector);
www=zeros(N_view,nn2);%??????????????????????????
for i=1:N_view
   fourier__gWeigthed = fft(gWeigthed(i,:),nn2);
   
   Q_temp = real(ifft(fourier__gWeigthed.*fourier__H_n'));  
   www(i,:)=Q_temp;%??????????????????????????
   Q(i,:)= Q_temp(1:N_detector);
   %Q(i,:)= Q_temp(nn2-N_detector+1:nn2);
end
Q = Q .* (uVector(end)-uVector(1))/N_detector;

%%%%%%%%%%%%%%%%%%%%%%%% Finished Weighted g(u,beta)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Finished Weighted g(u,beta)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaX=2*R_obj/imgSize;      deltaY=deltaX;

tempImg=zeros(imgSize,imgSize);
imgCenter=(imgSize+1)/2;
for k=1:downsample:N_view
    k
beta = ((k-1)/N_view)*360*Deg2Rad;
Lambda = pi/2 + beta;
    for imageRow=1:imgSize    
        for imageColumn=1:imgSize
        %%%%%%%%%%%%%%%%%%%%%  x & y   %%%%%%%%%%%%%%%%%%%%%%%%
        x = deltaX * (imageColumn-imgCenter)  ;
        %ie: (2,1) "row 2, column 1" is imgRow=2,imgCol=1
        y = -deltaY * (imageRow-imgCenter)  ;
        
        %%%%%%%%%%%%%%%%%%%%% r & phi %%%%%%%%%%%%%%%%%%%%%%%%%
        r=sqrt(x^2+y^2);
        
        phid = atan2d(y,x);
        phi  = phid*pi/180;  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% u Tilde %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%   first utilde  %%%%%%%%%OOOKKK,nott OK
%       uTilde=D*r*cos(beta-phi) / (D+r*sin(beta-phi));
        
%         %%%%%%   second utilde  %%%%%%%%%OOOKKK,, nott OK
%         xMinusa_LambdaD2 = -(x-D*cos(Lambda))*sin(Lambda) + (y-D*sin(Lambda))*cos(Lambda);
%         xMinusa_LambdaD1 = -(x-D*cos(Lambda))*cos(Lambda) - (y-D*sin(Lambda))*sin(Lambda);
%         uTilde= -D_sd * xMinusa_LambdaD2 / xMinusa_LambdaD1 ;
        

%         %%%%%%%   third utilde  %%%%%%%%%
        gammaTilde = atan( r*cos(beta-phi)/(D+r*sin(beta-phi)) );
        uTilde=tan(gammaTilde)*D_sd;
        
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D Tilde  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        DTilde2=  (x-D*cos(Lambda))^2 + (y-D*sin(Lambda))^2 ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interpolation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
            position = uTilde/L_Detector + N_detector/2;
            
            if (position <= 0 || position >= N_detector+1)
                c1 = 0;
                c2 = 0;
                
            elseif (position >= 1 && position <= N_detector)
            c1 = Q(k ,ceil(position)) *(position-floor(position));
            c2 = Q(k,floor(position)) *(ceil(position)-position) ; 
            
            elseif (position > N_detector && position < N_detector+1 )
                c1 = 0;
                c2 = Q(k,floor(position)) *(ceil(position)-position) ;
            
            elseif (position > 0 && position <1 )
                c1 = Q(k ,ceil(position)) *(position-floor(position));
                c2 = 0;
                
            elseif (ceil(position) == floor(position))
            c1 = Q(k,ceil(position)) ;
            c2 = 0 ;          
            end    
            
            tempImg(imageRow,imageColumn)= tempImg(imageRow,imageColumn) + (c1+c2)*(D_sd^2+uTilde^2)*D /(2*DTilde2*D_sd);
        
        end
    end
end

% numDetector should be based on the number of views
img = tempImg  *  (2*pi/N_view) *downsample  ;

% subplot(133)
% imshow(img,[])
% title('Reconstructed image (512*512)')

subplot(132)
imshow(rot90(img,4),[])
title('Reconstructed image rotated')

imshow(img,[])

%%%% SAVING
NoNoiseFBP_11282016_CH1=img;

codepath=pwd;
filesave = [codepath, '\Results4e3\NoNoiseFBP_11282016_CH1.mat'];
save(filesave,'NoNoiseFBP_11282016_CH1')
