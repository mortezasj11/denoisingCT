clear
clc

% bi manie in m file alan
 


Im = phantom('Modified Shepp-Logan',512);

params.NumofView     = 360;
params.NumofBin      = 512;       % number of detector bins
params.reconsize     = 400;
params.binsize       = 0.11;
params.Dsource2centr = 158;
params.Dsource2detec = 255;

doff = 0;%8.25
params.binshift = doff*params.binsize*0;  
ReconSize=512;
params.pixelsize =0.05;
paramsones = params;

%imshow(Im,[])
paramsones.im = Im;
Proj = projdd(paramsones);
%imshow(Proj,[])

PhantomProj=Proj;

%corrupting the detector bins by
DetMaxIndex=512-86;
DetStartIndex=86;

indexx=randi(DetMaxIndex,30,1)+DetStartIndex;

for i=1:length(indexx)
    Plus1OrMinus1=real(exp(sqrt(-1)*randi(2)*pi));
    PhantomProj(indexx(i),:)=PhantomProj(indexx(i),:)+ Plus1OrMinus1*( (4*randn/100)*PhantomProj(indexx(i),:) );
end
% figure(100)
% imshow(PhantomProj,[],'border','tight')


Proj=PhantomProj;
%save('ProjPhantom512N4.mat','Proj','indexx')
%save('ProjPhantom512N2.mat','Proj','indexx')
% load('ProjPhantom512N2.mat')

load('ProjPhantom512N4.mat')
figure(100)
imshow(Proj,[],'border','tight')
%%        
params.NumofView     = 360;
params.NumofBin      = 512;       % number of detector bins
params.reconsize     = 400;
params.binsize       = 0.11;
params.proj          = Proj;
params.Dsource2centr = 158;
params.Dsource2detec = 255;

doff = 8.25;
params.binshift = doff*params.binsize*0;  

ReconSize=512;

params.pixelsize =0.05;

        
%% weighting
projBadDetecWeightMC = ones(size(Proj));
projParkerWeightMC = ones(size(Proj));
projTotalWeightMC = projBadDetecWeightMC.*projParkerWeightMC;
relaxFactor = 1;         % relaxition facter, typically set to 1
projTotalWeightMC = relaxFactor*projTotalWeightMC;      % relaxition factor

reconsize = ReconSize;

imgInitMC = zeros(reconsize,reconsize);   % initalized image

params.imgInitMC = imgInitMC;
params.reconsize = reconsize;

%%
OS=32;
% Reserve memory for the auxillary variables
NumofView=params.NumofView;
reconsize=params.reconsize;
dxMC = zeros(reconsize,reconsize);      %512*512
dyMC = zeros(reconsize,reconsize);      %512*512
dTNN1MC = zeros(reconsize,reconsize);   %512*512
bxMC = zeros(reconsize,reconsize);      %512*512
byMC = zeros(reconsize,reconsize);      %512*512
bTNN1MC = zeros(reconsize,reconsize);   %512*512


paramsones = params;%%%%%%%%%%%%%%%%%%%???????????????
paramsones.im = ones(reconsize);
projones = projdd(paramsones);
paramsones.reconsize = reconsize;
imgRaypixsumMC = zeros(reconsize, reconsize);

% backprojection
paramsones.proj = projTotalWeightMC(:, :).*projones;%ones 512*512*8
bprojones = bprojdd(paramsones);
imgRaypixsumMC(:,:) = bprojones;

RecIm_M=zeros(reconsize,reconsize);

CorrectionProjection=zeros(size(Proj));

% OS_SART with split Bregman
muTV = 0;%0.02
imgNoiseStd=0.0006;%0.0012
%lambdaTV = muTV*imgNoiseStd ;

G_matrix=Proj;
max_iteration=12;

%% ONLY SART
RecIm_M0=RecIm_M;
params.projMC = Proj;
% iteration of each subset
for SART=1:10
    fprintf('SART=%d \n',SART);
    for iOS = 1:OS

        indViewOS = iOS:OS:NumofView;
        gradFidMC = - OSSARTCoreMCwdd(params, RecIm_M0, projTotalWeightMC, indViewOS);%512*512*8

        % OS-SART update
        RecIm_M0 = RecIm_M0 - NumofView/length(indViewOS)*gradFidMC./imgRaypixsumMC;
%         imshow(- NumofView/length(indViewOS)*gradFidMC./imgRaypixsumMC,[])
        % nonnegtivity
        RecIm_M0(RecIm_M0<0) = 0;

        % min TV & low rank

        RecIm_M0 = gsU2TNN1(RecIm_M0, RecIm_M0, dxMC, dyMC,dTNN1MC, bxMC, byMC, bTNN1MC, muTV, 0, reconsize, reconsize);

        [dxMC, dyMC] = gsSpace2(RecIm_M0, dxMC, dyMC, bxMC, byMC, imgNoiseStd, reconsize, reconsize);%[B_1, B_2]=

        bxMC = bregmanX(dxMC, RecIm_M0, bxMC, reconsize, reconsize);%  bxMC ~ W_1  , dxMC ~ B_1
        byMC = bregmanY(dyMC, RecIm_M0, byMC, reconsize, reconsize);%  byMC ~ W_2  , dyMC ~ B_2

        typeflag = 'Thresh';
        dTNN1MC = RecIm_M0+bTNN1MC;%dTNN1MC~D
        bTNN1MC = bTNN1MC + RecIm_M0 - dTNN1MC;% bTNN1MC~V            
    end
end
index=0;
RecIm_M0(RecIm_M0<0)=0;
figure(200),imshow(RecIm_M0,[0 .6],'border','tight')
title('SART Only')
%%
RecIm_M=zeros(reconsize,reconsize);
index=0;
for kIter = 1:max_iteration

    G_matrix=G_matrix+ CorrectionProjection;
    %G_matrix=reshape(G_matrix,[448 1024]);%???
    params.projMC = G_matrix;
    % iteration of each subset
    for SART=1:3
        fprintf('itter=%d , SART=%d    ,   zero_Coeffs= %d / %d \n',kIter,SART, index,365-68);
        
        for iOS = 1:OS        
            indViewOS = iOS:OS:NumofView;
            gradFidMC = - OSSARTCoreMCwdd(params, RecIm_M, projTotalWeightMC, indViewOS);%512*512*8

            % OS-SART update
            RecIm_M = RecIm_M - NumofView/length(indViewOS)*gradFidMC./imgRaypixsumMC;
    %         imshow(- NumofView/length(indViewOS)*gradFidMC./imgRaypixsumMC,[])
            % nonnegtivity
            RecIm_M(RecIm_M<0) = 0;

            % min TV & low rank

            RecIm_M = gsU2TNN1(RecIm_M, RecIm_M, dxMC, dyMC,dTNN1MC, bxMC, byMC, bTNN1MC, muTV, 0, reconsize, reconsize);

            [dxMC, dyMC] = gsSpace2(RecIm_M, dxMC, dyMC, bxMC, byMC, imgNoiseStd, reconsize, reconsize);%[B_1, B_2]=

            bxMC = bregmanX(dxMC, RecIm_M, bxMC, reconsize, reconsize);%  bxMC ~ W_1  , dxMC ~ B_1
            byMC = bregmanY(dyMC, RecIm_M, byMC, reconsize, reconsize);%  byMC ~ W_2  , dyMC ~ B_2

            typeflag = 'Thresh';
            dTNN1MC = RecIm_M+bTNN1MC;%dTNN1MC~D
            bTNN1MC = bTNN1MC + RecIm_M - dTNN1MC;% bTNN1MC~V            
        end
    end
    index=0;
    RecIm_M(RecIm_M<0)=0;
    %imshow(RecIm_M,[]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RING_TV  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if kIter<max_iteration-3 || kIter==max_iteration
        lambda2 =linspace(0.006,0.001,max_iteration); %weight of the TV term in the cost function (see eq(7) in the manuscript)
        alp = 3; % length of the major axis of the ellipse (minor axis is unit-length)
        Center=[256,256];% CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE
        itter=50;        % CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE
 
        RecIm_M_NRTV=RecIm_M;
        RecIm_M_RTV= CircleTVMorteza18_11_05(RecIm_M_NRTV,lambda2(kIter),alp,Center,itter); 
        
        %imshow([RecIm_M_NRTV RecIm_M_RTV],[]); for test if the lam is OK
        %imshow([RecIm_M_RTV-RecIm_M_NRTV],[]);
        if abs(max(RecIm_M_RTV(:)-RecIm_M(:)))==0
            fprintf('Warning:No Change by ring artifact!!!! \n')
        end
        
        RecIm_M=RecIm_M_RTV;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  finding k_d  %%%%%%%%%%%%%%%%%%%%       
        % I should decrease lambda ... ??????????????????????????
        params.im = RecIm_M ;
        Sino_temp = projdd(params);%512*360
        labmda1=1e-9*(0.99)^kIter ;% CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE
        for d=DetStartIndex-2:DetMaxIndex+2     % CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE 
            AdX = Sino_temp(d,:)';                 %1*1024
            AdX_Yd=AdX-G_matrix(d,:)'; % ?????????????????????????????????????????????????????
            Sum__AdX_Yd=sum(AdX_Yd);
            P=NumofView;
            
            NormCase='norm1';
            switch NormCase
                case 'norm2'
%                     labmda2=0.02;
%                     k=(yT_Ad_X - YT_Y)/(YT_Y+labmda2);
                case 'norm1'
                    if (Sum__AdX_Yd > labmda1/2)
                        k = AdX_Yd-labmda1/2;
                    elseif (abs(Sum__AdX_Yd)<=labmda1/2)
                        k = zeros(1,size(Proj,2));
                        index=index+1;
                    else
                        k = AdX_Yd+labmda1/2;
                    end
                    %hist_xAy_yy_yy(d)=xAy_yy_Over_yy;
            end %end switch
            CorrectionProjection(d,:) = k;

        end %end  "for d=68:365"
        
    else
        CorrectionProjection=zeros(size(Proj));
    end %end if
    if kIter==110 || kIter==100 || mod(kIter,max_iteration/4)==0
        %%%%%%%
        figure(kIter); 
        imshow([RecIm_M_NRTV RecIm_M],[0 .6],'border','tight')
        %%%%%%%
        figure(max_iteration+kIter);
        imshow([G_matrix Proj],[1 6.8],'border','tight')
        %title('The correction')
        pause(0.1);
    end
    pause(0.1);
end