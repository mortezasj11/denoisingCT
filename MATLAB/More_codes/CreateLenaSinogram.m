clear
clc

Im=imread('lena512.jpg');
Im=imread('Lena256.jpg');
% Im = phantom('Modified Shepp-Logan',256);

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

Im1=zeros(512,512);
Im1(129:512-128,129:512-128)=Im;
imshow(Im1,[])
paramsones.im = Im1;
Proj = projdd(paramsones);
imshow(Proj,[])

LenaProj=Proj;

%corrupting the detector bins by
indexx=randi(200,30,1)+150;

for i=1:length(indexx)
    Plus1OrMinus1=real(exp(sqrt(-1)*randi(2)*pi));
    LenaProj(indexx(i),:)=LenaProj(indexx(i),:)+ Plus1OrMinus1*( (2*randn/100)*LenaProj(indexx(i),:) );
end
imshow(LenaProj,[])

Proj=LenaProj;



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
max_iteration=20;

%% ONLY SART
RecIm_M0=RecIm_M;
params.projMC = Proj;
% iteration of each subset
for SART=1:25
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
figure,imshow(RecIm_M0,[])
%%
RecIm_M=zeros(reconsize,reconsize);
index=0;
for kIter = 1:max_iteration

    G_matrix=G_matrix+ 1*(G_matrix.*CorrectionProjection);
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
    if kIter~=max_iteration
        lam = 0.002; %weight of the TV term in the cost function (see eq(7) in the manuscript)
        alp = 3; % length of the major axis of the ellipse (minor axis is unit-length)
        Center=[256,256];
        itter=50;   
        RecIm_M_NRTV=RecIm_M;
        RecIm_M_RTV= CircleTVMorteza18_11_05(RecIm_M_NRTV,lam,alp,Center,itter); 
        
        %imshow([RecIm_M RecIm_M_RTV],[]);% for test if the lam is OK
        %imshow([RecIm_M_RTV-RecIm_M],[]);
        
        RecIm_M=RecIm_M_RTV;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  finding k_d  %%%%%%%%%%%%%%%%%%%%       
        % I should decrease lambda ... ??????????????????????????
        params.im = RecIm_M ;
        Sino_temp = projdd(params);%448*1024
        labmda1=0.001 ;
        for d=163:345  
            AdX = Sino_temp(d,:)';                 %1*1024
            yT_Ad_X = G_matrix(d,:)*AdX;       %should be a number
            YT_Y = G_matrix(d,:)*G_matrix(d,:)';
            if YT_Y==0
                xAy_yy_Over_yy=0;
            else
                xAy_yy_Over_yy = (yT_Ad_X - YT_Y)/YT_Y;
            end

            NormCase='norm1';

            switch NormCase
                case 'norm2'
                    labmda2=0.04;
                    k=(yT_Ad_X - YT_Y)/(YT_Y+labmda2);

                case 'norm1'

                    if (xAy_yy_Over_yy > labmda1)
                        k = xAy_yy_Over_yy-labmda1;
                    elseif (abs(xAy_yy_Over_yy)<=labmda1)
                        k = 0;
                        index=index+1;
                    else
                        k = xAy_yy_Over_yy+labmda1;
                    end
                    %hist_xAy_yy_yy(d)=xAy_yy_Over_yy;
            end %end switch
            CorrectionProjection(d,:) = k*ones(1,size(Proj,2));

        end %end  "for d=68:365"

    end %end if
    if kIter==110 || kIter==100 || mod(kIter,1)==0
        %%%%%%%
        figure(kIter); 
        imshow([RecIm_M RecIm_M0],[])
        %%%%%%%
        figure(max_iteration+kIter);
        %G_matrix(158,:)=1000;G_matrix(350,:)=1000;
        imshow([G_matrix Proj],[])
        title('The correction')
        pause(0.1);
    end

end