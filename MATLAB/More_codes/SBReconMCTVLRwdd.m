function imgReconMC = SBReconMCTVLRwdd(ProjEnergy_noisy, projTotalWeightMC, params, Niter_OSnum, fSaveName, SaveMidResultsMode)
%          
% SBReconMCTVTNN1.m
%
% TV and low rank based reconstruction for spectral CT,
% solved by Split-Bregman method.
%
% Input variables:
%           'ProjMC'                            spectral sinogram, 3D
%           'projTotalWeightMC'       projection weighting, 3D
%           'params'                            parameters, including initial image, regularization 
%                                                     coefficients,parameters of dictionary, etc.
%           'Niter_OSnum'                   iteration number and subset number
%           'pixelsize'                          real size of each pixel
%           'filepathproj'                     path of system matrix
%           'SaveMidResultsMode'     results saving types: 0-only save the final result; 
%                                                     >1- save the main loop results
%
% Output variables:
%           'imgReconMC'                  result
%
% Yanbo Zhang
% University of Massachusetts Lowell
% yanbozhang007@gmail.com
% 2015-07-13

%% parameters

reconsize = size(params.imgInitMC,1);          %M 512
[~, NumofView, nMC] = size(ProjEnergy_noisy);            %M 640,4

Niter = Niter_OSnum(1, :);      %M 25    % vector of main loops numbers
NiterSum = sum(Niter);          %M 25
OSnum = Niter_OSnum(2, :);     %M 20      % vectior of subsets numbers
OSnumIter = [];                           % subsets numbers of each iteration
for i = 1:size(Niter_OSnum,2)
    OSnumIter = [OSnumIter, OSnum(i)*ones(1, Niter(i))]; %M [20 20 ... 20] 50taa
end

lambdaTV = params.lambdaTV; %M 8*1e-5
muTV = params.muTV;         %M0.1
lambdaNuclear = params.lambdaNuclear; %M 0.025
muNuclear = params.muNuclear;         %M 0.25

imgkMC = params.imgInitMC;    % initial image

%% Reserve memory for the auxillary variables

dxMC = zeros(reconsize,reconsize, nMC);
dyMC = zeros(reconsize,reconsize, nMC);
dTNN1MC = zeros(reconsize,reconsize, nMC);
bxMC = zeros(reconsize,reconsize, nMC);
byMC = zeros(reconsize,reconsize, nMC);
bTNN1MC = zeros(reconsize,reconsize, nMC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramsones = params;
paramsones.im = ones(reconsize);
projones = projdd(paramsones);
paramsones.reconsize = reconsize;
imgRaypixsumMC = zeros(reconsize, reconsize, nMC);
for iMC = 1:nMC
    % backprojection
    paramsones.proj = projTotalWeightMC(:, :, iMC).*projones;
    bprojones = bprojdd(paramsones);
    imgRaypixsumMC(:,:,iMC) = bprojones;
end

%% main loop

params.projMC = ProjEnergy_noisy;
for kIter = 1:NiterSum
    kIter
    OS = OSnumIter(kIter);          % subset number of current iteration
    
    % iteration of each subset
    for iOS = 1:OS
        
        % compute numerator of fidelity part
%         iViews = iOS:OS:NumofView;
%         gradFidMC = OSSARTcoreMC(ProjMC, projTotalWeightMC, imgkMC, iViews, reconsize, pixelsize,  filepathproj);
        indViewOS = iOS:OS:NumofView;
        gradFidMC = - OSSARTCoreMCwdd(params, imgkMC, projTotalWeightMC, indViewOS);
        
        % OS-SART update
        imgkMC = imgkMC - NumofView/length(indViewOS)*gradFidMC./imgRaypixsumMC;
        
        % nonnegtivity
        imgkMC(imgkMC<0) = 0;
        
        % min TV & low rank
        for iMC = 1:nMC
            imgkMC(:,:,iMC) = gsU2TNN1(imgkMC(:,:,iMC), imgkMC(:,:,iMC), dxMC(:,:,iMC), dyMC(:,:,iMC), dTNN1MC(:,:,iMC), bxMC(:,:,iMC), byMC(:,:,iMC), bTNN1MC(:,:,iMC), muTV, muNuclear, reconsize, reconsize);
            [dxMC(:,:,iMC), dyMC(:,:,iMC)] = gsSpace2(imgkMC(:,:,iMC), dxMC(:,:,iMC), dyMC(:,:,iMC), bxMC(:,:,iMC), byMC(:,:,iMC), lambdaTV/muTV, reconsize, reconsize);
            bxMC(:,:,iMC) = bregmanX(dxMC(:,:,iMC), imgkMC(:,:,iMC), bxMC(:,:,iMC), reconsize, reconsize);
            byMC(:,:,iMC) = bregmanY(dyMC(:,:,iMC), imgkMC(:,:,iMC), byMC(:,:,iMC), reconsize, reconsize);
        end
        
        % low rank
        typeflag = 'Thresh';
        nNormThresh = lambdaNuclear/muNuclear;
        dTNN1MC = TNN1(imgkMC+bTNN1MC, nNormThresh, typeflag);
        bTNN1MC = bTNN1MC + imgkMC - dTNN1MC;
         
    end

    % save results
    if (SaveMidResultsMode>= 1 && mod(kIter,5)==0)
        save(strcat(fSaveName,strcat('imgkMC_',num2str(kIter),'-',num2str(iOS))),'imgkMC');
    end
        
end

imgReconMC = imgkMC;
%save(strcat(fSaveName, 'imgReconMC'),'imgReconMC');  % save results


