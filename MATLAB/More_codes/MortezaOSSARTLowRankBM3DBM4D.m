function imgReconMC = MortezaOSSARTLowRankBM3DBM4D(ProjEnergy_noisy, projTotalWeightMC, params, Niter_OSnum, fSaveName, SaveMidResultsMode,LowRank,BM4D)
%% parameters

reconsize = size(params.imgInitMC,1);          %M 512
          
[~, NumofView, nMC] = size(ProjEnergy_noisy);

Niter = Niter_OSnum(1, :);      %M 25    % vector of main loops numbers
NiterSum = sum(Niter);          %M 25
OSnum = Niter_OSnum(2, :);     %M 20      % vectior of subsets numbers
OSnumIter = [];                           % subsets numbers of each iteration
for i = 1:size(Niter_OSnum,2)
    OSnumIter = [OSnumIter, OSnum(i)*ones(1, Niter(i))]; %M [20 20 ... 20] 50taa
end
imgkMC = params.imgInitMC;    % initial image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramsones = params;
paramsones.im = ones(reconsize);
projones = projdd(paramsones);
paramsones.reconsize = reconsize;
imgRaypixsumMC = zeros(reconsize, reconsize);

imgRaypixsumMC = zeros(reconsize, reconsize, nMC);
for iMC = 1:nMC
    % backprojection
    paramsones.proj = projTotalWeightMC(:, :, iMC).*projones;
    bprojones = bprojdd(paramsones);
    imgRaypixsumMC(:,:,iMC) = bprojones;
end

lambdaNuclear = params.lambdaNuclear; %M 0.025
muNuclear = params.muNuclear;         %M 0.25
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
        
        
        if LowRank==1
        % low rank
        typeflag = 'Thresh';
        nNormThresh = lambdaNuclear/muNuclear;
        imgkMC = TNN1(imgkMC, nNormThresh, typeflag);
        end
             
    end  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          Low rank          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%     if LowRank==1
%     % low rank
%     typeflag = 'Thresh';
%     nNormThresh = lambdaNuclear/muNuclear;
%     imgkMC = TNN1(imgkMC, nNormThresh, typeflag);
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          BM3D          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if (BM4D==0 && mod(kIter,5)==0) %BM3D
        
        sigma=(26-kIter)/10;
        %sigma=1;
        for iii=1:size(ProjEnergy_noisy,3) 
        [~, imgkMC(:,:,iii)] = BM3D([], imgkMC(:,:,iii), sigma, 'np', 1);
        end
        
    elseif (BM4D==1 && mod(kIter,1)==0) %BM4D
        
        [imgkMC, ~] = bm4d(imgkMC, 'Gauss', 0, 'mp', 1, 0);        
        
    end
    
     
% save results
    if (SaveMidResultsMode>= 1 && mod(kIter,5)==0)
        save(strcat(fSaveName,strcat('imgkMC_',num2str(kIter),'-',num2str(iOS))),'imgkMC');
    end       
end
imgReconMC = imgkMC;
