clear all
close all
clc

%% PARAMETERS 

params.Dsource2centr = 132;
params.Dsource2detec = 180;
params.NumofView     = 640;
params.NumofBin      = 512;       % number of detector bins
params.pixelsize     = 0.075;
params.binsize       = 0.1;
params.binshift      = 0;         % detector shift, mm

params.reconsize     = 512;

%% materials

load MiuofH2O
load MiuofBLOOD_Whole_ICRU44
load MiuofI

% Iodine proportion by weight, assuming that volume of 
% water doesn't change when Iodine solutes in water

%% compute the spectral projection

load Spectrum20131219interp
Spectrumtemp = Spectrum20131219interp;
AttenuMode = 7;
%pixelsize = 0.075;       % 0.02mm/pix

maxEnergy = max(Spectrumtemp(:,1));
% Iodine k-edge = 33.17keV
energyBin = [16, 22, 25, 28, 31, 34, 37,41, maxEnergy];
nEnergy = length(energyBin) - 1;
Spectrum = Spectrumtemp(:,2)/sum(Spectrumtemp(:,2));  % normalization
clear Spectrumtemp Spectrum20131219interp
%% generate projection of mouse phantom and save

load('MousePhan_slice875_kev1_kev50.mat')

minEn = energyBin(1);
maxEn = energyBin(end);

projEn=zeros(params.reconsize,params.NumofView,50);

jjj=0;
for ii=minEn:maxEn
    jjj=jjj+1;
params.im = MousePhan_slice875_kev1_kev50(:,:,ii);
projEn(:,:,ii) = projdd(params);
end
%save( 'projEn8.mat','projEn')

%%%%%%%% generate the projection of iodine in the blood of mouse %%%%%%%%%
load('bwBlood.mat')
params.im = bwBlood;
projbwBlood = projdd(params)/10;
% save('projbwBlood.mat','projbwBlood')

%save('params_energyBin_Spectrum.mat','params','energyBin','Spectrum')

%%

clear all
load('projbwBlood');
load('projEn8');
load('params_energyBin_Spectrum.mat');
load MiuofI
%% compute the projection at each engery channel

AttenuMode=7;
nEnergy = length(energyBin) - 1;
densityIwater = 1;   ratioI = 0.012;  
[row, col] = size(projbwBlood);
ProjEnergy = zeros(row, col, nEnergy);
ProjEnergyBlankRatio = zeros(row, col, nEnergy);
for iEnBin = 1:nEnergy
    iEnBin
    for iEn = energyBin(iEnBin):energyBin(iEnBin+1)-1        
        temp = MiuofI(iEn,AttenuMode)*densityIwater*ratioI*projbwBlood + projEn(:,:,iEn);
        Ptmp = Spectrum(iEn)*exp(-temp);    % according to the spectrum ratio
        ProjEnergy(:,:, iEnBin) = ProjEnergy(:,:, iEnBin) + Ptmp;      
    end
    ProjEnergyBlankRatio(:,:, iEnBin) = sum(Spectrum(energyBin(iEnBin):energyBin(iEnBin+1)-1))*ones(row, col);
end
ProjEnergy = -log(ProjEnergy./ProjEnergyBlankRatio);

% save('ProjEnergy.mat','ProjEnergy','ProjEnergyBlankRatio')
%% Adding Noise   "ProjEnergy_noisy"

%%%% here we work with the "ProjEnergy" and "spectrum" and "energyBin"
channel=zeros(nEnergy,1);
for ii=1:nEnergy
    for jj = energyBin(ii):energyBin(ii+1)-1  
        channel(ii)=channel(ii)+Spectrum(jj);
    end
end
nSndEachCHL_percent = channel / sum(Spectrum(15:49));

nSendTotal=4e4;
nSend_CHNL = nSendTotal*nSndEachCHL_percent;


nReceive          =    zeros(512,640,4);
nReceive_hat      =    zeros(512,640,4);
ProjEnergy_noisy  =    zeros(512,640,4);
for ii=1:length(nSndEachCHL_percent) 
   nReceive(:,:,ii) = (  nSend_CHNL(ii)  ) * exp( -ProjEnergy(:,:,ii)  );
   nReceive_hat(:,:,ii) =   poissrnd(  nReceive(:,:,ii)  );
   ProjEnergy_noisy(:,:,ii) = - log(  nReceive_hat(:,:,ii)  /  nSend_CHNL(ii)  );
end   
% save('ProjEnergy_noisy4e4.mat','ProjEnergy','ProjEnergy_noisy','params');
% clear nReceive nReceive_hat ii nSndEachCHL_percent nSend_CHNL jj channel

save('ProjEnergy_noisy5e3.mat','ProjEnergy','ProjEnergy_noisy','params');
clear nReceive nReceive_hat ii nSndEachCHL_percent nSend_CHNL jj channel

%%
clc;clear all;
load 'ProjEnergy_noisy4e4.mat'

%Low Rank
 TNN1Thresh = 0.1;  % threshold of SVT
 params.muNuclear = 0.25;%0.1;
 params.lambdaNuclear = params.muNuclear*TNN1Thresh;
 
reconsize = 512; % size of image
nMC = size(ProjEnergy, 3);%number of channel
imgInitMC = zeros(reconsize,reconsize, nMC);   % initalized image

% Niter_OSnum is a two rows matrix, which specifies the numbers of main loops
% and corresponding subsets respectively of OS
Niter_OSnum = [25; 20];% [50; 40];%[10 5 3 2; 10 5 2 1];
NiterSum = sum(Niter_OSnum(1, :));        % number of total main loops

params.imgInitMC = imgInitMC;

%% Setting for results saving

codepath=pwd;
% saving directory of results
filesavepath0 = [codepath, '\Results4e4\BM3D\'];
t = clock;        % current time
tstr = strcat('_',num2str(t(1)),'-',num2str(t(2)),'-',num2str(t(3)),'-',num2str(t(4)),'-',num2str(t(5)),'-',num2str(floor(t(6))),'NoLR-NoBM34D');
fSaveName = strcat('Case', tstr);     % file name
fSaveName = strcat(filesavepath0, fSaveName,'\');
mkdir(fSaveName);       % create folder

% results saving types: 0-only save the final result; 1- save the main loop results; 2-save all results
SaveMidResultsMode =1;
%% iterative reconstruction

projTotalWeightMC = ones(size(ProjEnergy_noisy));

LowRank=0;
BM4D=2;%0:BM3D     , 1:BM4D   

imgReconMC = MortezaOSSARTLowRankBM3DBM4D(ProjEnergy_noisy, projTotalWeightMC, params, Niter_OSnum, fSaveName, SaveMidResultsMode,LowRank,BM4D);

