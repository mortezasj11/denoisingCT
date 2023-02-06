clear all%117 118 119
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

load MiuofH2O%120*8
load MiuofBLOOD_Whole_ICRU44%120*8
load MiuofI%120*8

% Iodine proportion by weight, assuming that volume of 
% water doesn't change when Iodine solutes in water

%% compute the spectral projection

load Spectrum20131219interp%50*2
Spectrumtemp = Spectrum20131219interp;
AttenuMode = 7;
%pixelsize = 0.075;       % 0.02mm/pix
maxEnergy = max(Spectrumtemp(:,1));
% Iodine k-edge = 33.17keV
energyBin = [15, 22, 27, 34, 50];%energyBin = [16, 22, 25, 28, 31, 34, 37,41, maxEnergy];energyBin = [16,31,50];
nEnergy = length(energyBin) - 1;
Spectrum = Spectrumtemp(:,2)/sum(Spectrumtemp(:,2));  % normalization
clear Spectrumtemp Spectrum20131219interp
%% generate projection of mouse phantom and save

load('MousePhan_slice875_kev1_kev50.mat')%512*512*50

minEn = energyBin(1);
maxEn = energyBin(end);

projEn=zeros(params.reconsize,params.NumofView,50);%512*640*50

jjj=0;
for ii=minEn:maxEn
    jjj=jjj+1;
    params.im = MousePhan_slice875_kev1_kev50(:,:,ii);
    projEn(:,:,ii) = projdd(params);
end
%save( 'projEn.mat','projEn')

%%%%%%%% generate the projection of iodine in the blood of mouse %%%%%%%%%
load('bwBlood.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%   MORTEZA shatranji   %%%%%%%%%%%%%%%%%%%%%%%%%%
%             N=7;
%             sigmaaa=3;
%             % N is grid size, sigma speaks for itself
%             [xxx ,yyy]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
%             t1=exp(-xxx.^2/(2*sigmaaa^2)-yyy.^2/(2*sigmaaa^2));
%             t1=repmat(t1,100);
%             t2=t1(1:512,1:512);
%             bwBlood=bwBlood.*t2;
%             imshow(bwBlood,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%   MORTEZA shatranji   %%%%%%%%%%%%%%%%%%%%%%%%%%

params.im = bwBlood;
projbwBlood = projdd(params)/10;
% save('projbwBlood.mat','projbwBlood')

%save('params_energyBin_Spectrum.mat','params','energyBin','Spectrum')

%%
clear all
load('projbwBlood');
load('projEn');
%load('params_energyBin_Spectrum.mat');
load('params_energyBin_Spectrum8.mat');
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

nSendTotal=5e3;
nSend_CHNL = nSendTotal*nSndEachCHL_percent;


nReceive          =    zeros(512,640,nEnergy);
nReceive_hat      =    zeros(512,640,nEnergy);
ProjEnergy_noisy  =    zeros(512,640,nEnergy);
for ii=1:length(nSndEachCHL_percent) 
   nReceive(:,:,ii) = (  nSend_CHNL(ii)  ) * exp( -ProjEnergy(:,:,ii)  );
   nReceive_hat(:,:,ii) =   poissrnd(  nReceive(:,:,ii)  );
   ProjEnergy_noisy(:,:,ii) = - log(  nReceive_hat(:,:,ii)  /  nSend_CHNL(ii)  );
end   
% save('ProjEnergy_noisy1e4.mat','ProjEnergy','ProjEnergy_noisy','params');
% clear nReceive nReceive_hat ii nSndEachCHL_percent nSend_CHNL jj channel
ProjEnergy_noisy_Simple7e3_8=ProjEnergy_noisy;
ProjEnergy_Orig_Simple7e3_8=ProjEnergy;
save('ProjEnergy_noisyOrig_Simple7e3_8.mat','ProjEnergy_noisy_Simple7e3_8','ProjEnergy_Orig_Simple7e3_8','params');
clear nReceive nReceive_hat ii nSndEachCHL_percent nSend_CHNL jj channel